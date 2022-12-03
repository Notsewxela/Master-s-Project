import numpy as np
from scipy.integrate import simps as integrate
import pylab
from math import pi

#define constants
m_p = 1.6726231e-27
c = 2.99792458e8
sigma_SB = 5.67e-8
sigma_T = 6.652e-29
G = 6.67259e-11
M_sun = 1.99e30
h = 1.05457266e-34 * 2*pi #reduced planck constant times 2pi
k_B = 1.380658e-23
L_sun = 3.828e26

def T(r, m, mdot_E, a, h_x, L_x, A, x_ray = 0):
    '''Calculates the temperature of the accretion disk at a radius r (in
    gravitational units) given the black hole mass (in solar masses) and the
    accretion rate (in eddington units), height of x ray source above disk
    (in gravitational units), the xray luminosity (in watts) and the albedo.
    The x-ray parameter is 1 for an illuminated disk, 0 otherwise.'''
    eta = eddington_efficiency(a)
    const = 3*m_p*c**5/(2*eta*sigma_T*sigma_SB*G*M_sun)
    extra_r_factor = 1-(r_isco(a)/r)**(1/2)
    T_old_forthpower = const * (mdot_E * extra_r_factor / (r**3*m))
    const2 = c**4/(4*pi*sigma_SB*G**2*M_sun**2)
    chi_0 = 1
    optical_depth = (1-np.exp(-chi_0*(1-(r_isco(a)/r)**(1/2))**(4/5)))
    xray_heating_term = x_ray * const2 * (1-A) * L_x * h_x * optical_depth / (m**2 * (h_x**2+r**2)**(3/2))
    return (T_old_forthpower + xray_heating_term)**(1/4)

def B_nu(Ts, nus):
    '''Calculates the blackbody function WHz^-1m^-2Sr^-1 and removes Nan values.
    T and nu can both be inputted as arrays!'''
    Ts_col = np.vstack(Ts)
    
    B = 2*h*nus**3/c**2 / (np.exp(h * nus / (k_B * Ts_col)) - 1)
    return np.nan_to_num(B)

def r_isco(a):
    '''Given the dimensionless spin parameter a from [-1, 1] calculates the
    innermost stable circular orbit in gravitational radii. This function cleverly
    uses special units to make equations simpler.'''
    if abs(a) > 1:
        raise Exception("a must lie on the closed interval [-1, 1]")
    Z1 = 1 + (1-a**2)**(1/3)*((1+a)**(1/3)+(1-a)**(1/3))
    Z2 = (3*a**2 + Z1**2)**(1/2)
    Z3 = np.sign(a)*((3-Z1)*(3+Z1+2*Z2))**(1/2) 
    return 3 + Z2 - Z3

def eddington_efficiency(a):
    '''Calculates the eddington efficiency dependent on the inner accretion disk
    radius which is assumed to be at the ISCO. The ISCO depends on the dimensionless
    spin parameter a.'''
    return 1-(1-2/(3*r_isco(a)))**(1/2)

def D(r, m, h_x):
    '''Calculates the extra distance light travels to a point on the disk and towards the observer'''
    return G*M_sun*m/c**2 * (np.sqrt(r**2+h_x**2) + h_x)

        
def response(r, m, mdot_E, a, h_x, L_x, A, freqs):
    '''Calculates the response function of our non-inclined disk.'''    
    T1 = T(r, m, mdot_E, a, h_x, L_x, A, x_ray = 0)
    T2 = T(r, m, mdot_E, a, h_x, L_x, A, x_ray = 1)
    
    '''Creates a 2x2 array of dimension (frequencies, radii)'''
    B2 = B_nu(T2, freqs)
    B1 = B_nu(T1, freqs)
    deltaB = B2-B1
    
    '''Integrates along the frequency axis to get Wm^-2Sr^-1 at each radius'''
    deltaLs = integrate(deltaB, x = freqs, axis = 1)
    
    '''Multiplies by all the other bits to get the correct final answer in Ws^-1! Hooray!'''
    deltaLs *= 4*pi * G*M_sun*m/c * np.sqrt(r**2 + h_x**2) * 2*pi
    return deltaLs #our response function

'''Main block of codes starts here:'''
mass_BH = 5e7 # mass in solar unit
mdot_eddington = 0.05 #accrestion rate in eddington units
spin = 1 #dimensionless spin parameter
r_in = r_isco(spin)+1e-15 #gravitational units
r_out = 10**(3) #outer edge of disk
'''We have to add a small number to the innermost radius as other wise our T(r=r_isco) generates a 0 value
which causes our _inside_of_integral function to throw an error due to dividing by zero in the
exponential function inside of the blackbody function. This change of r_min has no noticable effect on results.'''

radii_data_points = 1000
radii = np.geomspace(r_in+1e-15, r_out, num = radii_data_points) #in gravitational units

l_X = 0.08 #x ray luminosity as fraction of L_Bol
L_X = l_X*1e12*L_sun #here we take the total xray emission as a fraction of L_Bol taken to be ~1e12 suns
h_X = 3 #+ radii[0] #x ray source height above disk in gravitational units
albedo = 0.3 #albedo

freq_data_points = 1000
freq_low =  1e5 #Hz
freq_high = 10**(17.5)
frequencies = np.geomspace(freq_low, freq_high, num = freq_data_points)

responses = response(radii, mass_BH, mdot_eddington, spin, h_X, L_X, albedo, frequencies)
taus = D(radii, mass_BH, h_X) / c

print("Total extra luminosity from disk response spin " + str(spin) + ": {:e}".format(integrate(responses, x = taus) / L_sun) + " L_☉")
pylab.figure(figsize = (10,8))
pylab.title("Response Functions", fontsize = 16)
pylab.plot(taus, responses, label = "Response")
pylab.xlabel("$\\tau$ (s)", fontsize = 14)
pylab.ylabel("$\\frac{\\partial[\\Delta L(r)]}{\\partial\\tau}$ (Ws$^{-1}$)", fontsize = 14)
pylab.xscale("log")
#pylab.yscale("log")
pylab.xticks(fontsize = 12)
pylab.yticks(fontsize = 12)
pylab.legend(fontsize = 12)
pylab.savefig("response.pdf", bbox_inches = 'tight')
pylab.show()