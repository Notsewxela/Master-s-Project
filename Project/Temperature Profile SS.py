import numpy as np
from scipy.integrate import trapz as integrate
import pylab
from math import pi

m_p = 1.6726231e-27
c = 2.99792458e8
sigma_SB = 5.67e-8
sigma_T = 6.652e-29
G = 6.67259e-11
M_sun = 1.99e30
h = 1.05457266e-34*2*pi
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
    optical_depth = (1-np.exp(-(1-(r_isco(a)/r)**(1/2))**(4/5)))#set to 1 to ignore optical depth thing
    xray_heating_term = x_ray * const2 * (1-A) * L_x * h_x * optical_depth / (m**2 * (h_x**2+r**2)**(3/2))
    return (T_old_forthpower + xray_heating_term)**(1/4)

def r_isco(a):
    '''Given the dimensionless spin parameter a from [-1, 1] calculates the
    innermost stable circular orbit in gravitational radii. This function cleverly
    uses special units to make equations simpler.'''
    if abs(a) > 1:
        raise Exception("a must lie on the closed interval [-1, 1]")
    if a == 0:
        #non spinning black hole
        return 6
    Z1 = 1 + (1-a**2)**(1/3)*((1+a)**(1/3)+(1-a)**(1/3))
    Z2 = (3*a**2 + Z1**2)**(1/2)
    Z3 = ((3-Z1)*(3+Z1+2*Z2))**(1/2) 
    if a > 0:
        #prograde black hole
        return 3 + Z2 - Z3
    else:
        #retrograde black hole
        return 3 + Z2 + Z3

def eddington_efficiency(a):
    '''Calculates the eddington efficiency dependent on the inner accretion disk
    radius which is assumed to be at the ISCO. The ISCO depends on the dimensionless
    spin parameter a.'''
    return 1-(1-2/(3*r_isco(a)))**(1/2)

def B_nu(Ts, nus):
    '''Calculates the blackbody function WHz^-1m^-2Sr^-1 and removes Nan values.
    T and nu can both be inputted as arrays!'''
    Ts_col = np.vstack(Ts)
    
    B = 2*h*nus**3/c**2 / (np.exp(h * nus / (k_B * Ts_col)) - 1)
    return np.nan_to_num(B)

inclination = 0 #degrees
mass_BH = 5e7 # mass in solar unit
mdot_eddington = 0.05 #accretion rate in eddington units
a = 1 #dimensionless spin parameter
r_in = r_isco(a)+1e-15
r_out = 1e3 #outer edge of disk

radii_data_points = 1000
radii = np.geomspace(r_in, r_out, num = radii_data_points) #in gravitational units
'''We have to add a small number to the innermost radius as other wise our T(r=r_isco) generates a 0 value
which causes our _inside_of_integral function to throw an error due to dividing by zero in the
exponential function. This change of r_min has no noticable effect on results.'''

l_x = 0.08 #x ray luminosity as fraction of L_Bol
L_x = l_x*1e12*L_sun #here we take the total xray emission as a fraction of L_Bol taken to be ~1e12
h_x = 3 * radii[0] # x ray source height above disk in gravitational units
A = 0.3 #albedo

freq_data_points = 1000
freq_low =  1e5#1.3*10**(15) #
freq_high = 10**(17.5)#1.9*10**(15) #
frequencies = np.geomspace(freq_low, freq_high, num = freq_data_points)

pylab.figure(figsize = (10, 8))


temperatures = T(radii, mass_BH, mdot_eddington, a, h_x, L_x, A, 0)
temperatures_new = T(radii, mass_BH, mdot_eddington, a, h_x, L_x, A, 1)



pylab.title("Disk Temperature Profiles", fontsize = 16)
pylab.xlabel("Radius (R$_g$)", fontsize = 14)
pylab.ylabel("Temperature (K)", fontsize = 14)
pylab.plot(radii, temperatures, label = "SS Temperature Profile")
pylab.plot(radii, temperatures_new, label = "Illuminated Temperature Profile")
pylab.plot(radii, temperatures_new-temperatures, ls="--", color = "gold", label = "Difference")
pylab.xscale("log")
pylab.xticks(fontsize=12)
pylab.yticks(fontsize=12)
pylab.legend(fontsize = 14)
pylab.savefig("Temperature-Profiles.pdf")
pylab.show()

bb1 = 2*pi*integrate(B_nu(temperatures, frequencies), x = frequencies, axis = 1)
bb2 = 2*pi*integrate(B_nu(temperatures_new, frequencies), x = frequencies, axis = 1)
bb3 = 2*pi*integrate(B_nu(temperatures_new, frequencies)-B_nu(temperatures, frequencies), x = frequencies, axis = 1)

pylab.figure(figsize = (10, 8))
pylab.title("Surface Brightness Profiles", fontsize = 16)
pylab.xlabel("Radius (R$_g$)", fontsize = 14)
pylab.ylabel("Surface Brightness (Wm$^{-2}$)", fontsize = 14)
pylab.plot(radii, bb1, label = "SS Surface Brightness")
pylab.plot(radii, bb2, label = "Illuminated Surface Brightness")
pylab.plot(radii, bb3, ls="--", color = "gold", label = "Difference")
pylab.xscale("log")
pylab.xticks(fontsize=12)
pylab.yticks(fontsize=12)
pylab.legend(fontsize = 14)
pylab.savefig("Surface-Brightness-Profiles.pdf")
pylab.show()

