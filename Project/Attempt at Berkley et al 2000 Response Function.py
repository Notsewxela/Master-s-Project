import numpy as np
from scipy.integrate import simps as integrate
import pylab
from math import pi

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
    xray_heating_term = x_ray * const2 * (1-A) * L_x * h_x * (1-np.exp(-(1-(r_isco(a)/r)**(1/2))**(4/5))) / (m**2 * (h_x**2+r**2)**(3/2))
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

def D(r, m, h_x, phi, inc_deg):
    '''Calculates the extra distance light travels to a point on the disk and towards the observer'''
    #This function intends to receive many radii and a single phi 
    inc = inc_deg*pi/180
    return G*M_sun*m/c**2 * (np.sqrt(r**2+h_x**2) + h_x*np.cos(inc) - r*np.cos(phi)*np.sin(inc))

def area(r, m, h_x, inc_deg, tau):
    '''Calculates the effective area in the response integral for a single radius/radii and single tau'''
    tau_l = D(r, m, h_x, 0, inc_deg)/c #creates arrays of tau_l for each r inputted
    tau_t = D(r, m, h_x, pi, inc_deg)/c #creates arrays of tau_t for each r inputted
    
    '''Here we convert to a complex array so np.sqrt doesn't error. This way allows us to take a whole array
    and calculate area for more than one value. We set to complex then sqrt then remove imaginary parts which, for example
    sqrt(-9) = 3j. Then real(3j) = 0 so the area is 0!'''
    taus = np.zeros(len(r), dtype = np.complex)+tau #makes a complex data type array of value tau
    
    taus[(taus < tau_l) & (taus > tau_t)] = 0j #forcing the area to 0 where appropriate
    
    return np.real(2*G*M_sun*m/c**3 * r/np.sqrt((taus-tau_l)*(tau_t-taus)))
    
'''if tau_l < tau < tau_t:
        return 2*G*M_sun*m/c**3 * r/np.sqrt((tau-tau_l)*(tau_t-tau))
    else:
        return 0'''

def response(r, m, mdot_E, a, h_x, L_x, A, inc_deg, freqs, taus):
    '''Calculates responses for our inclined disk'''
    T1 = T(r, m, mdot_E, a, h_x, L_x, A, x_ray = 0)
    T2 = T(r, m, mdot_E, a, h_x, L_x, A, x_ray = 1)
    
    '''Creates a 2x2 array of dimension (frequencies, radii)'''
    B2 = B_nu(T2, freqs)
    B1 = B_nu(T1, freqs)
    deltaB = B2-B1
    
    '''Integrates along the frequency axis to get Wm^-2 at each radius'''
    deltaLs = 4*pi*integrate(deltaB, x = freqs, axis = 1)
    responses = np.zeros(len(taus))
    
    integrand = np.zeros((len(taus), len(r)))
    
    for i in range(len(taus)):
        integrand[i] = deltaLs*area(r, m, h_x, inc_deg, taus[i])
    '''for i in range(len(taus)):
        for j in range(len(r)):
            integrand[i][j] = deltaLs[j]*area(r[j], m, h_x, inc_deg, taus[i])'''
    responses = integrate(integrand, x = r, axis = 1)
    
    #return np.convolve(responses*G*M_sun*m/c**2, v = np.ones(1), mode = "same")
    return responses*G*M_sun*m/c**2
    
freq_data_points = 1000 #increase this to get the sum of the transfer function data points be
                            #closer to the area under the disk spectrum
                            
inclinations = [10, 20, 30, 40, 50, 60, 70, 80]
#inclination = 10 #0 <= degrees < 90

mass_BH = 5e7 # mass in solar unit
mdot_eddington = 0.05 #accretion rate in eddington units
spin = 0 #dimensionless spin parameter
r_in = r_isco(spin)+1e-15 #gravitational units
r_out = 1e3 #outer edge of disk
'''We have to add a small number to the innermost radius as other wise our T(r=r_isco) generates a 0 value
which causes our _inside_of_integral function to throw an error due to dividing by zero in the
exponential function. This change of r_min has no noticable effect on results.'''

freq_low = 10**(5) #5
freq_high = 10**(17.5)# #17.5
frequencies = np.geomspace(freq_low, freq_high, num = freq_data_points)

radii_data_points = 10000
radii = np.geomspace(r_in+1e-15, r_out, num = radii_data_points) #in gravitational units
#extra data point is so we can do r[i+1]-r[i] for all the data point we want to test. it is only for this use!

l_X = 0.08 #x ray luminosity as fraction of L_Bol
L_X = l_X*1e12*L_sun #here we take the total xray emission as a fraction of L_Bol taken to be ~1e12 suns
h_X = 3 * radii[0] #x ray source height above disk in gravitational units
albedo = 0.3 #albedo

tau_data_points = 10000



pylab.figure(figsize = (10, 8))
for inclination in inclinations:
    tau_min = D(radii[0], mass_BH, h_X, 0, inclination)/c
    tau_max = D(radii[-1], mass_BH, h_X, pi, inclination)/c #in seconds
    taus = np.linspace(tau_min, tau_max, num = tau_data_points)
    responses = response(radii, mass_BH, mdot_eddington, spin, h_X, L_X, albedo, inclination, frequencies, taus)
    pylab.plot(taus[responses < 8e24], responses[responses < 8e24], label = str(inclination) + "°")


pylab.title("Response Functions", fontsize = 16)
pylab.xlabel("Lag (s)", fontsize = 14)
pylab.ylabel("Response Function (Ws$^{-1}$)", fontsize = 14)

pylab.xscale("log")
#pylab.yscale("log")
pylab.xticks(fontsize=12)
pylab.yticks(fontsize=12)
pylab.legend(fontsize = 14)
#pylab.savefig("Response-Function.pdf")
pylab.show()

#print("Total observed extra luminosity from x-ray heating " + "{:e}".format(integrate(transfer_unnice, x=taus_unnice) / L_sun) + " L_☉")











