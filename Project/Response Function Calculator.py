import numpy as np
from scipy.integrate import simps as integrate
import pylab
from math import pi
import scipy

m_p = 1.6726231e-27
c = 2.99792458e8
sigma_SB = 5.67e-8
sigma_T = 6.652e-29
G = 6.67259e-11
M_sun = 1.99e30
h = 1.05457266e-34*2*pi
k_B = 1.380658e-23
L_sun = 3.828e26

def T(r, m, mdot_E, a, h_x, L_x, A):
    '''Calculates the temperature of the accretion disk at a radius r (in
    gravitational units) given the black hole mass (in solar masses) and the
    accretion rate (in eddington units).'''
    #THE LAST 3 PARAMETERS ARE IRRELEVANT AND ARE TO SAVE REDEFINITION OF L FOR DIFFERENT T
    #r is radius in gravitational units
    #m is the BH mass in solar masses
    #mdot_E is the accretion rate in eddington units
    eta = eddington_efficiency(a) #efficiency
    const = (3*m_p*c**5/(2*eta*sigma_T*sigma_SB*G*M_sun))**0.25
    extra_r_factor = 1-(r_isco(a)/r)**(1/2)
    return const * (mdot_E * extra_r_factor / (r**3*m))**0.25

def T_new(r, m, mdot_E, a, h_x, L_x, A):
    '''Calculates the new effective temperature of the disk after it has absorbed
    some x-rays from the xray source. Same inputs as T with the addition of
    height of x ray source above disk (in gravitational units), the xray luminosity
    (in watts) and the albedo.'''
    eta = eddington_efficiency(a)
    const = (3*m_p*c**5/(2*eta*sigma_T*sigma_SB*G*M_sun))
    extra_r_factor = 1-(r_isco(a)/r)**(1/2)
    T_old_forthpower = const * (mdot_E * extra_r_factor / (r**3*m))
    const2 = c**4/(4*pi*sigma_SB*G**2*M_sun**2)
    xray_heating_term = const2 * (1-A) * L_x * h_x / (m**2 * (h_x**2+r**2)**(3/2))
    return (T_old_forthpower + xray_heating_term)**(1/4)

def B_nu(T, nu):
    '''Calculates the blackbody function WHz^-1m^-2Sr^-1'''
    return 2*h*nu**3/c**2 / (np.exp(h * nu / (k_B * T)) - 1)

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

def D(r, m, h_x, phi, inc_deg):
    '''Calculates the extra distance light travels to a point on the disk and towards the observer'''
    #This function intends to receive many radii and a single phi 
    inc = inc_deg*pi/180
    return G*M_sun*m/c**2 * (np.sqrt(r**2+h_x**2) + h_x*np.cos(inc) - r*np.cos(phi)*np.sin(inc))

def area(radii, m, tau, tau_leading, tau_trailing):
    '''Calculates a projected area for a range of radii for a delay tau'''
    if tau_leading < tau < tau_trailing:
        #makes sure the lag is within the bounds of the min and max values
        return G*M_sun*m/c**3 * 2*radii/np.sqrt((tau-tau_leading)*(tau_trailing-tau))
    return 0
    
def Psi_nu(nu, radii, m, mdot_E, a, h_x, L_x, A, inc_deg, time):
    '''Calculates our response function at a frequency'''
    for r in radii:
        tau_l = D(r, m, h_x,)
    return 
    

freq_data_points = 1000 #increase this to get the sum of the transfer function data points be
                            #closer to the artea under the disk spectrum

inclination = 0 #degrees
mass_BH = 5e7 # mass in solar unit
mdot_eddington = 0.05 #accretion rate in eddington units
a = 1 #dimensionless spin parameter
r_in = r_isco(a)+1e-15
r_out = 1e5 #outer edge of disk

freq_low = 10**(15) #5
freq_high = 10**(16) #17.5
frequencies = np.geomspace(freq_low, freq_high, num = freq_data_points)

radii_data_points = freq_data_points
radii = np.geomspace(r_in+1e-15, r_out, num = radii_data_points) #in gravitational units
'''We have to add a small number to the innermost radius as other wise our T(r=r_isco) generates a 0 value
which causes our _inside_of_integral function to throw an error due to dividing by zero in the
exponential function. This change of r_min has no noticable effect on results.'''
phis = np.linspace(0, pi, radii_data_points) #as the disk is symmetric we only need to do 0->pi and double the result in our area function


taus = np.zeros(radii_data_points)

l_x = 0.08 #x ray luminosity as fraction of L_Bol
L_x = l_x*1e12*L_sun #here we take the total xray emission as a fraction of L_Bol taken to be ~1e12
h_x = 3 * radii[0] # x ray source height above disk in gravitational units
A = 0.3 #albedo

tau_leading = 










