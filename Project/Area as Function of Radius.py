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
    if T < 1e-5:
        #if nearly 0 pretend it's 0 to stop dividing by 0
        return 0
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

freq_data_points = 1000 #increase this to get the sum of the transfer function data points be
                            #closer to the artea under the disk spectrum

inclination = 0 #degrees
mass_BH = 5e7 # mass in solar unit
mdot_eddington = 0.05 #accretion rate in eddington units
a = 1 #dimensionless spin parameter
r_in = r_isco(a)+1e-15
r_out = 1e2 #outer edge of disk

freq_low = 10**(5) #5
freq_high = 10**(17.5) #17.5
frequencies = np.geomspace(freq_low, freq_high, num = freq_data_points)

radii_data_points = 1000
radii = np.geomspace(r_in+1e-15, r_out, num = radii_data_points + 1) #in gravitational units
'''We have to add a small number to the innermost radius as other wise our T(r=r_isco) generates a 0 value
which causes our _inside_of_integral function to throw an error due to dividing by zero in the
exponential function. This change of r_min has no noticable effect on results.'''
phi_data_points = 100
phis = np.linspace(0, pi, phi_data_points) #as the disk is symmetric we only need to do 0->pi and double the result


taus = np.zeros(radii_data_points*phi_data_points) #stored smartly xd
transfer = np.zeros(radii_data_points*phi_data_points)

l_x = 0.08 #x ray luminosity as fraction of L_Bol
L_x = l_x*1e12*L_sun #here we take the total xray emission as a fraction of L_Bol taken to be ~1e12
h_x = 3 * radii[0] #x ray source height above disk in gravitational units
A = 0.3 #albedo

area = np.zeros(radii_data_points)
B1 = np.zeros(radii_data_points)
B2 = np.zeros(radii_data_points)


for i in range(radii_data_points):
    T2 = T_new(radii[i], mass_BH, mdot_eddington, a, h_x, L_x, A)
    T1 = T(radii[i], mass_BH, mdot_eddington, a, h_x, L_x, A)
    B2[i] = integrate(B_nu(T2, frequencies))
    B1[i] = integrate(B_nu(T1, frequencies))
    #B2[i] = B_nu(T2, 1e15)
    #B1[i] = B_nu(T1, 1e15)
    area[i] = 2*pi*radii[i]

radii = radii[:-1]

pylab.figure(figsize = (10, 8))
pylab.title("Blackbody*Area", fontsize = 16)
pylab.xlabel("Radius (Rg)", fontsize = 14)
pylab.ylabel("Blackbody (W)", fontsize = 14)
pylab.plot(radii, B2*area, label="B2")
pylab.plot(radii, B1*area, label="B")
pylab.plot(radii, (B2-B1)*area, label="BDIFF")
pylab.xscale("log")
#pylab.yscale("log")
pylab.xticks(fontsize=12)
pylab.yticks(fontsize=12)
pylab.legend(fontsize = 14)
pylab.show()


#print("Total extra luminosity from x-ray heating " + "{:e}".format(integrate(sorted_transfer, x=sorted_taus) / L_sun) + " L_???")


