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

def T(r, m, mdot_E, a):
    '''Calculates the temperature of the accretion disk at a radius r (in
    gravitational units) given the black hole mass (in solar masses) and the
    accretion rate (in eddington units).'''
    #r is radius in gravitational units
    #m is the BH mass in solar masses
    #mdot_E is the accretion rate in eddington units
    eta = eddington_efficiency(a) #efficiency
    const = (3*m_p*c**5/(2*eta*sigma_T*sigma_SB*G*M_sun))**0.25
    extra_r_factor = 1-(r_isco(a)/r)**(1/2)
    return const * (mdot_E * extra_r_factor / (r**3*m))**0.25

def L_nu(nu, radii, inc_deg, m, mdot_E, a):
    '''Calculates the spectral flux at a frequency nu for our accretion disk
    contained within radii's min and max values given the accretion disk
    inclination, mass of BH in solar units and eddington accretion rate as well
    as the dimensionless spin parameter a.'''
    #nu is the frequency
    #inc_deg is the inclination in degrees, 0 is face on disk
    #radii are the array of radii we are measuring T at
    def _inside_of_integral(nu, r, m, mdot_E):
        '''Calculates the inside the spectral flux integral.'''
        return 2*pi*r/(np.exp(h*nu/(k_B*T(r, m, mdot_E, a)))-1)
    
    inc = inc_deg * pi/180
    non_integral_bit = 4*pi*2*h*nu**3*np.cos(inc)*G**2*M_sun**2*m**2 / c**6
    #https://personal.sron.nl/~kaastra/leiden2017/lnotes_part7.pdf for the extra 2pi on the end and explanation
    #2 from planck, 2pi from disk area, 2 from 2 sides, 2 from integration over half of physical space
    bit_to_integrate = _inside_of_integral(nu, radii, m, mdot_E)
    integral_bit = integrate(bit_to_integrate, x = radii)
    return non_integral_bit*integral_bit

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

inclination = 0 #degrees
mass_BH = 5e7 # mass in solar unit
mdot_eddington = 0.5 #accretion rate in eddington units
a = 1 #dimensionless spin parameter
r_in = r_isco(a)+1e-15
r_out = 1e5 #outer edge of disk

freq_low = 10**(5)
freq_high = 10**(17.5)
freq_data_points = 1000
frequencies = np.geomspace(freq_low, freq_high, num = freq_data_points)
specific_intensities = np.zeros(freq_data_points)

radii_data_points = 1000
radii = np.geomspace(r_in, r_out, num = radii_data_points) #in gravitational units
'''We have to add a small number to the innermost radius as other wise our T(r=r_isco) generates a 0 value
which causes our _inside_of_integral function to throw an error due to dividing by zero in the
exponential function. This change of r_min has no noticable effect on results.'''

#here we loop through all our frequencies we would like to find the strength of by at each one integrating over all radii
for i in range(len(specific_intensities)):
    specific_intensities[i] = L_nu(frequencies[i], radii, inclination, mass_BH, mdot_eddington, a)
    
pylab.figure(figsize = (10, 8))
pylab.title("Disk Spectrum with No X-Ray Reprocessing", fontsize = 16)
pylab.xlabel("$\\nu$ (Hz)", fontsize = 14)
pylab.ylabel("Spectral Flux (WHz$^{-1}$)", fontsize = 14)
pylab.loglog(frequencies, specific_intensities, label = "$L_\\nu$")
pylab.legend(fontsize = 14)

print("{:e}".format(integrate(specific_intensities, x = frequencies) / L_sun) + " L_☉")