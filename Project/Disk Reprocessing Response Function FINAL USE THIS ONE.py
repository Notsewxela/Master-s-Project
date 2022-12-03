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

def B_nu(T, nu):
    '''Calculates the blackbody function WHz^-1m^-2Sr^-1'''
    if T.any() < 1e-10:
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

def response(r, phis, m, mdot_E, a, h_x, L_x, A, inc_deg):
    transfer = np.zeros((radii_data_points, phi_data_points)) #transfer function for each segment
    taus = np.zeros((radii_data_points, phi_data_points)) #delay at each point

    dphi = 2*pi/len(phis)
    #arrays of the temperatures as functions of radius
    T1 = T(r, m, mdot_eddington, a, h_x, L_x, A, 0.0)
    T2 = T(r, m, mdot_eddington, a, h_x, L_x, A, 1.0)
    
    dr = np.append(np.delete(radii,0), radii[-1]*radii[-1]/radii[-2])-radii #dr = r[i+1]-r[i] but smrat
    
    #radius and angle are index variables in these loops
    for radius in range(len(transfer)):
        Bs_tot = integrate(B_nu(T2[radius], frequencies) - B_nu(T1[radius], frequencies), x = frequencies)
        for angle in range(len(transfer[radius])):
            area = (G*M_sun*m)**2/c**4 * r[radius] * dr[radius] * dphi
            
            lag_spread = abs(D(r[radius+1], m, h_x, phis[angle], inclination)\
                   - D(r[radius], m, h_x, phis[angle], inclination)) / c
                             
            transfer[radius][angle] = 4*pi*Bs_tot * area / lag_spread
            
        taus[radius][0:len(phis)] = D(r[radius], m, h_x, phis, inclination) / c
    
    transfer *= np.cos(inclination * pi/180)
    
    flat_taus = taus.flatten()
    flat_transfer = transfer.flatten()

    tausinds = flat_taus.argsort()
    sorted_taus = flat_taus[tausinds[::1]]
    sorted_transfer = flat_transfer[tausinds[::1]]
    return sorted_transfer, sorted_taus

def binify(r, phis, m, h_x, inc_deg, transfer, taus):
    #takes out response function for each segment and bins them to make plots 10/10
    binsplus1 = np.ceil(len(r))
    taumin = D(r[0], m, h_x, 0, inc_deg)/c
    taumax = D(r[-1], m, h_x, pi, inc_deg)/c
    nice_taus0 = np.geomspace(taumin, taumax*1.00001, binsplus1) #the edges of our bins
    nice_transfer0 = np.zeros(len(nice_taus0))
    indices = np.digitize(taus, nice_taus0, right=False) #bin[i-1] <= x < bin[i]
    for i in range(len(indices)):
        nice_transfer0[indices[i]] += transfer[i]/len(np.where(indices == indices[i])[0]) #[0] is because np.where is odd
    
    nice_transfer0 *= len(phis) #somehow correctly? normalises it

    nice_transfer = np.insert(np.delete(nice_transfer0, np.where(nice_transfer0 == 0)), 0, 0) #remove 0s
    nice_taus = np.insert(np.delete(nice_taus0, np.where(nice_transfer0 == 0)), 0, taumin)
    
    #N = 20
    #return np.convolve(nice_transfer, np.ones(N)/N, mode='valid'), np.convolve(nice_taus, np.ones(N)/N, mode='valid')
    return nice_transfer, nice_taus

freq_data_points = 1000 #increase this to get the sum of the transfer function data points be
                            #closer to the area under the disk spectrum

inclination = 0 #0 <= degrees < 90
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

radii_data_points = 1000
radii = np.geomspace(r_in+1e-15, r_out, num = radii_data_points + 1) #in gravitational units
#extra data point is so we can do r[i+1]-r[i] for all the data point we want to test. it is only for this use!

phi_data_points = 100 #needs to be >= to 1 for the code to work
angles = np.linspace(pi/(2*phi_data_points), pi*(2*phi_data_points-1)/(2*phi_data_points), phi_data_points)
#as the disk is symmetric we only need to do 0->pi and double the result
#angles are stored to be in the middle of segments

l_X = 0.08 #x ray luminosity as fraction of L_Bol
L_X = l_X*1e12*L_sun #here we take the total xray emission as a fraction of L_Bol taken to be ~1e12 suns
h_X = 3 * radii[0] #x ray source height above disk in gravitational units
albedo = 0.3 #albedo

transfer_unbinned, taus_unbinned = response(radii, angles, mass_BH, mdot_eddington, spin, h_X, L_X, albedo, inclination)

transfer, taus = binify(radii, angles, mass_BH, h_X, inclination, transfer_unbinned, taus_unbinned)

pylab.figure(figsize = (10, 8))
pylab.title("Response Functions", fontsize = 16)
pylab.xlabel("Lag (s)", fontsize = 14)
pylab.ylabel("Response Function (Ws$^{-1}$)", fontsize = 14)
pylab.plot(taus[2:], transfer[2:], label = "Inclination: " + str(inclination) + "°")
pylab.xscale("log")
#pylab.yscale("log")
pylab.xticks(fontsize=12)
pylab.yticks(fontsize=12)
pylab.legend(fontsize = 14)
#pylab.savefig("Response-Function.pdf")
pylab.show()

#print("Total observed extra luminosity from x-ray heating " + "{:e}".format(integrate(transfer_unnice, x=taus_unnice) / L_sun) + " L_☉")











