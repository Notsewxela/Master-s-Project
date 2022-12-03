import numpy as np
from scipy.integrate import simps as integrate
import pylab
from math import pi
import scipy
from matplotlib import ticker, cm, colors

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


inclination = 30 #0 <= degrees < 90
mass_BH = 5e7 # mass in solar unit
mdot_eddington = 0.05 #accrestion rate in eddington units
spin = 1 #dimensionless spin parameter
r_in = r_isco(spin)+1e-15 #gravitational units
r_out = 1e3 #outer edge of disk
'''We have to add a small number to the innermost radius as other wise our T(r=r_isco) generates a 0 value
which causes our _inside_of_integral function to throw an error due to dividing by zero in the
exponential function. This change of r_min has no noticable effect on results.'''


radii_data_points = 1000
radii = np.geomspace(r_in+1e-15, r_out, num = radii_data_points) #in gravitational units

phi_data_points = 100 # needs to be >= to 1 for the code to work
angles = np.linspace(0, 2 * pi, phi_data_points)
#as the disk is symmetric we only need to do 0->pi and double the result
#angles are stored to be in the middle of segments

l_X = 0.08 #x ray luminosity as fraction of L_Bol
L_X = l_X*1e12*L_sun #here we take the total xray emission as a fraction of L_Bol taken to be ~1e12 suns
h_X = 3 * radii[0] #x ray source height above disk in gravitational units
albedo = 0.3 #albedo

freq_data_points = 1000
freq_low = 10**(5) #5
freq_high = 10**(17.5) #17.5
frequencies = np.geomspace(freq_low, freq_high, num = freq_data_points)

xs = np.zeros((len(radii), len(angles)))
ys = np.zeros((len(radii), len(angles)))
taus = np.zeros((len(radii), len(angles)))
surface_brightness = np.zeros((len(radii), len(angles)))
response = np.zeros((len(radii), len(angles)))

T1 = T(radii, mass_BH, mdot_eddington, spin, h_X, L_X, albedo, 0.0)
T2 = T(radii, mass_BH, mdot_eddington, spin, h_X, L_X, albedo, 1.0) #temperature is function of radius only

dr = np.append(np.delete(radii,0), radii[-1]*radii[-1]/radii[-2])-radii
dphi = 2*pi/len(angles)

for i in range(len(xs)):
    xs[i][0:phi_data_points] = radii[i]*np.sin(angles)
    ys[i][0:phi_data_points] = -radii[i]*np.cos(angles)*np.cos(inclination*pi/180)
    taus[i][0:phi_data_points] = D(radii[i], mass_BH, h_X, angles, inclination)/c/1000
    surface_brightness[i][0:phi_data_points] = 4*pi*integrate(B_nu(T2[i], frequencies) - B_nu(T1[i], frequencies), x = frequencies)

    dP_r_phi = surface_brightness[i][0] * ((G*M_sun*mass_BH)**2/c**4) * radii[i] * dr[i] * dphi
        
    for i_phi in range(len(response[i])):
        tau_r_plus_dr_phi = D(radii[i] + dr[i], mass_BH, h_X, angles[i_phi], inclination) / c
        tau_r_phi = D(radii[i], mass_BH, h_X, angles[i_phi], inclination) / c
        #(1/c) * (G*M_sun*m/c**2) * (radii[i]/sqrt(radii[i]**2+h_X**2)) - np.cos(angles[i_phi])*np.sin(inclination))*dr[i]
        response[i][i_phi] = dP_r_phi / abs(tau_r_plus_dr_phi - tau_r_phi)
        #response[i][i_phi] = dP_r_phi / abs((1/c)*(G*M_sun*mass_BH/c**2) * ((radii[i]/np.sqrt(radii[i]**2+h_X**2)) - np.cos(angles[i_phi])*np.sin(inclination)*dr[i]))
                    
response *= np.cos(inclination * pi/180)


pylab.figure(figsize = (10, 8*np.cos(inclination*pi/180)))
h=pylab.contourf(xs, ys, taus, 40, cmap = pylab.cm.gnuplot_r)
pylab.xlabel("x (R$_{g}$)", fontsize = 14)
pylab.ylabel("Projectecd y (R$_{g}$)", fontsize = 14)
pylab.title("Lag Contour Plot", fontsize = 18)
pylab.arrow(0, 0, 0, -radii[-1]/7, head_width = radii[-1]/20*np.cos(inclination*pi/180), color="white")
pylab.text(radii[-1]*4/7, radii[-1]*np.cos(inclination*pi/180)*0.85, "Inclination: " + str(inclination) + "°", fontsize=12)
cb=pylab.colorbar(h)
cb.ax.tick_params(labelsize=14)
cb.set_label("Lag (ks)", fontsize=14)
#pylab.xscale("log")
#pylab.yscale("log")
pylab.xticks(fontsize=12)
pylab.yticks(fontsize=12)
#pylab.savefig("Response-Function.pdf")
pylab.show()

pylab.figure(figsize = (10, 8*np.cos(inclination*pi/180)))
i=pylab.contourf(xs, ys, np.log10(surface_brightness), 100, cmap = pylab.cm.gnuplot2)
pylab.xlabel("x (R$_{g}$)", fontsize = 14)
pylab.ylabel("Projectecd y (R$_{g}$)", fontsize = 14)
pylab.title("Surface Brightness Contour Plot", fontsize = 18)
pylab.arrow(0, 0, 0, -radii[-1]/7, head_width = radii[-1]/20*np.cos(inclination*pi/180), color="white")
pylab.text(radii[-1]*4/7, radii[-1]*np.cos(inclination*pi/180)*0.85, "Inclination: " + str(inclination) + "°", fontsize=12)
cb=pylab.colorbar(i)
cb.ax.tick_params(labelsize=14)
cb.set_label("log$_{10}$ Surface Brightness", fontsize=14)
#pylab.xscale("log")
#pylab.yscale("log")
pylab.xticks(fontsize=12)
pylab.yticks(fontsize=12)
#pylab.savefig("Response-Function.pdf")
pylab.show()

#print("Total observed extra luminosity from x-ray heating " + "{:e}".format(integrate(transfer_unnice, x=taus_unnice) / L_sun) + " L_☉")

pylab.figure(figsize = (10, 8*np.cos(inclination*pi/180)))
i=pylab.contourf(xs, ys, np.log10(response), 100, cmap = pylab.cm.gnuplot2)
pylab.xlabel("x (R$_{g}$)", fontsize = 14)
pylab.ylabel("Projectecd y (R$_{g}$)", fontsize = 14)
pylab.title("Response function contour plot", fontsize = 18)
pylab.arrow(0, 0, 0, -radii[-1]/7, head_width = radii[-1]/20*np.cos(inclination*pi/180), color="white")
pylab.text(radii[-1]*4/7, radii[-1]*np.cos(inclination*pi/180)*0.85, "Inclination: " + str(inclination) + "°", fontsize=12)
cb=pylab.colorbar(i)
cb.ax.tick_params(labelsize=14)
cb.set_label("log$_{10}$ Response (Ws-1)", fontsize=14)
#pylab.xscale("log")
#pylab.yscale("log")
pylab.xticks(fontsize=12)
pylab.yticks(fontsize=12)
#pylab.savefig("Response-Function.pdf")
pylab.show()

for i in range(30):
    pylab.plot(np.linspace(0, phi_data_points, phi_data_points), response[i], label = "Inclination: " + str(inclination) + "°")
pylab.xlabel("Angle as portions of pi", fontsize = 14)
pylab.ylim(0, 2e32)
pylab.show()

for i in range(30):
    pylab.plot(np.linspace(0, radii_data_points, radii_data_points), response.T[i], label = "Inclination: " + str(inclination) + "°")
pylab.xlabel("Radius portions", fontsize = 14)
pylab.ylim(0, 2e32)
pylab.show()


