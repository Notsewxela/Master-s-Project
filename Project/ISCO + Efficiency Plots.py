import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
                               
def r_isco(a):
    '''Given the dimensionless spin parameter a from [-1, 1] calculates the
    innermost stable circular orbit in gravitational radii. This function cleverly
    uses special units to make equations simpler.'''
    #if abs(a) > 1:
    #    raise Exception("a must lie on the closed interval [-1, 1]")
    Z1 = 1 + (1-a**2)**(1/3)*((1+a)**(1/3)+(1-a)**(1/3))
    Z2 = (3*a**2 + Z1**2)**(1/2)
    Z3 = np.sign(a)*((3-Z1)*(3+Z1+2*Z2))**(1/2) 
    return 3 + Z2 - Z3

def eddington_efficiency(a):
    '''Calculates the eddington efficiency dependent on the inner accretion disk
    radius which is assumed to be at the ISCO. The ISCO depends on the dimensionless
    spin parameter a.'''
    return 1-(1-2/(3*r_isco(a)))**(1/2)

a = np.linspace(-1, 1, 10000)

fig, ax = plt.subplots(2, sharex='col', sharey='row', figsize = (10, 8))

ax[0].plot(a, r_isco(a), ls = "--", color = "red")
ax[0].set_ylabel('$r_{in}$', fontsize = 16)
ax[0].set_ylim([0, 10])
ax[1].set_ylim([0, 0.5])
ax[0].yaxis.set_minor_locator(AutoMinorLocator())
ax[1].yaxis.set_minor_locator(AutoMinorLocator())
ax[1].xaxis.set_minor_locator(AutoMinorLocator())



ax[1].plot(a, eddington_efficiency(a))
ax[1].set_xlabel('$a$', fontsize = 14)
ax[1].set_ylabel('$\\eta$', fontsize = 14)

ax[0].tick_params(axis='both', which='major', labelsize=12)
ax[1].tick_params(axis='both', which='major', labelsize=12)
plt.savefig("ISCO+Efficiency.pdf", bbox_inches = 'tight')
plt.show()
