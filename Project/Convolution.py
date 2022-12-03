import numpy as np
import pylab
from scipy.integrate import simps as integrate
from math import pi

import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator, FormatStrFormatter

def timesortingout(t_x, x, t_r, r, N_r = 200):
    '''Assumes sorted lists in order of increasing time with same units for time! N_r is the number of equally
    spaced points (in time) we would like our response function to have. x stands for x-ray!'''
    t_max_x = t_x[-1]
    t_min_x = t_x[0]
    
    t_max_r = t_r[-1]
    t_min_r = t_r[0]
    
    length_ratio = (t_max_x-t_min_x)/(t_max_r-t_min_r)
    N_x = int(np.floor(length_ratio * N_r)) #have to floor for it to work, as long as N_r is big should't affect results noticably
    
    
    new_t_x = np.linspace(t_min_x, t_max_x, N_x)
    new_t_r = np.linspace(t_min_r, t_max_r, N_r)
    
    new_x = np.interp(new_t_x, t_x, x)
    new_r = np.interp(new_t_r, t_r, r)
    
    
    return new_t_x, new_x, new_t_r, new_r

'''
data = np.loadtxt("Data Files/xrays_0.5-10keV.dat")
xray_times = data.T[0]
xray_flux = data.T[1]
xray_flux_err = data.T[2]

pylab.figure(figsize = (10, 8))
pylab.errorbar(xray_times, xray_flux, xray_flux_err)
pylab.title("X-ray 0.5-10 keV Light Curve", fontsize = 16)
pylab.xlabel("Time (MJD - 50000, Days)", fontsize = 14)
pylab.ylabel("Flux (mJy)", fontsize = 14)
pylab.xticks(fontsize = 12)
pylab.yticks(fontsize = 12)
pylab.show()

resp = np.loadtxt("Data Files/response_spin_1.dat")
resp_times = resp[0] / (60 * 60 * 24) #times converted to days
resp_flux = resp[1] / np.amax(resp[1]) #responses normalised


pylab.figure(figsize = (10, 8))
pylab.plot(resp_times, resp_flux)
pylab.title("Response Function", fontsize = 16)
pylab.xlabel("Time (Days)", fontsize = 14)
pylab.ylabel("Response (Normalised By Area)", fontsize = 14)
pylab.xticks(fontsize = 12)
pylab.yticks(fontsize = 12)
pylab.show()

xray_times_interp, xray_flux_interp, resp_times_interp, resp_flux_interp = timesortingout(xray_times, xray_flux, resp_times, resp_flux)


pylab.figure(figsize = (10, 8))
pylab.plot(xray_times_interp, xray_flux_interp)
pylab.title("Interpolated X-ray 0.5-10 keV Light Curve", fontsize = 16)
pylab.xlabel("Time (MJD - 50000, Days)", fontsize = 14)
pylab.ylabel("Flux (mJy)", fontsize = 14)
pylab.xticks(fontsize = 12)
pylab.yticks(fontsize = 12)
pylab.show()



pylab.figure(figsize = (10, 8))
pylab.plot(resp_times_interp, resp_flux_interp)
pylab.title("Interpolated Response Function", fontsize = 16)
pylab.xlabel("Time (Days)", fontsize = 14)
pylab.ylabel("Response (Normalised to Peak)", fontsize = 14)
pylab.xticks(fontsize = 12)
pylab.yticks(fontsize = 12)
pylab.show()

#print(resp_times_interp[1] - resp_times_interp[0])
#print(xray_times_interp[1] - xray_times_interp[0])

model_UVW2 = np.convolve(resp_flux_interp, xray_flux_interp)
model_UVW2_times = np.linspace(xray_times[0], xray_times[-1], len(model_UVW2))


pylab.figure(figsize = (10, 8))
pylab.plot(model_UVW2_times, model_UVW2)
pylab.title("Predicted UVW2 Light Curve", fontsize = 16)
pylab.xlabel("Time (MJD - 50000, Days)", fontsize = 14)
pylab.ylabel("Flux (mJy)", fontsize = 14)
pylab.xticks(fontsize = 12)
pylab.yticks(fontsize = 12)
pylab.show()

np.savetxt("Data Files/Predicted UVW2 Light Curve.dat", (model_UVW2_times, model_UVW2))
'''




fig, axs = plt.subplots(6, sharex=True, sharey=False, gridspec_kw={'hspace': 0}, figsize = (10, 8))

xfile = "Data Files/xrays_0.5-10keV.dat"

data = np.loadtxt(xfile)
xray_times = data.T[0]
xray_flux = data.T[1]

axs[0].set_title("Predicted Light Curves", fontsize = 16)
axs[0].yaxis.set_minor_locator(AutoMinorLocator())
axs[0].yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
axs[-1].xaxis.set_minor_locator(AutoMinorLocator())

filenames = ["w2", "m2", "w1", "u", "b", "v"]
colors = ["midnightblue", "navy", "darkblue", "mediumblue", "blue", "green"]


for i, f in enumerate(filenames):
   
    resp = np.loadtxt("Data Files/response_" + f.upper() + ".dat")
    resp_times = resp[0] / (60 * 60 * 24) #times converted to days
    resp_flux = resp[1] / np.amax(resp[1]) #responses normalised

    
    xray_times_interp, xray_flux_interp, resp_times_interp, resp_flux_interp = timesortingout(xray_times, xray_flux, resp_times, resp_flux)
    model_flux = np.convolve(resp_flux_interp, xray_flux_interp)
    model_flux /= np.amax(model_flux)
    model_times = np.linspace(xray_times[0], xray_times[-1], len(model_flux))
    
    axs[i].plot(model_times, model_flux, ls = "", marker = "x", markersize = 0.2, color = colors[i])
    axs[i].set_ylabel(f.upper(), fontsize = 14)
    axs[i].yaxis.set_minor_locator(AutoMinorLocator())
    axs[i].yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    axs[i].set_ylim(-0.1, 1.1)
    np.savetxt("Data Files/predicted_lightcurve_" + f.upper() + ".dat", (model_times, model_flux))
# Hide x labels and tick labels for all but bottom plot.
for ax in axs:
    ax.label_outer()
axs[-1].set_xlabel("Time (MJD - 50000)", fontsize = 14)

#plt.savefig("predited_lightcurves.pdf", bbox_inches = 'tight')


