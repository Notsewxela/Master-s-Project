import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator, FormatStrFormatter



fig, axs = plt.subplots(7, sharex=True, sharey=False, gridspec_kw={'hspace': 0}, figsize = (10, 8))



xfile = "Data Files/xrays_0.5-10keV.dat"

data = np.loadtxt(xfile)
x = data.T[0]
y = data.T[1]
err = data.T[2]
axs[0].errorbar(x, y, yerr=err, ls = "", color = "black")
axs[0].set_ylabel("X-Rays", fontsize = 14)
axs[0].yaxis.set_minor_locator(AutoMinorLocator())
axs[0].yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
axs[-1].xaxis.set_minor_locator(AutoMinorLocator())


filenames = ["w2", "m2", "w1", "u", "b", "v"]
colors = ["midnightblue", "navy", "darkblue", "mediumblue", "blue", "green"]


for i, f in enumerate(filenames):
    data = np.loadtxt("Data Files/" + f + "_mjy_final.dat")
    x = data.T[0]
    y = data.T[1]
    err = data.T[2]
    axs[i+1].errorbar(x, y, yerr=err, ls="", color = colors[i])
    axs[i+1].set_ylabel(f.upper(), fontsize = 14)
    axs[i+1].yaxis.set_minor_locator(AutoMinorLocator())
    axs[i+1].yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
axs[-1].yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
# Hide x labels and tick labels for all but bottom plot.
for ax in axs:
    ax.label_outer()
axs[-1].set_xlabel("Time (MJD - 50000)", fontsize = 14)

plt.savefig("lightcurves.pdf", bbox_inches = 'tight')