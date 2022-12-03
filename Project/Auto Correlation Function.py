import numpy as np
import pylab
from math import pi
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator, FormatStrFormatter
from scipy.signal import chirp, find_peaks, peak_widths

def UDCCF_calc(x_times, x_data, y_times, y_data):
    '''Evaluate our constants'''
    x_mean = np.mean(x_data)
    y_mean = np.mean(y_data)
    x_std = np.std(x_data)
    y_std = np.std(y_data)
    
    '''Creating empty arrays to store our UDCCF and tau data'''
    UDCCF_non_normalised = np.zeros((len(x_data), len(y_data)))
    taus = np.zeros((len(x_data), len(y_data)))
    
    '''Now i've done in smratly so doen't need any loops and it works using clever numpy magic! Much faster for big arrays'''
    y_times =  np.array([y_times]).T
    
    y_diffs = np.array([y_data - y_mean]).T
    x_diffs = x_data - x_mean
    
    UDCCF_non_normalised = x_diffs*y_diffs
    taus = x_times-y_times
    
    '''Divide by the standard deviations of the data sets'''
    UDCCF = UDCCF_non_normalised / (x_std * y_std)
    
    '''Flattens our 2d arrays of lag differences and the UDCCFs into a long 1D array'''
    flat_taus = taus.flatten() 
    flat_UDCCF = UDCCF.flatten()
    
    '''Sorts our arrays based off the corresponding tau value in increasing tau'''
    taus_indices = flat_taus.argsort()
    sorted_taus = flat_taus[taus_indices[::1]]
    sorted_UDCCF = flat_UDCCF[taus_indices[::1]]
    
    return sorted_taus, sorted_UDCCF

def DCCF_calc(sorted_taus, sorted_UDCCF, num_bins):
    '''Sorts out our binning parameters'''
    tau_min = sorted_taus[0]
    tau_max = sorted_taus[-1]
    bin_width = (tau_max - tau_min) / num_bins
    
    DCCF = np.zeros(num_bins)
    DCCF_errs = np.zeros(num_bins)
    bin_lims = np.linspace(tau_min, tau_max, num_bins + 1)
    bin_lims[0] -= abs(tau_min) * 1e-7
    bin_lims[-1] += abs(tau_max) * 1e-7
    
    for i in range(num_bins):
        '''The first line gets all the values of UDCCF which have a tau in the correct range'''
        temp = sorted_UDCCF[(bin_lims[i] < sorted_taus) & (sorted_taus < bin_lims[i+1])]
        
        DCCF[i] = np.mean(temp)
        DCCF_errs[i] = np.std(temp)*np.sqrt(2)/np.sqrt(len(temp))
    
    xpositions = bin_lims[:-1] + bin_width/2
    return xpositions, DCCF, DCCF_errs

def FWHM(old_x, old_y, bins, trim, height = 0.5):
    #Main peak is obviously going to be at time=0 or the ACF has gone wrong
    '''This function "fixes" the data so that scipy.peak_widths works properly *eyeroll emoji*'''
    
    x = np.copy(old_x)#if you don't copy to new arrays it all breaks when trying FWHM multiple times with same data set
    y = np.copy(old_y)
    peak_index = np.array([(bins-1-trim-trim)//2])
    left_most = np.amin(x)
    x += abs(left_most)
    dx = x[1]-x[0]
    x /= dx

    #print(peak_index, xpositions[peak_index], DCCF[peak_index], dx)
    widths, width_heights, left_ips, right_ips = peak_widths(y, peak_index, rel_height=height)
    widths *= dx
    left_ips = left_ips * dx - abs(left_most) #"unfix" them
    right_ips = right_ips * dx - abs(left_most)    
    
    return widths, width_heights, left_ips, right_ips
    
    

filenames = ["w2", "m2", "w1", "u", "b", "v"]
colors = ["midnightblue", "navy", "darkblue", "mediumblue", "blue", "green"]

fig, axs = plt.subplots(6, sharex=True, sharey=False, gridspec_kw={'hspace': 0}, figsize = (10, 8))

axs[0].set_title("Auto Correlation Functions", fontsize = 16)
axs[0].yaxis.set_minor_locator(AutoMinorLocator())
axs[0].yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
axs[-1].xaxis.set_minor_locator(AutoMinorLocator())

bin_num1 = 101 #makes these odd number so there is bin at 0
bin_num2 = 101
trim = 10 #>=1
for i, f in enumerate(filenames):
    data = np.loadtxt("Data Files/" + f + "_mjy_final.dat")
    xdataunnice = data.T[0]
    ydataunnice = data.T[1]
    del data
    model = np.loadtxt("Data Files/predicted_lightcurve_" + f.upper() + ".dat")
    xmodelunnice = model[0]
    ymodelunnice = model[1]
    del model
    
    xdata = np.linspace(xdataunnice[0], xdataunnice[-1], int(np.floor(len(xdataunnice))))
    ydata = np.interp(xdata, xdataunnice, ydataunnice)
    
    xmodel = np.linspace(xmodelunnice[0], xmodelunnice[-1], len(xmodelunnice))
    ymodel = np.interp(xmodel, xmodelunnice, ymodelunnice)

    sorted_taus, sorted_UDCCF = UDCCF_calc(xdata, ydata, xdata, ydata)
    xpositions, DCCF, DCCF_errs = DCCF_calc(sorted_taus, sorted_UDCCF, bin_num1)
    xpositions = xpositions[trim:-trim]
    DCCF = DCCF[trim:-trim]
    DCCF_errs = DCCF_errs[trim:-trim]
    axs[i].errorbar(xpositions, DCCF, yerr=DCCF_errs, ls = "", marker = ".", markersize = 2, color = colors[i], label = "Data")
    
    heights = [0.25, 0.5, 0.75]
    labels = ["fwqm", "fwhm", "fw3qm"]
    for j, h in enumerate(heights):
        widths, width_heights, left_ips, right_ips = FWHM(xpositions, DCCF, bin_num1, trim, height = h)
        print(f + " data " + labels[j] + " = " + "{:.2f}".format(widths[0]) + " $\pm$ " +\
              "{:.2f}".format(FWHM(xpositions, DCCF, bin_num1, trim, height = h+0.05)[0][0] - FWHM(xpositions, DCCF, bin_num1, trim, height = h-0.05)[0][0]))
        axs[i].hlines(y = width_heights, xmin = left_ips, xmax = right_ips, color = colors[i])
    #widths, width_heights, left_ips, right_ips = FWHM(xpositions, DCCF, bin_num1, trim, height = 1)
    #print(f + " data FWM = " + "{:.2f}".format(widths[0]) + " $\pm$ " +\
    #          "{:.2f}".format(2*(FWHM(xpositions, DCCF, bin_num1, trim, height = h*1)[0][0] - FWHM(xpositions, DCCF, bin_num1, trim, height = h*0.90)[0][0])))
    axs[i].hlines(y = width_heights, xmin = left_ips, xmax = right_ips, color = colors[i])
    
    
    sorted_taus, sorted_UDCCF = UDCCF_calc(xmodel, ymodel, xmodel, ymodel)
    xpositions, DCCF, DCCF_errs = DCCF_calc(sorted_taus, sorted_UDCCF, bin_num2)
    xpositions = xpositions[trim:-trim]
    DCCF = DCCF[trim:-trim]
    DCCF_errs = DCCF_errs[trim:-trim]
    axs[i].errorbar(xpositions, DCCF, yerr=DCCF_errs, ls = "", marker = "x", markersize = 2, color = "red", label = "Model")
    
    for j, h in enumerate(heights):
        widths, width_heights, left_ips, right_ips = FWHM(xpositions, DCCF, bin_num2, trim, height = h)
        print(f + " model " + labels[j] + " = " + "{:.2f}".format(widths[0]) + " $\pm$ " +\
              "{:.2f}".format(FWHM(xpositions, DCCF, bin_num2, trim, height = h+0.05)[0][0] - FWHM(xpositions, DCCF, bin_num2, trim, height = h-0.05)[0][0]))
        axs[i].hlines(y = width_heights, xmin = left_ips, xmax = right_ips, color = "red")
    #widths, width_heights, left_ips, right_ips = FWHM(xpositions, DCCF, bin_num2, trim, height = 1)
    #print(f + " model FWM = " + "{:.2f}".format(widths[0]) + " $\pm$ " +\
    #          "{:.2f}".format(2*(FWHM(xpositions, DCCF, bin_num2, trim, height = h*1)[0][0] - FWHM(xpositions, DCCF, bin_num2, trim, height = h*0.90)[0][0])))
    axs[i].hlines(y = width_heights, xmin = left_ips, xmax = right_ips, color = colors[i])
    
    print("\n")
    
    
    axs[i].set_ylabel(f.upper(), fontsize = 14)
    axs[i].yaxis.set_minor_locator(AutoMinorLocator())
    axs[i].yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    axs[i].legend()
    
# Hide x labels and tick labels for all but bottom plot.
for ax in axs:
    ax.label_outer()
axs[-1].set_xlabel("Time (Days)", fontsize = 14)

#plt.savefig("autocorrelation.pdf", bbox_inches = 'tight')
