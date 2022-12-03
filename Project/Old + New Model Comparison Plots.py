import pylab
import numpy as np

data1 = np.loadtxt("Data Files/response_spin_1_1.00E+05-3.16E+17.dat")
data2 = np.loadtxt("Data Files/response_spin_1_1.00E+05-3.16E+17_OLD_MODEL.dat")

taus1 = data1[0]
resp1 = data1[1]

taus2 = data2[0] #should equal taus1
resp2 = data2[1]

pylab.figure(figsize = (10, 8))
pylab.plot(taus1, resp1, label = "Optical Depth Considered")
pylab.plot(taus2, resp2, label = "Optical Depth Not Considered")
pylab.xscale("log")
pylab.title("Response Function Comparison", fontsize = 16)
pylab.xlabel("$\\tau$ (s)", fontsize = 14)
pylab.ylabel("Response (Ws$^{-1}$)", fontsize = 14)
pylab.xticks(fontsize = 12)
pylab.yticks(fontsize = 12)
pylab.legend(fontsize = 12)
pylab.savefig("New-Old_Model_Comparison.pdf")
pylab.show()