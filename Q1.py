# Q1 is for plotting results from the critical_mass.f90 simulation

# Import packages
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from matplotlib.ticker import FixedFormatter

# Read the file and assign the data to numpy arrays
f = open("problem1flat19.dat", 'r')
generationf, neutronsf = np.loadtxt(f, usecols = (0, 1), unpack = True)
f.close()
f = open("problem1sub16.2.dat", 'r')
generationSub, neutronsSub = np.loadtxt(f, usecols = (0, 1), unpack = True)
f.close()
f = open("problem1super23.dat", 'r')
generationSup, neutronsSup = np.loadtxt(f, usecols = (0, 1), unpack = True)
f.close()

N0 = neutronsf[0]

# Plot
ax = plt.subplot(1, 1, 1)
plt.xlabel(r'neutron cycle')
plt.ylabel(r'neutrons / N(0)')
ax.set_yscale('log')
ax.set_ylim(0.001, 1000.0)
ax.yaxis.set_major_formatter(FixedFormatter(['0.001', '0.010', '0.100', '1.0', '10.0', '100.0', '1000.0']))
ax.set_yticks([0.001, 0.01, 0.1, 1.0, 10.0, 100.0, 1000.0])
ax.xaxis.set_minor_locator(MultipleLocator(10))
plt.plot(generationf, neutronsf/N0, ':', color='black')
plt.plot(generationSub, neutronsSub/N0, '--', color='black')
plt.plot(generationSup, neutronsSup/N0, '-', color='black')
plt.title(r'Critical Mass of Uranium Sphere')
plt.show()
