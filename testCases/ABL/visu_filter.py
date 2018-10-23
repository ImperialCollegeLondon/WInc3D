import numpy as np
import matplotlib.pyplot as plt

A=np.genfromtxt('filter.dat')
plt.plot(A[:,0],label='filtered')
plt.plot(A[:,1],label='unfiltered')
plt.legend()
plt.show()

