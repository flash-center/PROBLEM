import numpy as np
import matplotlib.pyplot as plt

flux = np.loadtxt('flux.dat', delimiter=',')
x = np.loadtxt('x.dat', delimiter=',')
y = np.loadtxt('y.dat', delimiter=',')
fig, ax = plt.subplots()

cs = ax.pcolorfast(x[:,0], y[0,:], flux, vmin=0, vmax=5500)
cb = plt.colorbar(cs, ax=ax)
ax.set_title('True Flux (count/bin)')
ax.set_xlabel('x (cm)')
ax.set_ylabel('y (cm)')
fig.savefig('true_flux.png')
