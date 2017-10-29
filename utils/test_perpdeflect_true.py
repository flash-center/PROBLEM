import prad.deflect
import numpy as np
import matplotlib.pyplot as plt

wBx = np.loadtxt('wBx.dat', delimiter=',')
wBy = np.loadtxt('wBy.dat', delimiter=',')
x = np.loadtxt('x.dat', delimiter=',')
y = np.loadtxt('y.dat', delimiter=',')
ri = 1.3
li = .2
rs = 20
v = 5.24e9

magnify = (rs + ri + li)/(ri)
phix = (rs/(magnify*v))*wBx + x
phiy = (rs/(magnify*v))*wBy + y

true_flux = np.loadtxt('flux.dat', delimiter=',')
flux0 = np.mean(true_flux) * np.ones(true_flux.shape) # A reference, background 'no-deflections' flux image based on true_flux
flux_total = np.sum(true_flux) # Total number of protons in the flux image

#flux1 = prad.deflect.fluximage(ri, li, rs, v, x, y, flux_total, wBx, wBy)
#print(np.max(flux1))
flux2 = prad.deflect.fluximage2(x, y, phix, phiy, flux0, scale_fact=5, scale_order=3)
fig, ax = plt.subplots()
cs = ax.pcolorfast(x[:,0], y[0,:], flux2)
cb = plt.colorbar(cs, ax=ax)
#bound=9000
#cs = ax.pcolorfast(x[:,0], y[0,:], flux2, vmin=0, vmax=bound)
#v = np.linspace(0, bound, 10)
#cb = plt.colorbar(cs, ax=ax, ticks=v)
ax.set_title('Calculated Flux (count/bin)')
ax.set_xlabel('x (cm)')
ax.set_ylabel('y (cm)')
fig.savefig('test_perpdeflect_true2.png')

