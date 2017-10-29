"""
Throw-away test of various flux generation methods. Created by Scott 2017-10-17
"""

import sys, os

try:
    import prad.deflect
except:
    sys.path.append(os.path.realpath(".."))
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

#magnify = (rs + ri + li)/(ri)
magnify = (rs+li+ri)/(ri+.5*li)


phix = (rs/(magnify*v))*wBx + x
phiy = (rs/(magnify*v))*wBy + y

flux_true = np.loadtxt('flux.dat', delimiter=',')
flux0 = np.mean(flux_true) * np.ones(flux_true.shape) # A reference, background 'no-deflections' flux image based on true_flux
nprotons = np.sum(flux_true) # Total number of protons in the flux image


## PERFORM THE THREE METHODS
#print("METHOD 1, GO!")
#flux1 = prad.deflect.fluximage(ri, li, rs, v, x, y, nprotons, wBx, wBy)

print("METHOD 2, GO!")
flux2 = prad.deflect.fluximage2(x, y, phix, phiy, flux0, scale_fact=20, scale_order=3)

print("METHOD 3, GO!")
flux3, flux3_0 = prad.deflect.fluximage3(ri, li, rs, v, x, y, nprotons, wBx, wBy, Ntest=5000000)
#

### PLOTS
vmin = 0
vmax = 5500
cmap = 'viridis'

fig = plt.figure(1)
fig.clf()
ax = fig.add_subplot(111)
cs = ax.pcolorfast(x[:,0], y[0,:], flux_true, vmin=vmin, vmax=vmax, cmap=cmap)
cb = plt.colorbar(cs, ax=ax)
ax.set_title('"True flux" (count/bin)')
ax.set_xlabel('x (cm)')
ax.set_ylabel('y (cm)')
ax.set_aspect('equal')
fig.savefig("FluxMethods_True.png")

fig = plt.figure(2)
fig.clf()
ax = fig.add_subplot(111)
cs = ax.pcolorfast(x[:,0], y[0,:], flux1, vmin=vmin, vmax=vmax, cmap=cmap)
cb = plt.colorbar(cs, ax=ax)
ax.set_title('Method 1 (count/bin)')
ax.set_xlabel('x (cm)')
ax.set_ylabel('y (cm)')
ax.set_aspect('equal')
fig.savefig("FluxMethods_1.png")


fig = plt.figure(3)
fig.clf()
ax = fig.add_subplot(111)
cs = ax.pcolorfast(x[:,0], y[0,:], flux2, vmin=vmin, vmax=vmax, cmap=cmap)
cb = plt.colorbar(cs, ax=ax)
ax.set_title('Method 2 (count/bin)')
ax.set_xlabel('x (cm)')
ax.set_ylabel('y (cm)')
ax.set_aspect('equal')
fig.savefig("FluxMethods_2.png")

fig = plt.figure(4)
fig.clf()
ax = fig.add_subplot(111)
cs = ax.pcolorfast(x[:,0], y[0,:], flux3, vmin=vmin, vmax=vmax, cmap=cmap)
cb = plt.colorbar(cs, ax=ax)
ax.set_title('Method 3 (count/bin)')
ax.set_xlabel('x (cm)')
ax.set_ylabel('y (cm)')
ax.set_aspect('equal')
fig.savefig("FluxMethods_3.png")


fig = plt.figure(5)
fig.clf()
ax = fig.add_subplot(111)
cs = ax.pcolorfast(x[:,0], y[0,:], flux1 / flux_true, vmin=0.0, vmax=2.0, cmap='RdBu')
cb = plt.colorbar(cs, ax=ax)
ax.set_title('Flux1 / True')
ax.set_xlabel('x (cm)')
ax.set_ylabel('y (cm)')
ax.set_aspect('equal')
fig.savefig("FluxMethods_1r.png")

fig = plt.figure(6)
fig.clf()
ax = fig.add_subplot(111)
cs = ax.pcolorfast(x[:,0], y[0,:], flux2 / flux_true, vmin=0.0, vmax=2.0, cmap='RdBu')
cb = plt.colorbar(cs, ax=ax)
ax.set_title('Flux2 / True')
ax.set_xlabel('x (cm)')
ax.set_ylabel('y (cm)')
ax.set_aspect('equal')
fig.savefig("FluxMethods_2r.png")

fig = plt.figure(7)
fig.clf()
ax = fig.add_subplot(111)
cs = ax.pcolorfast(x[:,0], y[0,:], flux3 / flux_true, vmin=0.0, vmax=2.0, cmap='RdBu')
cb = plt.colorbar(cs, ax=ax)
ax.set_title('Flux3 / True')
ax.set_xlabel('x (cm)')
ax.set_ylabel('y (cm)')
ax.set_aspect('equal')
fig.savefig("FluxMethods_3r.png")


vmin = -1.25
vmax = 1.25
fig = plt.figure(8)
fig.clf()
ax = fig.add_subplot(111)
cs = ax.pcolorfast(x[:,0], y[0,:], np.log10(flux_true/flux0), vmin=vmin, vmax=vmax, cmap='RdBu')
cb = plt.colorbar(cs, ax=ax, label="Log10(Flux/Reference)")
ax.set_title('Contrast: True')
ax.set_xlabel('x (cm)')
ax.set_ylabel('y (cm)')
ax.set_aspect('equal')
fig.savefig("FluxMethods_LogTrue.png")

fig = plt.figure(9)
fig.clf()
ax = fig.add_subplot(111)
cs = ax.pcolorfast(x[:,0], y[0,:], np.log10(flux1/flux0), vmin=vmin, vmax=vmax, cmap='RdBu')
cb = plt.colorbar(cs, ax=ax, label="Log10(Flux/Reference)")
ax.set_title('Contrast: 1')
ax.set_xlabel('x (cm)')
ax.set_ylabel('y (cm)')
ax.set_aspect('equal')
fig.savefig("FluxMethods_Log1.png")

fig = plt.figure(10)
fig.clf()
ax = fig.add_subplot(111)
cs = ax.pcolorfast(x[:,0], y[0,:], np.log10(flux2/flux0), vmin=vmin, vmax=vmax, cmap='RdBu')
cb = plt.colorbar(cs, ax=ax, label="Log10(Flux/Reference)")
ax.set_title('Contrast: 2')
ax.set_xlabel('x (cm)')
ax.set_ylabel('y (cm)')
ax.set_aspect('equal')
fig.savefig("FluxMethods_Log2.png")

fig = plt.figure(11)
fig.clf()
ax = fig.add_subplot(111)
cs = ax.pcolorfast(x[:,0], y[0,:], np.log10(flux3/flux0), vmin=vmin, vmax=vmax, cmap='RdBu')
cb = plt.colorbar(cs, ax=ax, label="Log10(Flux/Reference)")
ax.set_title('Contrast: 3')
ax.set_xlabel('x (cm)')
ax.set_ylabel('y (cm)')
ax.set_aspect('equal')
fig.savefig("FluxMethods_Log3.png")

plt.close('all')