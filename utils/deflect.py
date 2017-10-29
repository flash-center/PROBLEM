#!/usr/bin/python2.7
# -*- coding: utf-8 -*-
"""
deflect.py: Provides a description of charged particle deflections in magnetic fields

# TODO: Add electric fields(?) May not be reasonable genetic-algorithm-wise.

Created by Scott Feister on Sat Apr 29 11:55:54 2017
"""

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import os

import numpy as np
from numpy import random
import scipy.constants as sc
from scipy.interpolate import RectBivariateSpline
from scipy import misc

def propZ(Vxyz, xi, yi, zi, zf):
    """ Straight-line propagation of particles by fixed amount in z-dimension based on a velocity vector
    Inputs (for N particles):
       Vxyz:  (N x 3) 2d Numpy Array, Velocity x/y/z vector, in SI units
       xi:     N-element 1d Numpy Array, initial x positions of particles
       yi:     N-element 1d Numpy Array, initial y positions of particles       
       zi:     Float or N-element 1d Numpy Array, initial z position(s) of particles       
       zf:     Float, Final z positions of particles
    Outputs:
       xf:     N-element 1d Numpy Array, final x positions of particles
       yf:     N-element 1d Numpy Array, final y positions of particles
    """
    deltaZ = zf - zi # Distance along Z dimension, in SI units, from initial position to final position
    deltaT = deltaZ/Vxyz[:,2] # Parametrized times of flights from region to screen
    
    xf = xi + Vxyz[:,0] * deltaT # Calculate final x-position of particle
    yf = yi + Vxyz[:,1] * deltaT # Calculate final y-position of particle
    
    return xf, yf

def deflect(Vixyz, Bxyz, mass = sc.m_p, q = sc.e):
    """ Calculate particle deflection through a (relatively) thin sheet of magnetic field
    Inputs (for N particles):
       Vixyz:  (N x 3) 2d Numpy Array, Pre-deflection velocity x/y/z vector, in SI units (meters/sec)
       Bxyz:   (N x 3) 2d Numpy Array, Magnetic field times linear distance, in SI units (Tesla-meters), at locations of particles
       mass:    Float, particle mass (SI unit). Default: Proton mass.
       q:  Float, particle charge (SI unit). Default: Proton charge.
    Outputs:
       Vfxyz:  (N x 3) 2d Numpy Array, Post-deflection velocity x/y/z vector, in SI units (meters/sec)
    """
    
    Vimag = np.sqrt(np.sum(Vixyz**2, axis=-1))
    Bmag = np.sqrt(np.sum(Bxyz**2, axis=-1))
    Vihat = (Vixyz.T/Vimag.T).T # Velocity vectors normalized (v-hat)
    Bhat = (Bxyz.T/Bmag.T).T # Magnetic field * distance vectors normalized (B-hat)
    
    Bperp = np.sqrt(np.sum(np.cross(Bxyz, Vihat)**2, axis=-1)) # Perpendicular magnitudes of Magnetic field * distance
    
    gamma = 1.0/np.sqrt(1 - Vimag**2/sc.c**2) # Relativistic gamma, for each particle
    thetas = (q * Bperp) / (gamma * mass * Vimag) # Deflection angular magnitude
    Vfxyz = Vixyz + (thetas.T * np.cross(Bhat, Vixyz).T).T # Adjust velocity according to the angular deflection. WARNING: Small angle approx. TODO: Check/improve this math.
    #print(np.mean(np.rad2deg(thetas)))
    #print(np.std(np.rad2deg(thetas)))

    return Vfxyz

# Todo: Does this really need its own function?
def phist(xf, yf, bins=20, wid=1.0, plotdir=None, name='phist.png'):
    """ Custom histogram for particles
    Inputs:
        bins:   Number of bins in each dimension
        wid:  Determines histogram2d range (as half the width); same units as xf, yf
        plotdir: (Optional) If None, no plot is generated. If not none, plot is saved into "outdir"
        outname: (Optional) Only appplies if outdir is not "None". In that case, this is the name of the output histogram plot PNG.
    """
    range=((-wid/2.,wid/2.),(-wid/2.,wid/2.))
    H, xedges, yedges = np.histogram2d(xf, yf, bins=bins, range=range)
    
    if plotdir is not None:
        fig = plt.figure(1)
        fig.clear()
        ax = fig.add_subplot(111)
        ax.imshow(H.T, interpolation='nearest', origin='low', extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]])
        fig.savefig(os.path.join(plotdir, name))
        
    return H, xedges, yedges

def mycbeam(srcz, srcrad = 100.0e-6, apdegs = 20, nparts = 10000):
    """ Custom capsule-based beamed source of particles, aimed in the +z direction
    Inputs:
       srcz:   Source (e.g. capsule) center Z-coordinate, in SI units (meters)
       srcrad: Source radius of volume, in SI units (meters)
       apdegs: Beam aperture full-angle, in degrees. Note: apdegs must be less than 180 (half-capsule)
       nparts: Number of particles to place in the beam
    Outputs:
       x:      N-element 1d Numpy Array, x positions of particles, in SI units
       y:      N-element 1d Numpy Array, y positions of particles, in SI units
       z:      N-element 1d Numpy Array, z positions of particles, in SI units
       vx:     N-element 1d Numpy Array, x velocity of particles, in SI units
       vy:     N-element 1d Numpy Array, y velocity of particles, in SI units
       vz:     N-element 1d Numpy Array, z velocity of particles, in SI units
    """
    
    # Define particle position; uniformly distributed within sphere, by building bounding box and using "discard" method
    posxyz = 2*srcrad * (random.rand(3*nparts, 3) - 0.5)  # Create three times more particles than are needed (generously too many), since discard method will be used  

    ct = np.sum(posxyz**2, axis=-1) < srcrad**2
    if np.sum(ct) < nparts:
        raise Exception("Didn't generate enough particles within the sphere. Could just be terrible random luck -- try again.")
        
    posxyz = posxyz[ct,:]
    
    x = posxyz[:nparts,0]
    y = posxyz[:nparts,1]
    z = srcz + posxyz[:nparts,2]
        
    # Define particle ejection theta/phi ejection angles, constraining to the aperture angle
    #KEs_MeV = random.normal(3.0, 0.15, (nparts,))
    #print(KEs_MeV.shape)

    velmag = 0.1*sc.c*random.normal(1.0, 0.01, (nparts,))# Magnitude of velocity vector # TODO: Base this upon gaussian distribution of kinetic energies

    costhmin = np.cos(np.deg2rad(apdegs/2.0))
    costhmax = np.cos(0.0)
    costheta = (costhmax - costhmin) * random.rand(nparts) + costhmin # Cos(theta) varies uniformly from -1 to 1, for 360 degrees. Or, 0.9 to 1.0 for narrow beam, etc.

    vthet = np.arccos(costheta) # Theta of ejection can be between 0 and 2pi
    vphi = 2 * np.pi * random.rand(nparts) # Phi of ejection can be between 0 and apdegs/2.0
    
    vx = velmag * np.sin(vthet) * np.cos(vphi)
    vy = velmag * np.sin(vthet) * np.sin(vphi)
    vz = velmag * np.cos(vthet)
    
    return x, y, z, vx, vy, vz
    
def myradios(xgv, ygv, CBx, xd, yd, z_defl, Vixyz):
    """ Create the custom radiographs from a given magnetic field distribution
    Inputs:
       xgv:    1d list of x values for CBx array at deflection plane, in SI units
       ygv:    1d list of y values for CBx array at deflection plane, in SI units
       CBx     2d uniform grid of Magnetic field * distance values at deflection plane, in SI units
       xd:      N-element 1d Numpy Array, x positions of particles at deflection plane, in SI units
       yd:      N-element 1d Numpy Array, y positions of particles at deflection plane, in SI units
       z_defl:  Float, Z-coordinate of the deflection plane, in SI Units
       Vixyz:  (N x 3) 2d Numpy Array, Pre-deflection velocity x/y/z vector, in SI units (meters/sec)
    Outputs:
        Hs:     List of histograms, which can be quantitatively compared for fitness as desired
    """
    # Interpolate to get dots down onto point
    rbv = RectBivariateSpline(xgv, ygv, CBx)
    Bx = rbv.ev(xd, yd) # Interpolate to the points of interest
    By = np.zeros(Bx.shape)
    Bz = np.zeros(Bx.shape)
    
    Bxyz = np.stack((Bx, By, Bz), axis=-1)
    
    ## Deflect particles within the deflection plane
    Vfxyz = deflect(Vixyz, Bxyz)
    
    ## Propagate from deflection plane to one or various imaging planes
    z_screens = np.arange(20e-3, 70e-3, 20e-3) # Z coordinates of imaging planes, in SI (meters)
    Hs = [None]*len(z_screens)
    for i in range(len(z_screens)):
        xf, yf = propZ(Vfxyz, xd, yd, z_defl, z_screens[i]) # Propagate particles forward to screen
        Hs[i], _, _ = phist(xf, yf, wid=20e-3, bins=10)

    return Hs
    
    
if __name__ == "__main__":
    plotdir = r'/home/sfeister/myouts/scratch'

    ## Create particles at source
    nparts = 500000
    z_src = -1.84e-2 # Z coordinate of capsule center, in SI (meters)
    xi, yi, zi, vxi, vyi, vzi = mycbeam(z_src, apdegs = 50.0, nparts=nparts)
    Vixyz = np.array([vxi, vyi, vzi]).T
    
    # (Optional) Make a slice plot of the source
    ct = np.abs(zi - z_src) < 10.0e-6
    H, xedges, yedges = phist(xi[ct], yi[ct], wid=200.0e-6, bins=100, plotdir=plotdir, name='sourceslc.png')
    
    ## Propagate from source to deflection plane
    z_defl = 0.0 # Z coordinate of deflection plane, in SI (meters)
    xd, yd = propZ(Vixyz, xi, yi, zi, z_defl)
    
    H, xedges, yedges = phist(xd, yd, wid=10e-3, bins=100, plotdir=plotdir, name='defplane.png')
    
    ## Assign bogus magnetic fields at the deflection plane, onto the particles
    #Bxyz = np.array([0.56e-1 * random.rand(nparts), 0.16e-1 * random.rand(nparts), 0*random.rand(nparts)]).T # B*dist fields, in SI units (Tesla-meters)
    
    # Generate fake scalar field data
    #X, Y = np.mgrid[-1.5:1.55:0.02, -1.5:1.55:0.02]
    #C = np.sin(X + Y) + np.sin(2 * X - Y) + np.cos(3 * X + 4 * Y)    
    #X = 2.e-3 * X
    #Y = 2.e-3 * Y
    #C = 1.e-2 * C
    C1 = 0.00005 * misc.face()[:,:,0].astype('float')
    xgv = 5e-3 * float(C1.shape[0])/C1.shape[1] * np.linspace(-0.75, 0.75, C1.shape[0])
    ygv = 5e-3 * np.linspace(-1, 1, C1.shape[1])
    
    # Interpolate to get dots down onto point
    rbv = RectBivariateSpline(xgv, ygv, C1)
    Bx = rbv.ev(xd, yd)
    By = np.zeros(Bx.shape)
    Bz = np.zeros(Bx.shape)
    
    Bxyz = np.stack((Bx, By, Bz), axis=-1)
    
    ## Deflect particles within the deflection plane
    Vfxyz = deflect(Vixyz, Bxyz)
    
    ## Propagate from deflection plane to one or various imaging planes
    for z_screen in np.arange(10e-3, 80e-3, 10e-3):
        #z_screen = 10e-3 # Z coordinate of imaging plane, in SI (meters)
        xf, yf = propZ(Vfxyz, xd, yd, z_defl, z_screen) # Propagate particles forward to screen
        
        H, xedges, yedges = phist(xf, yf, wid=20e-3, bins=50, plotdir=plotdir, name=str(z_screen*1e3) + '.png')
        #H, xedges, yedges = phist(xf, yf, wid=20e-3, bins=50)

