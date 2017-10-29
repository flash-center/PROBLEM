#!/usr/bin/python2.7
# -*- coding: utf-8 -*-
"""
testmaps.py: description

Test magnetic field maps

Created by Scott Feister on Thu May 04 16:24:35 2017
"""

from scipy.ndimage.filters import gaussian_filter
import matplotlib.pyplot as plt
import numpy as np

def myCircle(radius = 5, mapwid = 20, nx = 1000, ny = 1001, xcent = 0.0, ycent = 0.0, radblur=1.0):
    """ Example radial field function, centered on x = y = 0
    Inputs:
        mapwid       Float, width of the entire map (e.g. in meters)
        radius       Float, radius of the circle (e.g. in meters)
        nx, ny       Floats, Number of pixels along each dimension
        xcent, ycent Floats, x and y values for the center of the circle (does not affect the map center)
        radblur      Radius of gaussian blur (e.g. in meters). If set to None, no blur is applied.
    Outputs:
        X         2D NumPy array containing X values (nx x ny)
        Y         2D NumPy array containing Y values (nx x ny)
        C         2D NumPy array containing scalar field values, between 0 and 1.0 (nx x ny)
    Example usage:
        # Make a circle with radius of 5 cm, total screen width of 50 cm, centered on (x, y) = (-15 cm, 10 cm)
        X, Y, C = myCircle(radius = 5, mapwid = 50, radblur=1, xcent=-15, ycent=10)
        Bx = 1.0e5 * C # Rescale map C from 0 to 1, changed to 0 to 1.0e5 G*cm
        plt.pcolormesh(X, Y, Bx)
    """
    xmin = ymin = -mapwid/2.
    xmax = ymax = mapwid/2.
    
    xgv = np.linspace(xmin, xmax, nx)
    ygv = np.linspace(ymin, ymax, ny)
    X, Y = np.meshgrid(xgv, ygv)
    C = np.zeros(X.shape)
    
    R = np.sqrt((X - xcent)**2 + (Y - ycent)**2)
    C[R < radius] = 1.0
    
    if radblur is not None:
        C = gaussian_filter(C, float(radblur)/mapwid * np.array([ny, nx]))
    
    return X, Y, C
    
if __name__ == "__main__":
    
    X, Y, Bx = myCircle(radius = 5, mapwid = 50, radblur=1, xcent=-15, ycent=10)
    
    fig, ax = plt.subplots()
    im = ax.pcolormesh(X, Y, Bx, vmin=0.0, vmax=1.0)
    ax.set_aspect('equal')
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_title("Example circle map plot")
    plt.colorbar(im, label="Scalar value, B$_x$")

