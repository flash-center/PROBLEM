import numpy as np
import scipy.sparse as sp

def first_derivatives(x, y):
    X, Y = np.meshgrid(x, y, indexing='ij')
    N_x = x.size
    N_y = y.size
    dx = np.mean(np.diff(x))
    dy = np.mean(np.diff(y))
    mult = np.ones((N_x * N_y, 1))
    
    elements_x = mult * [-1, 1]
    elements_x[X.flat == x[ 0], :] = 0
    elements_x[X.flat == x[-1], :] = 0

    elements_y = mult * [-1, 1]
    elements_y[Y.flat == y[ 0], :] = 0
    elements_y[Y.flat == y[-1], :] = 0

    Jx = sp.spdiags(elements_x.T, [N_y, -N_y], N_x*N_y, N_x*N_y).T / (dx*2)
    Jy = sp.spdiags(elements_y.T, [1, -1], N_x*N_y, N_x*N_y).T / (dy*2)

    return Jx, Jy

def second_derivatives(x, y, mult_xx=None, mult_yy=None, mult_xy=None):
    X, Y = np.meshgrid(x, y, indexing='ij')
    N_x = x.size
    N_y = y.size
    dx = np.mean(np.diff(x))
    dy = np.mean(np.diff(y))
    if mult_xx is None: mult_xx = np.ones(N_x * N_y)
    if mult_yy is None: mult_yy = np.ones(N_x * N_y)
    if mult_xy is None: mult_xy = np.ones(N_x * N_y)

    elements_xx = mult_xx[:, np.newaxis] * [1, -2, 1]
    elements_xx[X.flat==x[ 0], :] = mult_xx[X.flat==x[ 0], np.newaxis] * [0, -2, 2]
    elements_xx[X.flat==x[-1], :] = mult_xx[X.flat==x[-1], np.newaxis] * [2, -2, 0]

    elements_yy = mult_yy[:, np.newaxis] * [1, -2, 1]
    elements_yy[Y.flat==y[ 0], :] = mult_yy[Y.flat==y[ 0], np.newaxis] * [0, -2, 2]
    elements_yy[Y.flat==y[-1], :] = mult_yy[Y.flat==y[-1], np.newaxis] * [2, -2, 0]

    elements_xy = mult_xy[:, np.newaxis] * [1, -1, -1, 1]
    elements_xy[X.flat==x[ 0], :] = 0
    elements_xy[Y.flat==y[ 0], :] = 0
    elements_xy[X.flat==x[-1], :] = 0
    elements_xy[Y.flat==y[-1], :] = 0

    Jxx = sp.spdiags(elements_xx.T, [N_y, 0, -N_y], N_x*N_y, N_x*N_y).T / dx**2
    Jyy = sp.spdiags(elements_yy.T, [1, 0, -1], N_x*N_y, N_x*N_y).T / dy**2
    Jxy = sp.spdiags(elements_xy.T, [N_y + 1, N_y - 1, 1 - N_y, -1 - N_y], N_x*N_y, N_x*N_y).T / (4*dx*dy)

    return Jxx, Jyy, Jxy
