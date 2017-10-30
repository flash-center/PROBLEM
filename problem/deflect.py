"""
Takes the gradients of the solution to the screen mapping potential problem and
reconstructs the perpendicular deflection field.
"""

import numpy as np
import scipy as sp
import scipy.interpolate
import scipy.misc
import scipy.ndimage

from .constants import M_PROTON_G, ESU, C_CMS

def reconstruct(ri, li, rs, v, x, y, phix, phiy):
    """
    Takes x, y gradients to the solution to screen mapping potential problem and
    reconstructs the perpendicular deflection fields wBx and wBy.

    Args:
        ri (float): Distance from source to plasma (cm).
        li (float): Distance across plasma (cm).
        rs (float): Distance from plasma to screen (cm).
        v (float): Velocity of protons (cm/s).
        x (array): Plasma x-coordinates (cm). 
        y (array): Plasma x-coordinates (cm).
        phix (array): Gradient of screen mapping potential in x-direction.
        phiy (array): Gradient of screen mapping potential in y-direction.

    Returns:
        wBx (array)
        
    """
    # TODO Add in option for masking the path-int B field.
    
    # Input variables.
    magnify = (rs + ri + .5*li)/(ri+.5*li)
    map_pot_x = np.copy(phix)
    map_pot_y = np.copy(phiy)
    plasma_x = np.copy(x)
    plasma_y = np.copy(y)
    
    # We multiply the whole expression by magnify to put the perp-deflection
    # fields into screen coordinates.
    wBx = magnify*(v/rs)*(map_pot_x - plasma_x)
    wBy = magnify*(v/rs)*(map_pot_y - plasma_y)
    
    return(wBx, wBy)

def magpath(wBx, wBy):
    """
    Takes the perpendicular deflection field and reconstructs the path
    integrated magnetic field.

    Args:
        wBx (array): x-component perpendicular deflection field.
        wBy (array): y-component perpendicular deflection field.

    Returns:
        Bxpath (array): Path integrated magnetic field x-component. 
        Bypath (array): Path integrated magnetic field y-component.
    """
    
    Bxpath = -(M_PROTON_G*C_CMS/ESU)*wBy
    Bypath = (M_PROTON_G*C_CMS/ESU)*wBx


    return(Bxpath, Bypath)

def fluximage(ri, li, rs, v, x, y, N, wBx, wBy):
    """
    Creates a flux image out of a perpendicular deflection field. 

    Args:
        ri:
        li:
        rs:
        v:
        x (array): Perpendicular deflection field x-coordinates.
        y (array): Perpendicular deflection field y-coordinates.
        wBx (array): Perpendicular deflection field x-component.
        wBy (array): Perpendicular deflection field y-component.

    Returns:
        flux_image (array): Generated flux image.
    """
    # TODO Maybe change this to act on the reference flux.
    magnify = (rs+ri+.5*li)/(ri+.5*li)
    
    print('Creating interpolator functions...')
    
    #fx = sp.interpolate.RegularGridInterpolator((x[:,0],y[0,:]),x,
    #                                            bounds_error=False)
    #fy = sp.interpolate.RegularGridInterpolator((x[:,0],y[0,:]),y,
    #                                            bounds_error=False)
    fwBx = sp.interpolate.RegularGridInterpolator((x[:,0],y[0,:]),wBx,
                                                  bounds_error=False)
    fwBy = sp.interpolate.RegularGridInterpolator((x[:,0],y[0,:]),wBy,
                                                  bounds_error=False)
    
    print('DONE')

    prot_num = int(np.sqrt(N))
    dx = x[1,0] - x[0,0]
    dy = y[0,1] - y[0,0]
    # Need to fix this-- cuts off some of the protons when moving to the centers
    # of the bins.
    samp_x = np.linspace(x[0,0]+.5*dx, x[-1,0]-.5*dx, num=prot_num)
    samp_y = np.linspace(y[0,0]+.5*dy, y[0,-1]-.5*dy, num=prot_num)
    samp_x, samp_y = np.meshgrid(samp_x, samp_y, indexing='ij')
    
    print('Interpolating proton deflections...')
    
    # The sampling of the coordinates is useless.
    #samp_x = fx((samp_x, samp_y))
    #samp_y = fy((samp_x, samp_y))
    samp_wBx = fwBx((samp_x, samp_y))
    samp_wBy = fwBy((samp_x, samp_y))
    
    print('DONE')

    screen_x = magnify*samp_x + (rs/v)*samp_wBx
    screen_y = magnify*samp_y + (rs/v)*samp_wBy

    print('Histogramming protons...')

    flux_image = np.histogram2d(screen_x.ravel(), screen_y.ravel(),bins=x.shape)
    
    print('DONE')
    
    return(flux_image[0])


def fluximage2(x, y, phix, phiy, flux0, scale_fact=1, scale_order=3):
    """
    An alternative approach to creating a flux image out of a perpendicular deflection field. 
    
    Args:
        x (array): Plasma x-coordinates (cm). 
        y (array): Plasma x-coordinates (cm).
        phix (array): Gradient of screen mapping potential in x-direction.
        phiy (array): Gradient of screen mapping potential in y-direction.
        scale_fact: Integer factor by which to upscale arrays before analysis; a larger number slows the algorithm but fills out low-flux regions better
        scale_order: Order of the spline interpolation for scipy.ndimage.zoom
    Returns:
        flux_image (array): Generated flux image.
    """    
    
    xgv = x[:,0].flatten()
    ygv = y[0,:].flatten()
    
    if scale_fact != 1:
        print("Rescaling...")
        xgv = scipy.ndimage.zoom(xgv, scale_fact, order=scale_order)
        ygv = scipy.ndimage.zoom(ygv, scale_fact, order=scale_order)
        phix = scipy.ndimage.zoom(phix, scale_fact, order=scale_order)
        phiy = scipy.ndimage.zoom(phiy, scale_fact, order=scale_order)
        flux0 = scipy.ndimage.zoom(flux0, scale_fact, order=scale_order)
        
    dx = np.mean(np.diff(xgv))
    dy = np.mean(np.diff(ygv))
    x_edges = np.append(xgv - dx/2.0, xgv[-1] + dx/2.0)
    y_edges = np.append(ygv - dy/2.0, ygv[-1] + dy/2.0)
    
    print('Performing histogram...')

    flux_image, _, _ = np.histogram2d(phix.flatten(), phiy.flatten(), bins=[x_edges, y_edges], weights=flux0.flatten())
    
    if scale_fact != 1:
        print("Descaling...")
        flux_image = scipy.misc.imresize(flux_image, 1./scale_fact, mode='F')

    print('DONE')
    
    return(flux_image)

def fluximage3(ri, li, rs, v, x, y, N, wBx, wBy, Ntest):
    """
    A Monte Carlo approach to creating a flux image out of a perpendicular deflection field. 
    
    Args:
        ri:
        li:
        rs:
        v:
        N: Number of protons in reality
        x (array): Perpendicular deflection field x-coordinates.
        y (array): Perpendicular deflection field y-coordinates.
        wBx (array): Perpendicular deflection field x-component.
        wBy (array): Perpendicular deflection field y-component.
        Ntest: Number of test protons (Monte Carlo)

    Returns:
        flux_image (array): Generated flux image.
    """    
    
   # magnify = (rs + ri + li)/(ri)
    magnify = (rs+li+ri)/(ri+.5*li)

    xgv = x[:,0].flatten()
    ygv = y[0,:].flatten()
    xmin = np.min(xgv)
    xmax = np.max(xgv)
    ymin = np.min(ygv)
    ymax = np.max(ygv)

    dx = np.mean(np.diff(xgv))
    dy = np.mean(np.diff(ygv))
    x_edges = np.append(xgv - dx/2.0, xgv[-1] + dx/2.0)
    y_edges = np.append(ygv - dy/2.0, ygv[-1] + dy/2.0)
    
    # xd:      N-element 1d Numpy Array, x positions of particles at deflection plane, in SI units
    # yd:      N-element 1d Numpy Array, y positions of particles at deflection plane, in SI units
    xd = np.random.uniform(xmin, xmax, size=(Ntest,))
    yd = np.random.uniform(ymin, ymax, size=(Ntest,))
    
    xyd = np.stack((xd, yd), axis=1)
    #del xd, yd
    
    #wBx_rbv = sp.interpolate.RectBivariateSpline(xgv, ygv, wBx)
    #wBy_rbv = sp.interpolate.RectBivariateSpline(xgv, ygv, wBy)
    #wBxd = wBx_rbv.ev(xd, yd)
    #wByd = wBy_rbv.ev(xd, yd)
    
    wBxd = sp.interpolate.interpn((xgv, ygv), wBx, xyd, method='linear')
    wByd = sp.interpolate.interpn((xgv, ygv), wBy, xyd, method='linear')

    xfd = xd + rs/(magnify*v) * wBxd
    yfd = yd + rs/(magnify*v) * wByd
        
    print("Histogramming reference...")
    flux_ref, _, _ = np.histogram2d(xd, yd, bins=[x_edges, y_edges])
    flux_ref = flux_ref * N/Ntest
    
    print("Histogramming signal...")
    flux_image, _, _ = np.histogram2d(xfd, yfd, bins=[x_edges, y_edges])
    flux_image = flux_image * N/Ntest

    print('DONE')
    
    return(flux_image, flux_ref)


