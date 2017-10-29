import deflect
import numpy as np
import scipy.interpolate
"""
Utility functions to create proton radiography images from magnetic fields.
"""
def rad_array(nprot, 
              B_x, B_y, B_z,
              B_x_coord, B_y_coord, 
              d1, d2):
    """
    Create a radiography image from a magnetic field and distances from
    capsule to deflection point to screen.

    Args:
        nprot (int): Number of protons.
        field (array): 2D magnetic field array.
        d1 (float): Distance from proton source to magnetic field deflection
            point in SI units.
        d2 (float): Distance from magnetic field deflection to screen in 
            SI units.
    
    Returns:
        rad_arr (array): Binned proton positions on the screen.
    """
    # We assume the magnetic field is set at 0.
    z_defl = 0
    # Source to magnetic field.
    src_to_mag = -d1
    # Magnetic field to screen.
    mag_to_screen = d2
    # Initialize particles at the source.
    xi, yi, zi, vxi, vyi, vzi = deflect.mycbeam(src_to_mag, 
                                                apdegs=50.0,
                                                nparts=nprot)
    # Create array of proton velocities.
    Vixyz = np.array([vxi, vyi, vzi]).T
    
    # Propogate the protons from the source to the magnetic field.
    xd, yd = deflect.propZ(Vixyz, xi, yi, zi, z_defl)
    
    # Create interpolating functions to fit the magnetic field to
    # the proton positions.
    interp_x = scipy.interpolate.RectBivariateSpline(B_x_coord, B_y_coord, B_x)
    interp_y = scipy.interpolate.RectBivariateSpline(B_x_coord, B_y_coord, B_y)
    interp_z = scipy.interpolate.RectBivariateSpline(B_x_coord, B_y_coord, B_z)

    # Evaluate these functions at the proton positions.
    B_x = interp_x.ev(xd, yd)
    B_y = interp_y.ev(xd, yd)
    B_z = interp_z.ev(xd, yd)
    B_xyz = np.stack((B_x, B_y, B_z), axis=-1)

    # Vixyz hasn't changed after the propogation since velocities don't change.
    # Now we will deflect the velocities.
    Vfxyz = deflect.deflect(Vixyz, B_xyz)

    # Using these final velocities, we propogate to the screen from the
    # deflection plane.
    xf, yf = deflect.propZ(Vfxyz, xd, yd, z_defl, mag_to_screen)

    H, xedges, yedges = deflect.phist(xf, yf, wid=20e-3, bins=50,
                                      plotdir='/home/jtlaune/Desktop',
                                      name='histogram.png')
