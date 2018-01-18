"""
Main solving algorithm of the
Monge-Ampere screen mapping potential problem using
a finite difference scheme.
"""
import numpy as np
import scipy as sp
import scipy.interpolate
import os.path

def setup(ri, li, rs, dxS, image, mask):
    """
    Takes proton radiography distance parameters, a proton flux image, and a
    flux image selection mask and sets up the deflection-field
    potential problem.

    Args:
        ri (float): Distance from source to plasma (cm).
        li (float): Distance across plasma (cm).
        rs (float): Distance from plasma to screen (cm).
        dxS (float): Pixel size (cm)
        image (array): Proton flux image to analyze.
        mask (array): Bitwise mask of data to analyze in the flux image,
                      assuming 1 is masked and 0 is unmasked.

    Returns:
        plasma_x (array): Plasma x-coordinates (cm).
        plasma_y (array): Plasma x-coordinates (cm).
        flux0 (array): Initial flux image.
        flux (array): Flux image.
    """

    magnify = (rs+ri+.5*li)/(ri+.5*li)
    # Leave the original image unchanged.
    flux_image = np.copy(image)
    flux_mask = np.copy(mask)

    # Mask the flux_image and
    # fill the masked positions with the mean flux.
    #
    # TODO This may add gradient to the proton fluence and therefore
    # create artificial magnetic fields-- may need to find a different way to
    # do this. Possibly extrapolate gradients from the boundary?
    #flux = np.ma.array(flux_image, mask=np.logical_not(flux_mask))
    flux = flux_image * np.logical_not(flux_mask)
    flux_mean = np.nanmean(flux)
    flux = flux + flux_mean * flux_mask

    # Calculate the total flux throughout the selected region
    flux_tot = np.nansum(np.logical_not(flux_mask)*flux)

    # TODO
    # @usualgaussnumber applies a Gaussian filter after this step, see
    # pradreconstructionfluxdistsetup.m

    # Prepare the flux and
    # flux_0 for the primary reconstruction process.
    # flux_0 and flux_0 mask are the initial flux distributions,
    # which we choose to be uniform over the selected region while
    # preserving total flux.
    # Then add a shadow flux region to flux0 which is equal to the flux mean.
    flux0_mask = np.copy(flux_mask)
    flux0 = np.zeros(flux.shape)
    flux0 = (flux_tot * np.logical_not(flux0_mask))/(np.nansum(np.logical_not(flux0_mask)))
    flux0 = flux0 + flux0_mask * np.mean(flux0)

    # Variables:
    # Nx, N1 are the number of cells on the 0th and 1st axis
    # dx is the length of each cell (square)
    # L0, L1 are the lengths of the flux image in the 0th and 1st axis
    Nx, Ny = flux.shape
    Lx = dxS*Nx
    Ly = dxS*Ny

    # Create the coordinates for the flux samples for the centers of the
    # bins on the screen.
    screen_x = np.linspace((-Lx/2)+.5*dxS, (Lx/2)-.5*dxS, Nx)
    screen_y = np.linspace((-Ly/2)+.5*dxS, (Ly/2)-.5*dxS, Ny)
    screen_x, screen_y = np.meshgrid(screen_x, screen_y, indexing='ij')

    # TODO In @usualgaussnumber's code,
    # the flux_bndy was created when he selected a
    # region from the Matlab environment. Here, instead, we will
    # create the coordinates from the mask.
    # I'm not quite sure if this is used beyond simply blurring the
    # edges of the flux, as in pradreconstructionfluxdistsetup.m --@jtlaune

    # Scale the plasma coordinates according to the distances.
    plasma_x = screen_x/magnify
    plasma_y = screen_y/magnify

    return(plasma_x, plasma_y, flux0, flux)

def solve(X, Y, flux0, flux, dt, tol,
          chk=False, interval=1000,
          nan_exception=True, start_chk=False,
          start_phin=None, start_step=None,
          save_dir=None,
          Nstep=50000):
    """
    Main solving algorithm of the
    Monge-Ampere deflection-field potential problem using
    a finite difference scheme.

    Args:
        X (array): Plasma x-coordinates (cm).
        Y (array): Plasma y-coordinates (cm).
        flux (array): Proton flux image.
        flux0 (array): Initial flux image.
        dt (float): Timestep for Euler method.
        tol: Relaxation constant for Euler method.
        chk=False (bool): Save checkpoint files to output text files at a
                          timestep interval specified by interval.
        interval=1000 (int): Interval with which to save output files.
        nan_exception=True (bool): Optional. Set to true if you would like the
                                   algorithm to ignore NaN values.
        start_chk=False (bool): Start the solver from a specified
                                     checkpoint file.
        start_phin=None (str): 'phin' checkpoint file to use. Must be specified
                               if save_chk=True.
        start_step=None (int): Starting step of the 'phin' file. Must be
                               specified if save_chk=True.
        save_dir=None (str): If chk=True, this specifies the directory to save
                             the checkpoint files to. If this is left as None,
                             the default is the current directory.
        Nstep=50000 (int): Maximum number of steps before stopping.

    Returns:
        phin (numpy.ndarray): Screen mapping potential in the plasma
                              coordinates.
        phix (numpy.ndarray): x-gradient of the screen mapping potential in
                              the plasma coordinates.
        phiy (numpy.ndarray): x-gradient of the screen mapping potential in
                              the plasma coordinates.
    """
    plasma_x = np.copy(X)
    plasma_y = np.copy(Y)
    flux0 = np.copy(flux0)
    flux = np.copy(flux)
    flux_mean = np.mean(flux)

    if start_chk:
        start_phin = os.path.abspath(start_phin)
    if save_dir is not None:
        save_dir = os.path.abspath(save_dir)

    # Make sure that the files/directories specified are valid.
    if start_chk and (start_phin == None or start_step == None):
        raise Exception('Please specify a checkpoint phin file and '
                        'starting step.')
    if start_chk:
        if not os.path.isfile(start_phin):
            raise Exception('Please specify a valid starting checkpoint file '
                            'for phin.')
        if type(start_step) != int:
            raise Exception('Please provide a valid starting step.')

    if save_dir is not None:
        if not os.path.isdir(save_dir):
            raise Exception('Please specify a valid save directory')

    ############################
    # BEGINNING OF MAIN SOLVER #
    ############################
    if start_chk:
        # If we are starting from a checkpoint file then we load it into phin to
        # iterate.
        phin = np.loadtxt(start_phin, delimiter=',')
    if not start_chk:
        # Set up plasma parameters.
        # Initial coordinate perturbation.
        phin = .5 * (plasma_x**2 + plasma_y**2)

    # Read the size of the grid
    plasma_Nx, plasma_Ny = plasma_x.shape

    # Create plasma parameters
    dx = plasma_x[1,0] - plasma_x[0,0]
    dy = plasma_y[0,1] - plasma_y[0,0]

    plasma_Lx = plasma_x[-1,0] - plasma_x[0,0]
    plasma_Ly = plasma_y[0,-1] - plasma_y[0,0]

    # Create the derivative arrays phix, phiy, phixx, phiyy, phixy.
    # They will be the same size as phin and will be replaced every
    # timestep.
    phix = np.zeros((plasma_Nx, plasma_Ny))
    phiy = np.zeros((plasma_Nx, plasma_Ny))
    phixx = np.zeros((plasma_Nx, plasma_Ny))
    phiyy = np.zeros((plasma_Nx, plasma_Ny))
    phixy = np.zeros((plasma_Nx, plasma_Ny))

    RMSfluxerror = np.zeros((Nstep))
    pointfluxerror = np.zeros((Nstep))
    Nxrand = int(np.floor(np.random.ranf()*(plasma_Nx+1)))
    Nyrand = int(np.floor(np.random.ranf()*(plasma_Ny+1)))

    if chk:
        if save_dir==None:
            # Save the files to the current directory.
            np.savetxt(os.path.join(os.getcwd(), 'plasma_x.txt'),
                       plasma_x, delimiter=',')
            np.savetxt(os.path.join(os.getcwd(), 'plasma_y.txt'),
                       plasma_y, delimiter=',')
        else:
            # Save the files to the specified save directory.
            np.savetxt(os.path.join(save_dir, 'plasma_x.txt'),
                       plasma_x, delimiter=',')
            np.savetxt(os.path.join(save_dir, 'plasma_y.txt'),
                       plasma_y, delimiter=',')

    # Create the interpolator function.
    flux_interp = sp.interpolate.RegularGridInterpolator(
                                    (plasma_x[:,0], plasma_y[0,:]),
                                    flux,
                                    method='linear',
                                    bounds_error=False,
                                    fill_value=flux_mean)

    ###################################
    # Begin finite difference scheme. #
    ###################################
    #
    # loop Nstep times
    # 1. Check for instabilities in phin
    # 2. Create arrays for phix, phiy, phixx, phiyy, phixy
    #    - These are the derivative functions
    # 3. Finite difference scheme in 9 regions:
    #    - corners (4), edges(4), middle
    #    - get the first and second derivatives in each of these regions
    # 4. Calculate quantities and update phin
    #

    for timestep in range(Nstep):
        if start_chk:
            # If we are starting from a checkpoint file then we want to make
            # sure the timestep is accurate. Note that this will overwrite any
            # other files made after the checkpoint file we are using.
            timestep = start_step + timestep + 1

        print(timestep)

        # Test whether algorithm is stable
        if nan_exception and np.any(np.isnan(phin)):
            raise Exception('The reconstruction algorithm is unstable.')

        # Create gradient arrays phix, phiy
        phix, phiy = np.gradient(phin, dx, dy)

        # Begin adjusting gradient arrays for Neumann BCs.

        # phix corners.
        #phix[0,0] = plasma_x[0,0]
        #phix[0,plasma_Ny-1] = plasma_x[0,plasma_Ny-1]
        #phix[plasma_Nx-1,0] = plasma_x[plasma_Nx-1,0]
        #phix[plasma_Nx-1,plasma_Ny-1] = plasma_x[plasma_Nx-1,plasma_Ny-1]
        #
        ## phiy corners.
        #phiy[0,0] = plasma_y[0,0]
        #phiy[0,plasma_Ny-1] = plasma_y[0,plasma_Ny-1]
        #phiy[plasma_Nx-1,0]  = plasma_y[plasma_Nx-1,0]
        #phiy[plasma_Nx-1,plasma_Ny-1] = plasma_y[plasma_Nx-1,plasma_Ny-1]

        # phix edges.
        phix[0,:] = plasma_x[0,:]
        phix[plasma_Nx-1,:] = plasma_x[plasma_Nx-1,:]

        # phiy edges.
        phiy[:,0] = plasma_y[:,0]
        phiy[:,plasma_Ny-1] = plasma_y[:,plasma_Ny-1]

        # Lets try to just match up the coordinates. Seems to be a problem on
        # the boundary with the non neumann BC sides. --@jtlaune
        #phix[0,:] = plasma_x[0,:]
        #phix[-1,:] = plasma_x[-1,:]
        #phix[:,0] = plasma_x[:,0]
        #phix[:,-1] = plasma_x[:,-1]
        #
        #phiy[0,:] = plasma_y[0,:]
        #phiy[-1,:] = plasma_y[-1,:]
        #phiy[:,0] = plasma_y[:,0]
        #phiy[:,-1] = plasma_y[:,-1]

        # Create gradient arrays phixx, phiyy, phixy after we have set
        # BCs for phix and phiy.
        phixx, phixy = np.gradient(phix, dx, dy)
        _, phiyy = np.gradient(phiy, dx, dy)

        # Boundary conditions for phixy.
        # TODO Not sure what this is?? -- @jtlaune
        phixy[0,:] = 0
        phixy[plasma_Nx-1,:] = 0
        phixy[:,0] = 0
        phixy[:,plasma_Ny-1] = 0

        # Now calculate the derived quantities.
        det = phixx*phiyy - phixy**2

        # Interpolate the flux at initial points.
        # We fill points outside of the grid with flux_mean since that is what
        # we are assuming the protons are doing outside of the interaction
        # region.
        fluxn = flux_interp((phix,phiy))

        # Calculate derivative quantity.
        # TODO @usualgaussnumber's code had absolute det? Not sure why. This is
        # not mentioned in the Sulman numerical solution paper.
        Fn = np.log(fluxn*np.absolute(det)/flux0)


        # Create new phi_derived array to do error checking.
        phi_derived = phin+dt*Fn

        # Calculate the flux error.
        pointfluxerror[timestep] = (phi_derived[Nxrand,Nyrand]
                                    - phin[Nxrand,Nyrand])

        RMSfluxerror[timestep] = np.sqrt(np.mean((phi_derived-phin)**2))

        # Finally update phin for the next timestep.
        phin = phi_derived

        # TODO figure out a better way to do the tolerance checking for steady
        # state solution.
        if np.linalg.norm(Fn) < tol:
            break

        if chk:
            if (timestep % interval) == 0:
                if save_dir==None:
                    # Save the files to the current directory.
                    np.savetxt('phin' + str(timestep) + '.txt',
                               phin, delimiter=',')
                    np.savetxt('phix' + str(timestep) + '.txt',
                               phix, delimiter=',')
                    np.savetxt('phiy' + str(timestep) + '.txt',
                               phiy, delimiter=',')
                else:
                    # Save the files to the specified save directory.
                    np.savetxt(save_dir+'/phin' + str(timestep) + '.txt',
                               phin, delimiter=',')
                    np.savetxt(save_dir+'/phix' + str(timestep) + '.txt',
                               phix, delimiter=',')
                    np.savetxt(save_dir+'/phiy' + str(timestep) + '.txt',
                               phiy, delimiter=',')

        #################################
        # End finite difference scheme. #
        #################################

    return(phin, phix, phiy)
