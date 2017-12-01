import problem.screen
import numpy as np

# Load the appropriate flux image.
flux = np.loadtxt('flux.dat', delimiter=',')

# Load the plasma coordinates.
plasma_x = np.loadtxt('x.dat', delimiter=',')
plasma_y = np.loadtxt('y.dat', delimiter=',')

# Create an initial flux (i.e., the flux image had there been no deflections).
# In this case, it is appropriate to approximate the initial flux as merely a
# uniform distribution of the mean flux.
flux0 = np.zeros(flux.shape) + np.nanmean(flux)

# We pick a solver step of 1e-8, a steady state relaxation constant of 1e-3, and
# a maximum number of steps as 400,000.
solvestep = 1e-8
relax = 1e-3
Nstep = 400000

# Run the solving algorithm.
# `chk=True` tells the solver to save checkpoint files at the specified
# interval.
# `interval=1000` tells the solver that the checkpoint file interval is 1000.
# `nan_exception=True` tells the program to halt if there are any NaN values
# encountered in the solving algorithm.
# `save_dir=solve` tells the solver to save checkpoint files to the `solve/`
# directory.
phin, phix, phiy = problem.screen.solve(plasma_x, plasma_y,
                                        flux0, flux, solvestep, relax,
                                        Nstep=Nstep,
                                        chk=True, interval=1000,
                                        nan_exception=True,
                                        save_dir='solve')
