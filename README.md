# PROBLEM

A proton radiography reconstruction tool in Python.

This provides a library of functions to reconstruct the path integrated magnetic
field present in a plasma based on a proton radiography flux image.

## Disclaimer


## Installation

Requirements:

* [`numpy`](http://www.numpy.org/)
* [`scipy`](https://www.scipy.org/)

```bash
pip install numpy scipy
git clone https://github.com/flash-center/PROBLEM.git
cd problem
python setup.py install
```

## Usage

In order to reconstruct a magnetic field from a flux image, the following data
are needed:

* the relevant distances for the problem, i.e. the distance from the proton
    source to the plasma, the distance across the plasma, and the distance from
    the plasma to the detector
* the detector pixel size
* the flux image
* the velocity of the protons

The standard process of reconstructing the path-integrated magnetic fields from
the flux image is as follows:

* set up the problem parameters with `problem.screen.setup()`
* solve for the screen mapping potential with `problem.screen.solve()`
* reconstruct the perpendicular deflection fields from the screen mapping
    potential with `problem.deflect.reconstruct()`
* reconstruct the path-integrated magnetic fields in the x & y directions using
    `problem.deflect.magpath()`

### Example Problem

```python
import problem.screen
import problem.deflect

ri = 1 # Distance from proton source to plasma.
li = .1 # Distance across plasma.
rs = 20 # Distance from plasma to screen.
dxS =.02 # Pixel size of detector screen.
v = 5.24e9 # Proton velocity from source.
flux = np.loadtxt('fluximage.out') # Load the flux image from a file.
# Load the mask from a file.
# The mask is a bitwise mask that will be used to select a subset of the flux 
# image for reconstruction purposes. 1's are masked points, 0's are unmasked.
mask = np.loadtxt('mask.out') 

# Set up the initial conditions for the reconstruction problem.
plasma_x, plasma_y, flux0, flux = problem.screen.setup(ri, li, rs, dxS, image, mask)

dt = 1e-8 # Timestep for solver.
tol = 1e-3 # Tolerance level for steady state solution.

# Initiate the solver.
# chk=True means the solver will save checkpoint files at the desired interval
# in the directory specified by save_dir.
phin, phix, phiy = problem.screen.solve(plasma_x, plasma_y, flux0, flux,
                                     dt, tol, chk=True, interval=1000,
                                     save_dir='solve/')

# Reconstruct the perpendicular deflection fields.
wBx, wBy = problem.deflect.reconstruct(ri, li, rs, v, 
                                    plasma_x, plasma_y, 
                                    phix, phiy)

# Reconstruct the path-integrated magnetic fields in x & y directions.
Bxpath, Bypath = problem.deflect.magpath(wBx, wBy)
```

## Future Work

Future goals include

* perform more validation and test cases
* create user-friendly command line tools
* create compatibility with reader functions from 
    [`pradreader`](https://github.com/jtlaune/problemreader)
