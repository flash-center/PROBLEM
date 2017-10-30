import numpy as np
import matplotlib.pyplot as plt
import problem.deflect

# Load the plasma coordinates.
X = np.loadtxt('x.dat', delimiter=',')
Y = np.loadtxt('y.dat', delimiter=',')

# Checkpoint interval.
interval = 1000

# Plasma parameters.
ri = 1. # Distance from proton source to plasma.
li = 0.1 # Distance across plasma.
rs = 30. # Distance from plasma to screen.
v = 2.51e9 # Velocity of protons.

# Plot all of the solve steps for this problem, saving them as
# images/magBpath######.png.
for i in range(400):
    phix = np.loadtxt('solve/phix'+str(interval*i)+'.txt', delimiter=',')
    phiy = np.loadtxt('solve/phiy'+str(interval*i)+'.txt', delimiter=',')
    phin = np.loadtxt('solve/phin'+str(interval*i)+'.txt', delimiter=',')
    wBx, wBy = problem.deflect.reconstruct(ri, li, rs, v, X, Y, phix, phiy)
    Bxpath, Bypath = problem.deflect.magpath(wBx, wBy)
    
    fig, ax = plt.subplots()
    
    ticks = np.linspace(0,5500,12)
    magB = np.sqrt(Bxpath**2 + Bypath**2)
    cs = ax.pcolorfast(X[:,0], Y[0,:], magB, vmin=0, vmax=5500)
    ax.set_xlabel('x (cm)')
    ax.set_ylabel('y (cm)')
    cb = plt.colorbar(cs, ax=ax,ticks=ticks)

    ax.set_title('Calculated Magnitude Path-Int Magnetic Field (G*cm) Step: '
                 '{}'.format(interval*i))

    fig.savefig('images/magBpath'+str(interval*i).zfill(6)+'.png')
    print('Step: {}'.format(i))
