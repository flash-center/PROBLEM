import pradreader
import argparse
import numpy as np
import os
import problem.screen
import problem.adaptive_fixed_point

def get_input():
    parser = argparse.ArgumentParser(
                description="This script is used to reconstruct proton"
                            "radiography data from a PRadReader data file.")

    parser.add_argument("input_file",
                        action="store", type=str,
                        help="Input file 1.")

    parser.add_argument("--tol", default=1e-3,
                        action="store", type=float,
                        help="Solver tolerance. Default: 1e-3")

    parser.add_argument("--step", default=1e-8,
                        action="store", type=float,
                        help="Solver step magnitude. Default: 1e-8.")

    parser.add_argument("--nstep", default=400000,
                        action="store", type=int,
                        help="Maximum number of solver steps. Default: 400000")

    parser.add_argument("--xstep",
                        action="store", type=float,
                        help="X-coordinate step size.")

    parser.add_argument("--ystep",
                        action="store", type=float,
                        help="Y-coordinate step size.")

    parser.add_argument("--chk", default=1000,
                        action="store", type=int,
                        help="Checkpoint interval (number of steps). " \
                             "Default: 1000 steps")

    parser.add_argument("--save_dir", default="solve",
                        action="store", type=str,
                        help="Directory to save checkpoint files. " \
                             "Default: solve/")

    parser.add_argument("--method", default="PMA",
                        action="store", type=str,
                        help="Method to use. " \
                             "Choose from PMA (parabolic Monge-Ampere) or AFP (adaptive fixed-point). " \
                             "Default: PMA")

    parser.add_argument("-p", "--plot",
                        action="store_true",
                        help="Use matplotlib to plot the results at each checkpoint? " \
                             "Default: False/")

    args = parser.parse_args()

    return(args)

def setup_input():
    args = get_input()
    tol = args.tol
    step = args.step
    nstep = args.nstep
    chk_interval = args.chk
    save_dir = args.save_dir
    method = args.method
    plot = args.plot
    if args.xstep is None:
        dx = float(input("Please specify the x-coordinate step size (in cm): "))
    else:
        dx = args.xstep
    if args.ystep is None:
        dy = float(input("Please specify the y-coordinate step size (in cm): "))
    else:
        dy = args.ystep
    prad = pradreader.reader.loadPRR(ifile=args.input_file)
    prad.show()
    return prad, tol, step, nstep, dx, dy, chk_interval, save_dir, method, plot

def calc_coord(dx, dy, shape):
    lenx = shape[0]
    leny = shape[1]
    totx = (lenx-1)*dx
    toty = (leny-1)*dy
    print(totx)
    x = np.linspace(-totx/2, totx/2, num=lenx, endpoint=True)
    y = np.linspace(-toty/2, toty/2, num=leny, endpoint=True)
    x, y = np.meshgrid(x, y, indexing="ij")
    return(x, y)

def ensure_save_dir(save_dir):
    dirpath = os.path.abspath(save_dir)
    if not os.path.exists(dirpath):
        os.makedirs(dirpath)

def reconstruct():
    prad, tol, step, nstep, dx, dy, chk_interval, save_dir, method, plot = setup_input()
    flux = prad.flux2D
    flux0 = prad.flux2D_ref
    x, y = calc_coord(dx, dy, flux.shape)
    ensure_save_dir(save_dir)
    print(x,y)
    chk = True
    if method == "PMA":
        phin, phix, phiy = problem.screen.solve(x, y, flux0, flux,
                                                step, tol, Nstep=nstep,
                                                chk=bool(chk_interval),
                                                interval=chk_interval,
                                                nan_exception=True,
                                                save_dir=save_dir,
                                                plot=plot)
    elif method == "AFP":
        AFP = problem.adaptive_fixed_point.AFP(x[:,0], y[0,:])
        phin, phix, phiy, J = AFP.AFP(flux, flux0, nstep,
                save_interval=chk_interval,
                save_dir=save_dir, plot=plot)
    else:
        raise NotImplementedError(f'Method {method} not known')


if __name__=="__main__":
    reconstruct()

