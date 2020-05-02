import numpy as np

import scipy as sp
import scipy.interpolate as interp
import scipy.sparse as sps
import scipy.sparse.linalg as la

import problem.sparse_jacobians as spj
import problem.plotting as pplot

class AFP:
    def __init__(self, x, y, epsilon=1e-5):
        self.x = x
        self.y = y
        self.X0, self.Y0 = np.meshgrid(x, y, indexing='ij')
        self.shape_2d = self.X0.shape
        self.shape_1d = (self.X0.size, 1)
        self.jac_x, self.jac_y = spj.first_derivatives(x, y)
        self.jac_xx, self.jac_yy, self.jac_xy = spj.second_derivatives(x, y)
        self.epsilon = epsilon
        A = sps.eye(self.X0.size, format='lil')
        B = A[1:, :]
        B[:, 0] = -1
        self.B = B.tocsc()
        self.C = A[:, 1:].tocsc()
        print(self.B.shape, self.C.shape)
        print(self.jac_x.shape, self.jac_y.shape)

    def AFP(self, F, F0=None, N=10, save_interval=None, save_dir='.', plot=False, epsilon=None):
        epsilon = epsilon or self.epsilon
        phi = np.zeros(F.size)
        F_mean = F.mean()
        if F0 is None:
            F0 = F_mean
        # It is essential for stability of the method that mean(F0) == mean(F)
        F0 *= np.mean(F) / np.mean(F0)
        F0 = F0.flat

        if plot:
            plotter = pplot.PlotPROBLEM(self.X0, self.Y0, F)

        F_interp = interp.RegularGridInterpolator(
            (self.x, self.y),
            F,
            method='linear',
            bounds_error=False,
            fill_value=F_mean
        )

        def get_dphi(phi):
            dX, dY = self.gradient(phi)
            dXx, dYy, dXy = self.hessian(phi)
            H = (1 + dXx) * (1 + dYy) - dXy**2
            L_min = 1 + (dXx + dYy)/2 - np.sqrt((dXx - dYy)**2 + 4*dXy**2)/2;
            gamma = np.maximum(epsilon - L_min, 0)
            j_xx, j_yy, j_xy = spj.second_derivatives(self.x, self.y,
                1 + dYy + gamma,
                1 + dXx + gamma,
                dXy
            )
            J = j_xx + j_yy 
            X = self.X0.flat + dX
            Y = self.Y0.flat + dY
            F_defl = F_interp((X, Y))
            S = self.B @ (F0 / F_defl - H)
            BJC = self.B @ J @ self.C
            dphi = la.spsolve(BJC, S, use_umfpack=True)
            dphi = np.insert(dphi, 0, 0)
            dphi -= dphi.mean()
            return dphi

        for step in np.arange(N):
            phi = phi + get_dphi(phi)
            if save_interval and step % save_interval == 0:
                ph, X, Y, J = self.get_outputs(phi)
                AFP.save_outputs(ph, X, Y, J, save_dir, step)
                if plot:
                    plotter.update(X, Y, f'Step {step}')

        return self.get_outputs(phi)

    def gradient(self, A_flat):
        return self.jac_x @ A_flat, self.jac_y @ A_flat

    def hessian(self, A_flat):
        return self.jac_xx @ A_flat, self.jac_yy @ A_flat, self.jac_xy @ A_flat

    def get_outputs(self, phi_flat):
        dX, dY = self.gradient(phi_flat)
        X = self.X0 + dX.reshape(self.shape_2d)
        Y = self.Y0 + dY.reshape(self.shape_2d)
        dXx, dYy, _ = self.hessian(phi_flat)
        J = (dXx + dYy).reshape(self.shape_2d)
        return phi_flat.reshape(self.shape_2d), X, Y, J

    @staticmethod
    def save_outputs(phi, X, Y, J, save_dir, step):
        # Save the files to the specified save directory.
        np.savetxt(f'{save_dir}/Phi{step}.txt',
                   phi, delimiter=',')
        np.savetxt(f'{save_dir}/X{step}.txt',
                   X, delimiter=',')
        np.savetxt(f'{save_dir}/Y{step}.txt',
                   Y, delimiter=',')
        np.savetxt(f'{save_dir}/J{step}.txt',
                   J, delimiter=',')

