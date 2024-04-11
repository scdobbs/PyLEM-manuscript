import numpy as np
import matplotlib.pylab as plt
from heapq import heappush, heappop
from .pyas import area_dinf as area
import pickle as p

def calc_K_U_D(l, L, Rf, time_to_steady_state, Pe, ka, h, m):

    K = np.power(time_to_steady_state,-1.0)*np.power(ka,-m)*np.power(1-h*m,-1)*(np.power(L,1-h*m)-np.power(l,1-h*m))  # Calculate K
    U = Rf * K *np.power(ka,m)*(1-h*m)*np.power(np.power(L,1-h*m)-np.power(l,1-h*m),-1)  # Calculate U
    D = K * np.power(l,h*m+1) * np.power(ka,m) / Pe  # Calculate D
    Final_relief_check = (U/K)*np.power(ka,-m)*np.power(1-h*m,-1)*(np.power(L,1-h*m)-np.power(l,1-h*m))

    return K, U, D

import numpy as np

def build_model_dzdt(size, dx, l, L, Rf, time_to_steady_state, Pe, ka, h, m, hook = None, renoise = None, return_dt = False):

    K, U, D = calc_K_U_D(l, L, Rf, time_to_steady_state, Pe, ka, h, m)

    (ny, nx) = size
    build_model_dzdt.counter = 0

    if hook is not None:
        hook.model_data = {"dx": dx,
                           "l": l,
                           "L": L,
                           "Rf": Rf,
                           "time_to_steady_state": time_to_steady_state,
                           "Pe": Pe,
                           "ka": ka,
                           "h": h,
                           "m": m,
                           "size": size,
                           "K": K,
                           "U": U,
                           "D": D,
                           "renoise": renoise}

    def dzdt(t, y):
        from numpy.random import rand
        z = np.reshape(y, (ny, nx))
        if renoise is not None:
            z += np.reshape(rand(nx*ny), (ny, nx))*renoise
        Qx = -D*np.diff(np.hstack((np.reshape(z[:,-1],(ny, 1)), z, np.reshape(z[:,0], (ny, 1)))), axis = 1)/dx
        Qy = -D*np.diff(z, axis = 0)/dx
        Qy = np.vstack((Qy[0,:]-np.ones((1,nx))*U*dx, Qy, np.ones((1,nx))*U*dx+Qy[-1,:]))
        dzdt_diffusion = -np.diff(Qx, axis = 1)/dx - np.diff(Qy, axis = 0)/dx
        a, s = area(z, dx)
        build_model_dzdt.counter += 1
        dzdt_erosion = -K*np.power(a,m)*s
        dzdt = U + dzdt_diffusion + dzdt_erosion
        dzdt[0,:] = 0.0
        dzdt[-1,:] = 0.0
        if hook is not None:
            hook.register_new_step(t, y, dzdt, dx, U)
        return np.reshape(dzdt, (ny*nx,))

    return dzdt

def randomized_grid(shape, noise_level = 1.0, slope = 0.0):
    from numpy.random import rand
    from numpy.matlib import repmat
    (ny, nx) = shape
    z0 = np.reshape(rand(nx*ny), (ny, nx))*noise_level - slope*np.abs(repmat(np.reshape(np.array(list(range(ny))) - 0.5*ny, (ny, 1)), 1, nx)) + ny*slope / 2
    z0[0,:] = 0.0
    z0[-1,:] = 0.0
    return z0

def unfreeze_from_checkpoint_file(filename):

    (t, y, checkpointer) = p.load(open(filename + '_checkpoint.p', 'rb'))
    md = checkpointer.model_data
    (dx, K, U, D, m, ny, nx, tss) = (md['dx'], md['K'], md['U'], md['D'], md['m'], md['size'][0], md['size'][1], md['time_to_steady_state'])

    unfreeze_from_checkpoint_file.counter = 1

    def dzdt(t, y):

        z = np.reshape(y, (ny, nx))
        Qx = -D*np.diff(np.hstack((np.reshape(z[:,-1],(ny, 1)), z, np.reshape(z[:,0], (ny, 1)))), axis = 1)/dx
        Qy = -D*np.diff(z, axis = 0)/dx
        Qy = np.vstack((Qy[0,:]-np.ones((1,nx))*U*dx, Qy, np.ones((1,nx))*U*dx+Qy[-1,:]))
        dzdt_diffusion = -np.diff(Qx, axis = 1)/dx - np.diff(Qy, axis = 0)/dx
        a, s = area(z, dx)
        unfreeze_from_checkpoint_file.counter += 1
        dzdt_erosion = -K*np.power(a,m)*s
        dzdt = U + dzdt_diffusion + dzdt_erosion
        dzdt[0,:] = 0.0
        dzdt[-1,:] = 0.0
        if checkpointer is not None:
            checkpointer.register_new_step(t, y, dzdt, dx, U)

        return np.reshape(dzdt, (ny*nx,))

    return dzdt, checkpointer, t, y, tss, (ny,nx)
