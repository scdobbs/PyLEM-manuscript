slope = 5E-4

m = 0.1
ka = 0.379
h = 1.895

filename = '/home/groups/hilley/scdobbs_hilley_lem/lem_outputs/conc_01_Rf05percent'
from pylem import build_model_dzdt, randomized_grid

dx = 10
nx = 4000
ny = 2000

l = 30.0
L = 10000.0
Rf = 1000.0
time_to_steady_state = 2.0E6
Pe = 3.0
scaling_factor = 1.7081925

import numpy as np
t_eval = np.arange(10.0, time_to_steady_state*1.5, time_to_steady_state*1.5/150.0)

from pylem.utils import Checkpointer
checkpointer = Checkpointer(filename)

f_dzdt = build_model_dzdt((ny, nx), dx, l, L, Rf, time_to_steady_state, Pe, ka, h, m, hook = checkpointer, renoise = 1E-4,renoise_interval = 1000.0, scaling_factor = scaling_factor)
z0 = randomized_grid((ny, nx), slope = slope*dx, noise_level = 0.0)
import pickle as p
z0 = z0 + p.load(open('/home/groups/hilley/scdobbs_hilley_lem/lem_inputs/random_grid.p', 'rb'))

from scipy.integrate import solve_ivp
import numpy as np

results = solve_ivp(f_dzdt, y0=np.reshape(z0, (ny*nx,)), t_span =(0.0, time_to_steady_state*1.5), method='DOP853', t_eval = t_eval)


p.dump(results, open(filename + '_results.p', 'wb'))

quit()
