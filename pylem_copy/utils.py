import numpy as np
import pickle as p
from pylem import unfreeze_from_checkpoint_file

class Checkpointer(object):

    def __init__(self, filename, checkpoint_every = 10, output_every = 0.1, model_data = None):
        self.filename = filename
        self.checkpoint_every = checkpoint_every
        self.output_every = output_every
        self.flux = []
        self.model_times = []
        self.counter = 0
        self.last_flux_output = 0
        self.model_data = model_data

    def __calculate_flux(self, dzdt, U):
        return 1 - np.mean(dzdt/U)

    def __checkpoint_model(self, t, y):
        p.dump((t, y, self), open(self.filename + "_checkpoint.p", 'wb'))
        print("Time is: ", t, flush = True)

    def __output_model(self, t, y):
        p.dump((t, y), open("{filename}_{fraction:.1f}.p".format(filename = self.filename, fraction = self.last_flux_output+self.output_every), 'wb'))
        self.last_flux_output += self.output_every

    def register_new_step(self, t, y, dzdt, dx, U):
        flux_fraction = self.__calculate_flux(dzdt, U)
        self.model_times += [t]
        self.flux += [flux_fraction]
        if flux_fraction >= (self.last_flux_output + self.output_every):
            self.__output_model(t, y)
        self.counter += 1
        if self.counter > self.checkpoint_every:
            self.counter = 0
            self.__checkpoint_model(t, y)

def restart_model(filename, method, max_step = np.inf):

    f_dzdt, checkpointer, t0, y0, time_to_steady_state, (ny, nx) = unfreeze_from_checkpoint_file(filename)

    t_eval = np.array(list(range(10)))*time_to_steady_state*1.5/10.0 + 10.0
    i = np.where(t_eval > t0)
    t_eval = t_eval[i]

    from scipy.integrate import solve_ivp

    results = solve_ivp(f_dzdt, y0=np.reshape(y0, (ny*nx,)), t_span =(t0, time_to_steady_state*1.5), method=method, t_eval = t_eval, max_step = max_step)

    import pickle as p
    p.dump(results, open(filename + '_results.p', 'wb'))
