from pylem.pyas import area_with_fill_dinf as area
import numpy as np
import random
call_count = 0

def create_optimization_function(dem, ks, concavity, dx=None, interpolation='quintic'):
    if dx is not None:
        dem = dem.resample(dx, interpolation=interpolation)
        dem._griddata = dem._griddata.copy()

    m, n = dem._griddata.shape
    dx = dem._georef_info.dx
    create_optimization_function.nt = 0

    def misfit(values):
        global call_count
        call_count += 1

        # Calculate the original misfit
        grid = np.reshape(values, (m, n))
        grid = np.pad(grid, ((1, 1), (0, 0)), 'constant', constant_values=0)
        a, s, filled_dem = area(grid, dx=dx)  # Use the new function with filled DEM
        ks_obs_full = np.power(a, concavity) * s
        ks_obs = np.reshape(ks_obs_full[1:-1,:], (m * n,))
        sse = np.sum(np.power(ks_obs - ks, 2))

        # Apply penalty only where the grid is lower than the filled DEM
        filled_dem_penalty = np.sum(np.abs(np.minimum(grid - filled_dem, 0))**4)
        penalty_weight = 1e6  # Adjust this weight as needed

        # Adjusted misfit calculation
        total_misfit = np.sqrt(sse / (m*n)) + penalty_weight * filled_dem_penalty

        create_optimization_function.nt += 1
        if create_optimization_function.nt == 2000:
            print(
                f"Iteration {call_count}: Total Misfit = {total_misfit}, Mean ks_obs / ks = {np.mean(ks_obs) / ks}, "
                f"SSE = {np.sqrt(sse / (m * n)) / ks}", flush=True)

        return total_misfit

    return np.reshape(dem._griddata, (m*n, )), misfit, (m, n), dem
