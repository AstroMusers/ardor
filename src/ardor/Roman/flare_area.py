from ardor.Flares.aflare import flare_total_area
from ardor.Flares.allesfitter_priors import flare_energy
from ardor.SPI_Forward_Models.SPI_Simulation import amp_log_normal, FWHM_uniform
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from scipy.ndimage import gaussian_filter1d
from scipy.interpolate import CubicSpline, interp1d
from scipy.integrate import quad
from ardor.Utils.planck_law import planck_integrator

n_radius_bins = 30


def run_flare_area_analysis():
    radius_file = '/ugrad/whitsett.n/Roman/teff_binned_smoothed_radius.csv'
    if not os.path.exists(radius_file):
        print(f"Data file not found: {radius_file}")
        return
    radii = pd.read_csv(radius_file)
    test = radii['teff_3500_3600_K'].to_numpy()
    rad = radii['radius_bin_center'].to_numpy()
    function = CubicSpline(rad, test, extrapolate=False)
    # Integrate over actual data bounds, not extrapolated range
    norm = function.integrate(np.min(rad), np.max(rad))
    # Normalize by scaling coefficients
    function.c /= norm
    # Use actual radius bounds from data to avoid NaNs from extrapolation
    x_grid = np.linspace(np.min(rad), np.max(rad), 1000)
    pdf_vals = function(x_grid)
    color_factor = planck_integrator(600e-6, 1000e-6, 3550)/planck_integrator(600e-6, 1000e-6, 9000)
    cdf = np.cumsum(pdf_vals)
    cdf /= cdf[-1]  # Normalize CDF to range [0, 1]
    inv_cdf = interp1d(cdf, x_grid, bounds_error=False, fill_value=(x_grid[0], x_grid[-1]))
    u = np.random.rand(100)  # 10k samples
    samples = inv_cdf(u)

    flare_file = '/ugrad/whitsett.n/Roman/Flare_Catalogs/All_Exoplanet_MCMC_Flares.csv'
    if os.path.exists(flare_file):
        flare_data = pd.read_csv(flare_file)
        print(flare_data.loc[(flare_data['FWHM'] > 0.02), ['Host_ID', '#']])


if __name__ == '__main__':
    run_flare_area_analysis()



# temps = np.arange(3500, 6000 + 500, step=500)
# radius_by_temp = {temp: [] for temp in temps}

# # Gather stellar radii for each Teff bin across all chunks.
# for chunk in pd.read_csv(tic_path, chunksize=chunk_size):
#     quality_mask = (chunk.iloc[:, 66] < 5) & (chunk.iloc[:, 66] > 4)
#     teff = chunk.iloc[:, 64]
#     r_stellar = chunk.iloc[:, 70]

#     for temp in temps:
#         temp_mask = quality_mask & (teff >= temp) & (teff < temp + 100)
#         selected = r_stellar.loc[temp_mask].to_numpy()
#         selected = selected[np.isfinite(selected)]
#         if selected.size > 0:
#             radius_by_temp[temp].append(selected)
#         print(temp)

# all_radii = [arr for arrays in radius_by_temp.values() for arr in arrays]
# if not all_radii:
#     raise RuntimeError('No stellar radii found for the requested Teff bins.')

# all_radii = np.concatenate(all_radii)

# # Build one shared radius grid and reuse it for every Teff bin.
# radius_grid_edges = np.linspace(np.nanmin(all_radii), np.nanmax(all_radii), n_radius_bins + 1)
# radius_grid_centers = 0.5 * (radius_grid_edges[:-1] + radius_grid_edges[1:])

# output_df = pd.DataFrame({'radius_bin_center': radius_grid_centers})

# for temp in temps:
#     col_name = f'teff_{temp}_{temp + 100}_K'
#     if radius_by_temp[temp]:
#         radii = np.concatenate(radius_by_temp[temp])
#         counts, _ = np.histogram(radii, bins=radius_grid_edges, density=True)
#         # Slightly stronger smoothing than before to reduce noise in sparse bins.
#         smoothed = gaussian_filter1d(counts, sigma=2.0)
#     else:
#         smoothed = np.full(radius_grid_centers.shape, np.nan)

#     output_df[col_name] = smoothed

# output_df.to_csv(output_path, index=False)
# print(f'Saved smoothed Teff-bin data to {output_path}')
