# -*- coding: utf-8 -*-
"""
Created on Wed May 29 16:08:22 2024

@author: Nate Whitsett
"""
#%%
import os
from matplotlib import pyplot as plt
from matplotlib import font_manager
import numpy as np
import pandas as pd
import glob
from astropy.io import fits
import collections as c
import lightkurve as lk
from ardor.Utils.Utils import find_nearest
# Configure Latin Modern Roman font
font_path = '/ugrad/whitsett.n/fonts/latin-modern-roman/lmroman10-regular.otf'
font_manager.fontManager.addfont(font_path)
plt.rcParams['font.family'] = 'Latin Modern Roman'
plt.rcParams['mathtext.fontset'] = 'cm'  # Computer Modern for math text

# %%

def TESS_data_extract_local(fits_lc_file, PDCSAP_ERR=True, apply_quality_filter=True):
    """
    Extract time and PDCSAP_FLUX from a TESS light curve FITS file.
    
    This is a local version to avoid importing the full ardor.Flares.Flare module.
    """
    lc = lk.read(fits_lc_file, flux_column='pdcsap_flux').remove_nans()
    
    # Print diagnostic information about data filtering
    
    # Apply quality flag filtering (optional)
    if apply_quality_filter:
        quality_filtered_lc = lc[(lc.quality == 0) | (lc.quality == 64)]
        print(len(lc[lc.quality == 1]))
        removed_by_quality = len(lc) - len(quality_filtered_lc)
        print(f"  After quality filtering: {len(quality_filtered_lc)} points (removed {removed_by_quality} with quality != 0)")
        if removed_by_quality > 0:
            unique_quality_flags = np.unique(lc.quality)
            print(f"  Quality flags present: {unique_quality_flags}")
        lc = quality_filtered_lc
    else:
        print("  Quality filtering disabled - keeping all data regardless of quality flags")
    
    # Normalize flux to median
    flux = lc.flux / np.median(lc.flux)
    error = lc.flux_err / np.median(lc.flux)
    
    # Convert time to numeric values (BTJD, already BJD - 2457000)
    time_values = lc.time.value
    flux_values = flux.value
    error_values = error.value
    
    if PDCSAP_ERR:
        LightCurve = c.namedtuple('LightCurve', ['time', 'flux', 'error'])
        lc_out = LightCurve(time_values, flux_values, error_values)
    else:
        LightCurve = c.namedtuple('LightCurve', ['time', 'flux'])
        lc_out = LightCurve(time_values, flux_values)
    
    return lc_out


def plot_flares_from_catalog(csv_file, host_id, fits_dir, flare_lengths=None, 
                              window_points=50, output_file=None, plot_title=None,
                              max_flares_per_page=30, apply_quality_filter=True, 
                              load_all_data=False):
    """
    Plot flares from a CSV catalog centered on flare events.
    
    Creates multi-page plots if there are more than max_flares_per_page flares.
    
    Parameters
    ----------
    csv_file : str
        Path to the CSV catalog file containing flare information.
        Expected to have columns for host name (first column) and flare epoch (4th column).
    host_id : str
        Name of the star/host as it appears in the first column of the catalog.
    fits_dir : str
        Directory containing the FITS light curve files for the specified host.
    flare_lengths : list or None, optional
        List of flare lengths (in number of data points) for each flare.
        If provided, highlights subsequent data points after the flare peak.
        Must be the same length as the number of flares for the host. Default is None.
    window_points : int, optional
        Number of data points to show before and after the flare center (±window_points).
        Default is 50.
    output_file : str or None, optional
        If provided, saves the plot to this file path. For multi-page plots, 
        appends "_page1", "_page2", etc. before the extension. Default is None (display only).
    plot_title : str or None, optional
        Custom title for the plot. If None, uses "{host_id} - Flare Events". Default is None.
    max_flares_per_page : int, optional
        Maximum number of flares to show on each page. Default is 30.
    apply_quality_filter : bool, optional
        Whether to apply TESS quality flag filtering (quality == 0). If False, uses all data
        regardless of quality flags. Default is True.
    load_all_data : bool, optional
        If True, loads both fast and regular data regardless of Fast_Bool values. 
        Useful for seeing full context. Default is False.
    
    Returns
    -------
    figures : list of matplotlib.figure.Figure
        List of figure objects (one per page).
    all_axes : list of numpy.ndarray
        List of axes arrays (one array per page).
    
    Notes
    -----
    - Time values are shown as BJD - 2457000 to save space.
    - Flare peaks are automatically detected and highlighted in red triangles.
    - If flare_lengths is provided, subsequent points after the peak are highlighted in green squares.
    - Flare numbers are shown as annotations in the upper right of each panel.
    - Y-axis ticks are moved to the right side for the rightmost column to minimize white space.
    - The plot uses a standard page size (8.5 x 11 inches) with 1-inch margins.
    """
    
    
    # Tier3 format: Has header + units rows (skiprows=2), 19 columns
    catalog = pd.read_csv(csv_file)

    # Sort by the flare number column ('#') to ensure flares are plotted in catalog order
    if '#' in catalog.columns:
        catalog = catalog.sort_values('#').reset_index(drop=True)

    # Convert to numeric
    catalog['Epoch_TESS'] = pd.to_numeric(catalog['Epoch_TESS'], errors='coerce')
    catalog['FWHM'] = pd.to_numeric(catalog['FWHM'], errors='coerce')
    catalog['Amplitude'] = pd.to_numeric(catalog['Amplitude'], errors='coerce')
    
    has_fwhm = True
    has_amplitude = True
    # Filter for the specified host
    host_flares = catalog[catalog['Host_ID'] == host_id].copy()
    
    if len(host_flares) == 0:
        raise ValueError(f"No flares found for host '{host_id}' in the catalog.")
    
    # Get flare epochs (use the BJD column - 4th column) and subtract 2457000
    flare_epochs = host_flares['Epoch_TESS'].values
    
    # Get amplitudes and FWHM if available
    if has_amplitude:
        flare_amplitudes = host_flares['Amplitude'].values
    else:
        flare_amplitudes = None
    
    if has_fwhm:
        fwhm_values = host_flares['FWHM'].values
    else:
        fwhm_values = None
    
    n_flares = len(flare_epochs)
    
    # Validate flare_lengths if provided
    if flare_lengths is not None:
        if len(flare_lengths) != n_flares:
            raise ValueError(f"flare_lengths must have {n_flares} elements (one per flare), "
                           f"but got {len(flare_lengths)}.")
    
    # Find all FITS files for this host
    all_fits_files = sorted(glob.glob(os.path.join(fits_dir, '*.fits')))
    
    if len(all_fits_files) == 0:
        raise ValueError(f"No FITS files found in directory: {fits_dir}")
    
    # Separate fast and regular FITS files
    fast_fits = [f for f in all_fits_files if f.endswith('a_fast-lc.fits')]
    regular_fits = [f for f in all_fits_files if f.endswith('_lc.fits') and not f.endswith('a_fast-lc.fits')]
    
    print(f"\nFile discovery summary:")
    print(f"  Total FITS files found: {len(all_fits_files)}")
    print(f"  Fast files found: {len(fast_fits)}")
    for f in fast_fits:
        print(f"    {os.path.basename(f)}")
    print(f"  Regular files found: {len(regular_fits)}")
    for f in regular_fits:
        print(f"    {os.path.basename(f)}")
    
    # Create a list of files to use based on Fast_Bool from catalog
    fits_files_to_use = []
    
    # Load ALL available data regardless of Fast_Bool values
    fits_files_to_use.extend(fast_fits)
    fits_files_to_use.extend(regular_fits)
    
    # Sort the final list
    fits_files_to_use = sorted(fits_files_to_use)
    print(f"Final list of FITS files to use ({len(fits_files_to_use)} files):")
    for f in fits_files_to_use:
        print(f"  {os.path.basename(f)}")
    
    # Read all light curve data and track file types
    all_time = []
    all_flux = []
    all_error = []
    all_is_fast = []  # Track which data points come from fast files
    
    total_fast_points = 0
    total_regular_points = 0
    
    for fits_file in fits_files_to_use:
        try:
            lc = TESS_data_extract_local(fits_file, PDCSAP_ERR=True, apply_quality_filter=apply_quality_filter)
            is_fast_file = fits_file.endswith('a_fast-lc.fits')
            
            all_time.extend(lc.time)
            all_flux.extend(lc.flux)
            all_error.extend(lc.error)
            all_is_fast.extend([is_fast_file] * len(lc.time))  # Track file type for each point
            
            if is_fast_file:
                total_fast_points += len(lc.time)
            else:
                total_regular_points += len(lc.time)
        except Exception as e:
            print(f"Warning: Could not read {fits_file}: {e}")
            continue
    
    if len(all_time) == 0:
        raise ValueError("No valid data could be read from FITS files.")
    
    # Convert to numpy arrays and sort by time
    all_time = np.array(all_time)
    all_flux = np.array(all_flux)
    all_error = np.array(all_error)
    all_is_fast = np.array(all_is_fast)
    
    sort_idx = np.argsort(all_time)
    all_time = all_time[sort_idx]
    all_flux = all_flux[sort_idx]
    all_error = all_error[sort_idx]
    all_is_fast = all_is_fast[sort_idx]  # Also sort the file type tracking array
    
    print(f"\nCombined dataset summary:")
    print(f"  Total data points: {len(all_time)}")
    print(f"  Fast data points: {total_fast_points}")
    print(f"  Regular data points: {total_regular_points}")
    print(f"  Time range: {all_time.min():.2f} to {all_time.max():.2f}")
    
    # Calculate cadence (average time spacing between data points) in days
    time_diffs = np.diff(all_time)
    cadence_days = np.median(time_diffs)  # Median spacing between points
    
    # Convert FWHM (in days) to number of data points for highlighting
    # flare_lengths will be used for decay highlighting
    if flare_lengths is None and fwhm_values is not None:
        # Auto-calculate from FWHM: convert days to number of data points
        flare_lengths = np.round(fwhm_values / cadence_days).astype(int)
    
    # Determine number of pages needed
    n_pages = int(np.ceil(n_flares / max_flares_per_page))
    
    # Set plot title
    if plot_title is None:
        plot_title = f'{host_id} - Flare Events'
    
    # Store all figures and axes
    figures = []
    all_axes = []
    
    # Create plots page by page
    for page_num in range(n_pages):
        # Determine flares for this page
        start_flare = page_num * max_flares_per_page
        end_flare = min((page_num + 1) * max_flares_per_page, n_flares)
        page_flares = end_flare - start_flare
        page_epochs = flare_epochs[start_flare:end_flare]
        page_amplitudes = flare_amplitudes[start_flare:end_flare] if flare_amplitudes is not None else None
        page_fwhm_values = fwhm_values[start_flare:end_flare] if fwhm_values is not None else None
        page_flare_lengths = flare_lengths[start_flare:end_flare] if flare_lengths is not None else None
        
        # Calculate grid layout for this page
        n_cols = min(3, page_flares)  # Maximum 3 columns
        n_rows = int(np.ceil(page_flares / n_cols))
        
        # Create figure with dynamic size based on number of rows
        # Constrained to fit on 7x11 inch page with 0.5 inch margins
        fig_width = 6.0  # 7 - 2*0.5 = 6 inches available width
        # Calculate height: base height per row + space for title and axis labels
        height_per_row = 2.0  # Height per row of subplots (reduced for more compact plots)
        base_height = 1.0     # Space for title and axis labels (reduced)
        fig_height = max(2.5, min(10.0, n_rows * height_per_row + base_height))  # Max 10 inches (11 - 2*0.5)
        
        fig, axes = plt.subplots(n_rows, n_cols, figsize=(fig_width, fig_height))
        
        # Adjust margins dynamically based on figure height
        # For smaller figures, use larger relative margins for labels
        bottom_margin = max(0.06, 0.8 / fig_height)  # Ensure space for x-label
        top_margin = max(0.06, 0.8 / fig_height)     # Ensure space for title
        plt.subplots_adjust(left=0.08, right=0.97, 
                           bottom=bottom_margin, top=1.0 - top_margin, 
                           hspace=0, wspace=0)
        
        # Flatten axes array for easier indexing
        if page_flares == 1:
            axes = np.array([axes])
        else:
            axes = axes.flatten()
        
        # Plot each flare on this page
        for i, (epoch, ax) in enumerate(zip(page_epochs, axes[:page_flares])):
            global_flare_num = start_flare + i + 1  # Global flare number
        
            filtered_time = all_time
            filtered_flux = all_flux
            filtered_error = all_error
            print(f"Flare {global_flare_num}: Using all data ({len(filtered_time)} points)")
            
            # Find data points closest to the flare epoch in filtered data
            center_idx = find_nearest(filtered_time, epoch)
            # Define window around flare
            start_idx = max(0, center_idx - window_points)
            end_idx = min(len(filtered_time), center_idx + window_points + 1)
            
            # Extract window data
            time_window = filtered_time[start_idx:end_idx]
            flux_window = filtered_flux[start_idx:end_idx]
            error_window = filtered_error[start_idx:end_idx]
            
            # Plot light curve
            ax.errorbar(time_window, flux_window, yerr=error_window, 
                       fmt='-o', markersize=3, color='black', alpha=0.7, 
                       elinewidth=0.5, capsize=0, label='Data')
            
            # Highlight the center point (epoch) with an orange diamond
            center_time = epoch
            # Find the closest data point to the epoch in filtered data
            center_flux_idx = np.argmin(np.abs(filtered_time - epoch))
            
            # Find the maximum flux point within a 5-point range around the epoch
            # Create a range of ±2 points around the center (5 points total)
            range_start = max(0, center_flux_idx - 2)
            range_end = min(len(filtered_flux), center_flux_idx + 3)  # +3 because range() is exclusive
            search_indices = list(range(range_start, range_end))
            
            # Find the index with maximum flux in this 5-point range
            if len(search_indices) > 0:
                peak_idx = search_indices[np.argmax([filtered_flux[i] for i in search_indices])]
            else:
                peak_idx = center_flux_idx  # Fallback to center if range is invalid
            
            center_flux = filtered_flux[peak_idx]
            
            ax.plot(filtered_time[peak_idx], center_flux, marker='D', markersize=6, 
                    color='orange', label='Flare Peak', zorder=5)
            
            # Use amplitude from catalog if available, otherwise calculate from flux
            if page_amplitudes is not None:
                flare_amplitude = page_amplitudes[i]
            else:
                flare_amplitude = center_flux - 1.0
            
            # Highlight flare decay if flare_lengths is provided
            if page_flare_lengths is not None and page_flare_lengths[i] > 0:
                decay_length = 2
                # Find the epoch point in the filtered light curve
                epoch_idx = np.argmin(np.abs(filtered_time - epoch))
                decay_end_idx = min(epoch_idx + decay_length, len(filtered_time))
                
                if decay_end_idx > epoch_idx + 1:
                    decay_time = filtered_time[epoch_idx+1:decay_end_idx]
                    decay_flux = filtered_flux[epoch_idx+1:decay_end_idx]
                    ax.plot(decay_time, decay_flux, marker='o', markersize=2.5, 
                           color='#00CED1', linestyle='', alpha=0.9, 
                           label='Flare Decay', zorder=4)
            
            # Add flare number as annotation (upper right)
            ax.text(0.95, 0.95, f'{global_flare_num}', 
                   transform=ax.transAxes, fontsize=8, 
                   verticalalignment='top', horizontalalignment='right',
                   bbox=dict(boxstyle='round', facecolor='white', alpha=0.7, 
                            edgecolor='gray', linewidth=0.5))
            
            # Add time and amplitude annotations (upper left)
            time_text = f'$t_\\text{{mid}} = {center_time:.2f}$'
            amp_text = f'$A_{{\\text{{peak}}}} = {flare_amplitude:.3f}$'
            
            # Build annotation text based on available data
            annotation_parts = [time_text, amp_text]
            if page_fwhm_values is not None:
                fwhm_value = page_fwhm_values[i]
                fwhm_text = f'$t_\\text{{FWHM}} = {1440*fwhm_value:.3f} \\, \\mathrm{{[min]}}$'
                annotation_parts.append(fwhm_text)
            
            ax.text(0.02, 0.95, '\n'.join(annotation_parts), 
                   transform=ax.transAxes, fontsize=6.5, 
                   verticalalignment='top', horizontalalignment='left',
                   bbox=dict(boxstyle='round', facecolor='white', alpha=0.7, 
                            edgecolor='gray', linewidth=0.5),
                   family='monospace')
            
            # Formatting - add ticks pointing inward
            ax.tick_params(direction='in', length=3, width=0.5, colors='black')
            # ax.set_xticklabels([])  # No tick labels
            # ax.set_yticklabels([])  # No tick labels
            ax.grid(True, alpha=0.3, linestyle='--', linewidth=0.5)
            
            # Set y-axis limits based on data in the window
            if len(flux_window) > 0:
                max_flux = np.max(flux_window)
                min_flux = np.min(flux_window)
                # Set upper limit to 120% of maximum, lower limit with some padding
                y_upper = max_flux * 1.001
                y_lower = min_flux - min_flux * 0.001  # 10% padding below minimum
                ax.set_ylim(y_lower, y_upper)
            
            # Move y-axis to right side for rightmost column
            row = i // n_cols
            col = i % n_cols
            if col == n_cols - 1:  # Rightmost column
                ax.yaxis.tick_right()
                ax.yaxis.set_label_position("right")
            
            # Only show legend on first subplot of first page
            # if page_num == 0 and i == 0:
            #     ax.legend(fontsize=6, loc='upper left', framealpha=0.9)
        
        # Hide unused subplots
        for i in range(page_flares, len(axes)):
            axes[i].set_visible(False)
        
        # Add unified axis labels
        fig.text(0.5, 0.015, 'Time (BJD - 2457000)', ha='center', fontsize=11, weight='bold')
        fig.text(0.02, 0.5, 'Rel. Flux', va='center', rotation='vertical', 
                fontsize=11, weight='bold')
        
        # Add title with page number if multiple pages
        if n_pages > 1:
            title = f'{plot_title} (Page {page_num + 1}/{n_pages})'
        else:
            title = plot_title
        fig.suptitle(title, fontsize=13, weight='bold', y=0.98)
        
        # Save or show
        if output_file:
            if n_pages > 1:
                # Add page number to filename
                base, ext = os.path.splitext(output_file)
                page_file = f"{base}_page{page_num + 1}{ext}"
            else:
                page_file = output_file
            
            plt.savefig(page_file, dpi=300, bbox_inches='tight')
            print(f"Plot saved to: {page_file}")
        
        figures.append(fig)
        all_axes.append(axes)
    
    return figures, all_axes


def list_catalog_hosts(csv_file, min_flares=1):
    """
    List all hosts in a flare catalog with their flare counts.
    
    Parameters
    ----------
    csv_file : str
        Path to the CSV catalog file.
    min_flares : int, optional
        Only show hosts with at least this many flares. Default is 1.
    
    Returns
    -------
    host_summary : pandas.DataFrame
        DataFrame with columns: 'Host_ID', 'Num_Flares'
        Sorted by number of flares (descending).
    
    Examples
    --------
    >>> summary = list_catalog_hosts('All_Flare_Parameters.csv', min_flares=5)
    >>> print(summary)
    """
    catalog = pd.read_csv(csv_file, header=None)
    host_counts = catalog[0].value_counts().reset_index()
    host_counts.columns = ['Host_ID', 'Num_Flares']
    
    # Filter by minimum flares
    host_counts = host_counts[host_counts['Num_Flares'] >= min_flares]
    
    return host_counts


def get_host_flare_epochs(csv_file, host_id):
    """
    Get all flare epochs for a specific host from the catalog.
    
    Parameters
    ----------
    csv_file : str
        Path to the CSV catalog file.
    host_id : str
        Name of the host/star.
    
    Returns
    -------
    epochs : numpy.ndarray
        Array of flare epochs in BJD for the specified host.
    
    Raises
    ------
    ValueError
        If no flares found for the specified host.
    
    Examples
    --------
    >>> epochs = get_host_flare_epochs('All_Flare_Parameters.csv', '51Peg')
    >>> print(f"Found {len(epochs)} flares")
    """
    catalog = pd.read_csv(csv_file, header=None)
    host_flares = catalog[catalog[0] == host_id]
    
    if len(host_flares) == 0:
        raise ValueError(f"No flares found for host '{host_id}' in the catalog.")
    
    # Column 3 (index 3) is Epoch_BJD
    epochs = host_flares[3].values
    
    return epochs

# Test code (commented out)
figures, axes = plot_flares_from_catalog(
csv_file='/ugrad/whitsett.n/Induced_Flares/Flare_Data/Tier_3/Grand_Lists/All_Exoplanet_MCMC_Flares.csv',
host_id='Gl49',
fits_dir='/data2/whitsett.n/TESS/Hosts/Gl49',
output_file='/ugrad/whitsett.n/Gl49.png',
plot_title='Gliese 49 - Flare Observations',
window_points=50, max_flares_per_page=30, apply_quality_filter=False
)