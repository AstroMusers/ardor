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
import matplotlib.patches as mpatches
import textwrap
# Configure Latin Modern Roman font
font_path = '/ugrad/whitsett.n/fonts/latin-modern-roman/lmroman10-regular.otf'
font_manager.fontManager.addfont(font_path)
plt.rcParams['font.family'] = 'Latin Modern Roman'
plt.rcParams['mathtext.fontset'] = 'cm'  # Computer Modern for math text


def place_column_text(ax, text, xy, wrap_n, shift, bbox=False, **kwargs):
    """ Creates a text annotation with the text in columns.
    The text columns are provided by a list of strings.
    A surrounding box can be added via bbox=True parameter.
    If so, FancyBboxPatch kwargs can be specified.
    
    The width of the column can be specified by wrap_n,
    the shift parameter determines how far apart the columns are.
    The axes are specified by the ax parameter.

    Requires:
    import textwrap
    import matplotlib.patches as mpatches
    """
    # place the individual text boxes, with a bbox to extract details from later
    x,y = xy
    n = 0
    text_boxes = []
    for i in text:
        text = textwrap.fill(i, wrap_n)
        box = ax.text(x = x + n, y = y, s=text, va='top', ha='left',
                         bbox=dict(alpha=0, boxstyle='square,pad=0'))
        text_boxes.append(box)
        n += shift
    
    if bbox == True: # draw surrounding box
        # extract box data
        plt.draw() # so we can extract real bbox data
        # first let's calulate the height of the largest bbox
        heights=[]
        for box in text_boxes:
            heights.append(box.get_bbox_patch().get_extents().transformed(ax.transData.inverted()).bounds[3])
        max_height=max(heights)
        # then calculate the furthest x value of the last bbox
        end_x = text_boxes[-1].get_window_extent().transformed(ax.transData.inverted()).xmax
        # draw final
        width = end_x - x
        fancypatch_y = y - max_height
        rect = mpatches.FancyBboxPatch(xy=(x,fancypatch_y), width=width, height=max_height, **kwargs)
        ax.add_patch(rect)

def TESS_data_extract_local(fits_lc_file, PDCSAP_ERR=True, apply_quality_filter=True):
    """
    Extract time and PDCSAP_FLUX from a TESS light curve FITS file.
    
    This is a local version to avoid importing the full ardor.Flares.Flare module.
    """
    lc = lk.read(fits_lc_file, flux_column='pdcsap_flux').remove_nans()
    
    # Print diagnostic information about data filtering
    
    # Apply quality flag filtering (optional)
    if apply_quality_filter:
        quality_filtered_lc = lc[(lc.quality == 0)]
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
                              fast = False, plot_both_cadences=False, cadence_list=None):
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
    window_points : int or list, optional
        Number of data points to show before and after the flare center (±window_points).
        Can be either:
        - int: Same window size for all flares (e.g., 50)
        - list: Individual window size for each flare (must match number of flares)
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
    fast : bool, optional
        If True and plot_both_cadences is False, uses only fast cadence data.
        If False and plot_both_cadences is False, uses only regular cadence data. Default is False.
    plot_both_cadences : bool, optional
        If True, plots both regular and fast cadence data on the same plots with different colors.
        Regular data shown in black, fast data shown in blue. Default is False.
    cadence_list : list of str or None, optional
        List specifying cadence type for each flare. Must be same length as number of flares.
        Each element should be 'r' or 'regular' for regular cadence, 'f' or 'fast' for fast cadence.
        Example: ['r', 'f', 'r', 'r', 'f', 'f'] for 6 flares.
        If provided, overrides the 'fast' and 'plot_both_cadences' parameters. Default is None.
    
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
    catalog.sort_values(by='Epoch_BJD', inplace=True)
    catalog['Epoch_BJD'] = pd.to_numeric(catalog['Epoch_BJD'], errors='coerce')
    catalog['FWHM'] = pd.to_numeric(catalog['FWHM'], errors='coerce')
    catalog['Amplitude'] = pd.to_numeric(catalog['Amplitude'], errors='coerce')
    catalog['dlogZ'] = pd.to_numeric(catalog['dlogZ'], errors='coerce')
    
    has_fwhm = True
    has_amplitude = True
    has_dlogz = True
    # Filter for the specified host
    host_flares = catalog[catalog['Host_ID'] == host_id].copy()
    
    if len(host_flares) == 0:
        raise ValueError(f"No flares found for host '{host_id}' in the catalog.")
    
    # Get flare epochs (use the BJD column - 4th column) and subtract 2457000
    flare_epochs = host_flares['Epoch_BJD'].values - 2457000
    
    # Get amplitudes and FWHM if available
    if has_amplitude:
        flare_amplitudes = host_flares['Amplitude'].values
    else:
        flare_amplitudes = None
    
    if has_fwhm:
        fwhm_values = host_flares['FWHM'].values
    else:
        fwhm_values = None

    if has_dlogz:
        dlogz_values = host_flares['dlogZ'].values
    else:
        dlogz_values = None
    
    n_flares = len(flare_epochs)
    
    # Validate and normalize window_points
    if isinstance(window_points, (list, np.ndarray)):
        if len(window_points) != n_flares:
            raise ValueError(f"window_points list must have {n_flares} elements (one per flare), "
                           f"but got {len(window_points)}.")
        window_points_list = list(window_points)
    else:
        # Single value - use for all flares
        window_points_list = [window_points] * n_flares
    
    # Validate flare_lengths if provided
    if flare_lengths is not None:
        if len(flare_lengths) != n_flares:
            raise ValueError(f"flare_lengths must have {n_flares} elements (one per flare), "
                           f"but got {len(flare_lengths)}.")
    
    # Validate and normalize cadence_list if provided
    if cadence_list is not None:
        if len(cadence_list) != n_flares:
            raise ValueError(f"cadence_list must have {n_flares} elements (one per flare), "
                           f"but got {len(cadence_list)}.")
        # Normalize cadence_list to 'r' or 'f'
        normalized_cadence = []
        for i, cadence in enumerate(cadence_list):
            if cadence.lower() in ['r', 'regular']:
                normalized_cadence.append('r')
            elif cadence.lower() in ['f', 'fast']:
                normalized_cadence.append('f')
            else:
                raise ValueError(f"cadence_list[{i}] = '{cadence}' is invalid. "
                               f"Must be 'r', 'regular', 'f', or 'fast'.")
        cadence_list = normalized_cadence
    
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
    
    # Create a list of files to use based on parameters
    fits_files_to_use = []
    
    if cadence_list is not None:
        # Load both types when cadence_list is specified
        print(f"\n📊 Using per-flare cadence specification")
        print(f"   Loading both regular and fast cadence data")
        print(f"   Regular files: {len(regular_fits)}")
        print(f"   Fast files: {len(fast_fits)}")
        fits_files_to_use.extend(regular_fits)
        fits_files_to_use.extend(fast_fits)
    elif plot_both_cadences:
        print(f"\n📊 Loading both regular and fast cadence data for comparison")
        print(f"   Regular files: {len(regular_fits)}")
        print(f"   Fast files: {len(fast_fits)}")
        fits_files_to_use.extend(regular_fits)
        fits_files_to_use.extend(fast_fits)
    elif fast == False:
        print(f"\n📊 Using regular cadence data only ({len(regular_fits)} files)")
        fits_files_to_use.extend(regular_fits)
    elif fast == True:
        print(f"\n📊 Using fast cadence data only ({len(fast_fits)} files)")
        fits_files_to_use.extend(fast_fits)
    
    # Sort the final list
    fits_files_to_use = sorted(fits_files_to_use)
    
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
        page_dlogz_values = dlogz_values[start_flare:end_flare] if dlogz_values is not None else None
        page_flare_lengths = flare_lengths[start_flare:end_flare] if flare_lengths is not None else None
        page_window_points = window_points_list[start_flare:end_flare]  # Get window points for this page
        
        # Calculate grid layout for this page
        n_cols = min(3, page_flares)  # Maximum 3 columns
        n_rows = int(np.ceil(page_flares / n_cols))
        
        # Create figure with dynamic size based on number of rows
        # Constrained to fit on 7x11 inch page with 0.5 inch margins
        fig_width = 7.2 # 7 - 2*0.5 = 6 inches available width
        # Calculate height: base height per row + space for title and axis labels
        height_per_row = 1.5  # Height per row of subplots (reduced for more compact plots)
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
            current_window_points = page_window_points[i]  # Get window size for this specific flare
        
            # Filter data based on cadence_list if provided
            if cadence_list is not None:
                cadence_type = cadence_list[start_flare + i]
                if cadence_type == 'r':
                    # Use only regular cadence data
                    mask = ~all_is_fast
                    filtered_time = all_time[mask]
                    filtered_flux = all_flux[mask]
                    filtered_error = all_error[mask]
                    print(f"Flare {global_flare_num}: Using regular cadence data ({len(filtered_time)} points)")
                else:  # cadence_type == 'f'
                    # Use only fast cadence data
                    mask = all_is_fast
                    filtered_time = all_time[mask]
                    filtered_flux = all_flux[mask]
                    filtered_error = all_error[mask]
                    print(f"Flare {global_flare_num}: Using fast cadence data ({len(filtered_time)} points)")
            else:
                filtered_time = all_time
                filtered_flux = all_flux
                filtered_error = all_error
                print(f"Flare {global_flare_num}: Using all data ({len(filtered_time)} points)")
            
            # Find data points closest to the flare epoch in filtered data
            center_idx = find_nearest(filtered_time, epoch)
            # Define window around flare using the flare-specific window size
            start_idx = max(0, center_idx - current_window_points)
            end_idx = min(len(filtered_time), center_idx + current_window_points + 1)
            
            # Extract window data
            time_window = filtered_time[start_idx:end_idx]
            flux_window = filtered_flux[start_idx:end_idx]
            error_window = filtered_error[start_idx:end_idx]
            is_fast_window = all_is_fast[start_idx:end_idx]

            # Local median baseline using a wider ±100-point window, excluding the central flare window
            baseline_window_points = 100
            base_start_idx = max(0, center_idx - baseline_window_points)
            base_end_idx = min(len(filtered_time), center_idx + baseline_window_points + 1)
            base_indices = np.arange(base_start_idx, base_end_idx)
            # Exclude the central plotting window to avoid flare bias
            exclusion_mask = (base_indices < start_idx) | (base_indices >= end_idx)
            base_flux = filtered_flux[base_start_idx:base_end_idx][exclusion_mask]
            if len(base_flux) > 0:
                local_median = np.nanmedian(base_flux)
            else:
                local_median = np.nan
            
            # Plot light curve - separate by cadence type if both are present
            if plot_both_cadences and cadence_list is None:
                # Plot regular cadence data in black
                regular_mask = ~is_fast_window
                if np.any(regular_mask):
                    ax.errorbar(time_window[regular_mask], flux_window[regular_mask], 
                               yerr=error_window[regular_mask], 
                               fmt='-o', markersize=3, color='black', alpha=0.7, 
                               elinewidth=0.5, capsize=0, label='Regular')
                
                # Plot fast cadence data in blue
                fast_mask = is_fast_window
                if np.any(fast_mask):
                    ax.errorbar(time_window[fast_mask], flux_window[fast_mask], 
                               yerr=error_window[fast_mask], 
                               fmt='-o', markersize=3, color='#1E90FF', alpha=0.7, 
                               elinewidth=0.5, capsize=0, label='Fast')
            else:
                # Plot all data in single color (when using cadence_list, data is pre-filtered)
                ax.errorbar(time_window, flux_window, yerr=error_window, 
                           fmt='-o', markersize=3, color='black', alpha=0.7, 
                           elinewidth=0.5, capsize=0, label='Data')

            # Plot dashed line at local median baseline
            if np.isfinite(local_median):
                ax.axhline(local_median, color='gray', linestyle='--', linewidth=0.8, alpha=0.8, label='Local Median')
            
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
            
            ax.plot(filtered_time[peak_idx], center_flux, marker='D', markersize=3, 
                    color='orange', label='Flare Peak', zorder=5)
            
            # Use amplitude from catalog if available, otherwise calculate from flux
            if page_amplitudes is not None:
                flare_amplitude = page_amplitudes[i]
            else:
                flare_amplitude = center_flux - 1.0
            
            # Highlight flare decay if flare_lengths is provided
            if page_flare_lengths is not None and page_flare_lengths[i] > 0:
                decay_length = page_flare_lengths[i]
                # Start from the peak index and extend for decay_length points
                decay_start_idx = peak_idx + 1
                decay_end_idx = min(peak_idx + decay_length + 1, len(filtered_time))
                
                if decay_end_idx > decay_start_idx:
                    decay_time = filtered_time[decay_start_idx:decay_end_idx]
                    decay_flux = filtered_flux[decay_start_idx:decay_end_idx]
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
            time_text = f'$t_\\text{{mid}} = {center_time:.2f} [\\text{{d}}]$'
            amp_text = f'$A_{{\\text{{peak}}}} = {flare_amplitude:.3f}$'
            dlogz_text = f'$\\Delta \\log{{Z}} = {page_dlogz_values[i]:.3f}$' if page_dlogz_values is not None else ''
            
            # Build annotation text based on available data
            annotation_parts = [time_text, amp_text, dlogz_text]
            if page_fwhm_values is not None:
                fwhm_value = page_fwhm_values[i]
                fwhm_text = f'$t_\\text{{FWHM}} = {1440*fwhm_value:.3f} \\, \\mathrm{{[min]}}$'
                annotation_parts.append(fwhm_text)
            
            ax.text(0.02, 0.95, '\n'.join(annotation_parts), 
                   transform=ax.transAxes, fontsize=7.5, 
                   verticalalignment='top', horizontalalignment='left',
                   bbox=dict(boxstyle='round', facecolor='white', alpha=0.7, 
                            edgecolor='gray', linewidth=0.5),
                   family='monospace')
            
            # Formatting - add ticks pointing inward
            ax.tick_params(direction='in', length=3, width=0.5, colors='black')
            ax.set_xticklabels([])  # No tick labels
            ax.set_yticklabels([])  # No tick labels
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
        fig.text(0.5, 0.05, 'Time (BJD - 2457000)', ha='center', fontsize=11, weight='bold')
        fig.text(0, 0.5, 'Rel. Flux', va='center', rotation='vertical', 
                fontsize=11, weight='bold')
        
        # Add title with page number if multiple pages
        if plot_title is not None:
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


def plot_long_timeseries_with_gaps(csv_file, host_id, fits_dir, gap_days=1.0,
                                   max_lines_per_page=8, max_points_per_line=400,
                                   apply_quality_filter=True, fast=False,
                                   output_file=None, plot_title=None,
                                   highlight_flares=True):
    """Plot long-form time series split into lines at day-scale gaps.

    Parameters
    ----------
    csv_file : str
        Path to flare catalog (needs Host_ID and Epoch_BJD columns) for marking flares.
    host_id : str
        Host/star name as it appears in the catalog and FITS folder name.
    fits_dir : str
        Directory containing the FITS light curve files for the host. All .fits files
        in this directory (and subdirectories) are ingested, similar to plot_flares_from_catalog.
    gap_days : float, optional
        Threshold (in days) to start a new line when time gaps exceed this value. Default 1.0.
    max_lines_per_page : int, optional
        Maximum number of lines (segments) to show per page. Default 8.
    max_points_per_line : int, optional
        If a continuous segment exceeds this many points, it is split into multiple
        contiguous lines (no thinning/decimation). Default 400.
    apply_quality_filter : bool, optional
        Apply TESS quality filtering (quality == 0). Default True.
    fast : bool, optional
        If True, use fast cadence files (a_fast-lc.fits). If False, use regular cadence files. Default False.
    output_file : str or None, optional
        If provided, save plots to this path (adds _pageN for multi-page). Default None (show only).
    plot_title : str or None, optional
        Title for the figure. Default "{host_id} - Long Time Series".
    highlight_flares : bool, optional
        If True, draw vertical highlights at flare epochs (if present). Default True.

    Returns
    -------
    figures : list of matplotlib.figure.Figure
    all_axes : list of numpy.ndarray
    """

    # Read catalog to find flares for this host
    catalog = pd.read_csv(csv_file)
    host_flares = catalog[catalog['Host_ID'] == host_id].copy()
    flare_epochs = None
    if len(host_flares) > 0 and 'Epoch_BJD' in host_flares.columns:
        flare_epochs = host_flares['Epoch_BJD'].values - 2457000.0  # convert to BTJD like FITS

    # Discover FITS files
    all_fits_files = sorted(glob.glob(os.path.join(fits_dir, '**', '*.fits'), recursive=True))
    if len(all_fits_files) == 0:
        raise ValueError(f"No FITS files found in directory: {fits_dir}")

    if fast:
        fits_files_to_use = [f for f in all_fits_files if f.endswith('a_fast-lc.fits')]
        cadence_label = "fast"
    else:
        fits_files_to_use = [f for f in all_fits_files if f.endswith('_lc.fits') and not f.endswith('a_fast-lc.fits')]
        cadence_label = "regular"

    if len(fits_files_to_use) == 0:
        raise ValueError(f"No {cadence_label} cadence FITS files found in {fits_dir}")

    # Load light curve data
    time_all = []
    flux_all = []
    err_all = []

    print(f"\n📊 Loading {len(fits_files_to_use)} {cadence_label} cadence files for {host_id}")
    for fits_file in fits_files_to_use:
        try:
            lc = TESS_data_extract_local(fits_file, PDCSAP_ERR=True, apply_quality_filter=apply_quality_filter)
            time_all.extend(lc.time)
            flux_all.extend(lc.flux)
            err_all.extend(lc.error)
        except Exception as e:
            print(f"Warning: Could not read {fits_file}: {e}")
            continue

    if len(time_all) == 0:
        raise ValueError("No valid data could be read from FITS files.")

    time_all = np.array(time_all)
    flux_all = np.array(flux_all)
    err_all = np.array(err_all)

    # Sort by time
    order = np.argsort(time_all)
    time_all = time_all[order]
    flux_all = flux_all[order]
    err_all = err_all[order]

    # Split into segments when gaps exceed gap_days
    gaps = np.where(np.diff(time_all) > gap_days)[0]
    breakpoints = list(gaps + 1)
    indices = [0] + breakpoints + [len(time_all)]
    segments = []
    for start, end in zip(indices[:-1], indices[1:]):
        seg_time = time_all[start:end]
        seg_flux = flux_all[start:end]
        seg_err = err_all[start:end]
        if len(seg_time) == 0:
            continue

        # Split long continuous segments into multiple lines without thinning
        if max_points_per_line is not None and max_points_per_line > 0 and len(seg_time) > max_points_per_line:
            num_chunks = int(np.ceil(len(seg_time) / max_points_per_line))
            threshold = 0.25 * max_points_per_line
            for chunk_idx in range(num_chunks):
                c_start = chunk_idx * max_points_per_line
                c_end = min((chunk_idx + 1) * max_points_per_line, len(seg_time))
                chunk_size = c_end - c_start
                
                # If this is the last chunk and it's too small, merge with previous segment
                if chunk_idx == num_chunks - 1 and chunk_size < threshold and len(segments) > 0:
                    # Merge with the last segment
                    prev_time, prev_flux, prev_err = segments[-1]
                    merged_time = np.concatenate([prev_time, seg_time[c_start:c_end]])
                    merged_flux = np.concatenate([prev_flux, seg_flux[c_start:c_end]])
                    merged_err = np.concatenate([prev_err, seg_err[c_start:c_end]])
                    segments[-1] = (merged_time, merged_flux, merged_err)
                else:
                    segments.append((seg_time[c_start:c_end], seg_flux[c_start:c_end], seg_err[c_start:c_end]))
        else:
            segments.append((seg_time, seg_flux, seg_err))

    if len(segments) == 0:
        raise ValueError("No segments could be formed from the light curve data.")

    # Pagination
    n_pages = int(np.ceil(len(segments) / max_lines_per_page))
    if plot_title is None:
        plot_title = f"{host_id} - Long Time Series"

    figures = []
    all_axes = []

    for page in range(n_pages):
        start_idx = page * max_lines_per_page
        end_idx = min((page + 1) * max_lines_per_page, len(segments))
        page_segments = segments[start_idx:end_idx]

        n_rows = len(page_segments)
        fig_width = 6.0
        height_per_row = 2.0
        base_height = 1.0
        fig_height = max(2.5, min(10.0, n_rows * height_per_row + base_height))

        fig, axes = plt.subplots(n_rows, 1, figsize=(fig_width, fig_height), sharex=False)
        if n_rows == 1:
            axes = np.array([axes])

        plt.subplots_adjust(left=0.08, right=0.97, bottom=0.08, top=0.95, hspace=0)

        for idx, (ax, segment) in enumerate(zip(axes, page_segments)):
            seg_time, seg_flux, seg_err = segment

            ax.errorbar(seg_time, seg_flux, yerr=seg_err, fmt='-o', markersize=2.0,
                        color='black', alpha=0.8, elinewidth=0.5, capsize=0)

            if highlight_flares and flare_epochs is not None:
                in_seg = (flare_epochs >= seg_time.min()) & (flare_epochs <= seg_time.max())
                for fe in flare_epochs[in_seg]:
                    ax.axvspan(fe, fe+0.01, color='orange', alpha=0.3, zorder=2)

            ax.grid(True, alpha=0.3, linestyle='--', linewidth=0.5)
            ax.set_ylabel('Rel. Flux', fontsize=9)
            ax.tick_params(direction='in', length=3, width=0.5, colors='black')

            # Remove x-tick labels from all but the last subplot
            if idx < len(page_segments) - 1:
                ax.set_xticklabels([])

            # Add segment label on the right
            ax.text(1.01, 0.5, f'Segment {start_idx + idx + 1}', transform=ax.transAxes,
                    va='center', ha='left', fontsize=8)

        # Shared x-label
        axes[-1].set_xlabel('Time (BJD - 2457000)', fontsize=10)

        # Page title
        if n_pages > 1:
            title = f"{plot_title} (Page {page + 1}/{n_pages})"
        else:
            title = plot_title
        fig.suptitle(title, fontsize=12, weight='bold', y=0.995)

        if output_file:
            base, ext = os.path.splitext(output_file)
            page_file = f"{base}_page{page + 1}{ext}" if n_pages > 1 else output_file
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
