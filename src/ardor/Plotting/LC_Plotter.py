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
from ardor.Flares.Flare import tier0, tier1, tier2
import pandas as pd
# Configure Latin Modern Roman font
font_path = '/ugrad/whitsett.n/fonts/latin-modern-roman/lmroman10-regular.otf'
font_manager.fontManager.addfont(font_path)
plt.rcParams['font.family'] = 'Latin Modern Roman'

def plot_time_slices(time, flux, indices, durations, window_size=25, max_subplots=None, 
                    figsize=(7, 5), highlight_color='red', highlight_marker='o', 
                    normal_color='blue', normal_marker='.', show_subplot_titles=False):
    """
    Plot multiple time slices of light curve data with highlighted regions.
    
    Parameters:
    -----------
    time : array-like
        Time array (x-axis data)
    flux : array-like
        Flux array (y-axis data)
    indices : array-like
        List of indices along the time array to center slices on
    durations : array-like
        Array of durations for highlighting (same length as indices)
    window_size : int, optional
        Total number of points to show around each center point (default: 25)
    max_subplots : int, optional
        Maximum number of subplots to create. If None, uses all indices
    figsize : tuple, optional
        Figure size (width, height) in inches
    highlight_color : str, optional
        Color for highlighted data points
    highlight_marker : str, optional
        Marker style for highlighted data points
    normal_color : str, optional
        Color for normal data points
    normal_marker : str, optional
        Marker style for normal data points
    show_subplot_titles : bool, optional
        Whether to show individual subplot titles (default: False)
    
    Returns:
    --------
    fig : matplotlib.figure.Figure
        The created figure object
    axes : array
        Array of subplot axes
    """
    
    # Convert inputs to numpy arrays
    time = np.array(time)
    flux = np.array(flux)
    indices = np.array(indices)
    durations = np.array(durations)
    
    # Validate input lengths
    if len(indices) != len(durations):
        raise ValueError("indices and durations arrays must have the same length")
    
    # Determine number of subplots
    if max_subplots is None:
        n_subplots = len(indices)
    else:
        n_subplots = min(max_subplots, len(indices))
    
    # Calculate subplot grid dimensions
    n_cols = int(np.ceil(np.sqrt(n_subplots)))
    n_rows = int(np.ceil(n_subplots / n_cols))
    
    # Create figure and subplots
    fig, axes = plt.subplots(n_rows, n_cols, figsize=figsize)
    if n_subplots == 1:
        axes = [axes]
    elif n_rows == 1 or n_cols == 1:
        axes = axes.flatten()
    else:
        axes = axes.flatten()
    
    half_window = window_size // 2
    
    for i in range(n_subplots):
        ax = axes[i]
        center_idx = indices[i]
        duration = durations[i]
        
        # Calculate slice boundaries
        start_idx = max(0, center_idx - half_window)
        end_idx = min(len(time), center_idx + half_window + 1)
        
        # Extract slice data
        time_slice = time[start_idx:end_idx]
        flux_slice = flux[start_idx:end_idx]
        
        # Plot normal data points
        ax.plot(time_slice, flux_slice, color=normal_color, marker=normal_marker, 
                linestyle='-', markersize=4, alpha=0.7, label='Data')
        
        # Determine highlight region based on number of data points
        center_time = time[center_idx]
        highlight_start_idx = center_idx
        highlight_end_idx = min(len(time), center_idx + int(duration))
        
        # Convert global indices to slice-relative indices
        slice_start_in_global = start_idx
        slice_end_in_global = end_idx
        
        # Find which part of the highlight region overlaps with this slice
        highlight_start_in_slice = max(0, highlight_start_idx - slice_start_in_global)
        highlight_end_in_slice = min(len(time_slice), highlight_end_idx - slice_start_in_global)
        
        # Create highlight mask for the slice
        if highlight_start_in_slice < highlight_end_in_slice and highlight_start_in_slice >= 0:
            highlight_mask = np.zeros(len(time_slice), dtype=bool)
            highlight_mask[highlight_start_in_slice:highlight_end_in_slice] = True
        else:
            highlight_mask = np.zeros(len(time_slice), dtype=bool)
        
        if np.any(highlight_mask):
            ax.plot(time_slice[highlight_mask], flux_slice[highlight_mask], 
                   color=highlight_color, marker=highlight_marker, linestyle='None',
                   markersize=6, label='Highlighted region')
        
        # Highlight the center point specifically
        center_time_in_slice = center_time
        if center_time_in_slice in time_slice:
            center_flux = flux[center_idx]
            ax.plot(center_time_in_slice, center_flux, color='orange', marker='*', 
                   markersize=10, label='Center point')
        
        # Formatting
        if show_subplot_titles:
            ax.set_title(f'Slice {i+1} (Center: t={center_time:.3f})')
        ax.grid(True, alpha=0.3)
        
        # Add legend only to first subplot to avoid clutter
        if i == 0:
            ax.legend(fontsize=8)
        
        # Set smaller tick label font sizes
        ax.tick_params(labelsize=8)
    
    # Hide unused subplots
    for j in range(n_subplots, len(axes)):
        axes[j].set_visible(False)
    
    # Add unified, centered, and larger axis labels with better positioning
    fig.supxlabel('Time (BJD)', fontsize=12, y=0.01)
    fig.supylabel('Normalized Flux', fontsize=16, x=0.01)
    
    # Minimize whitespace
    plt.tight_layout()
    plt.subplots_adjust(left=0.08, bottom=0.08, right=0.98, top=0.95, hspace=0.25, wspace=0.20)
    
    return fig, axes[:n_subplots]

def Epoch_Extractor(host_name, csv_path):
    """
    Extract epochs from a CSV file for a given host name.
    
    Parameters:
    -----------
    host_name : str
        The name of the host to extract epochs for.
    csv_path : str
        Path to the CSV file containing epoch data.
    
    Returns:
    --------
    epochs : list
        List of epochs corresponding to the host name.
    """
    df = pd.read_csv(csv_path)
    host_data = df[df['Host_ID'] == host_name]
    epochs = host_data['Flare_Epoch'].tolist()
    return epochs
csv = '/ugrad/whitsett.n/Induced_Flares/Flare_Data/Tier_3/Grand_Lists/All_Exoplanet_MCMC_Flares.csv'
epochs = Epoch_Extractor('TOI-1062', csv)
print(epochs)
files = os.listdir('/data2/whitsett.n/TESS/Hosts/TOI-1062')
print(len(files))
data = os.path.join('/data2/whitsett.n/TESS/Hosts/TOI-1062', files[8])
lc = tier0(data)

# Convert astropy MaskedQuantity to numpy array for compatibility
detrended_flux_array = np.array(lc.detrended_flux.value)
time_array = np.array(lc.time)

flares = tier1(time_array, detrended_flux_array, sigma=3)
fig, axes = plot_time_slices(lc.time, lc.flux, flares.index, flares.length, window_size=100, max_subplots=10)
plt.show()
# %%
