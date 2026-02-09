# ============================================================================
# NTA PLOTTING FUNCTIONS - ORIGINAL CODE FROM JUPYTER NOTEBOOK
# ============================================================================
# This is the ACTUAL plotting code that generates publication-quality plots
# with lognormal fits, error bars, D-value lines, and confidence intervals
# ============================================================================

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.lines import Line2D
from matplotlib.ticker import LogLocator, LogFormatter
from scipy import stats as scipy_stats
from scipy.optimize import curve_fit



# ======================================================================
# CELL 9
# ======================================================================





def lognormal_pdf(x, mu, sigma, amplitude):
    """Calculate lognormal probability density function."""
    return amplitude * (1 / (x * sigma * np.sqrt(2 * np.pi))) * np.exp(-(np.log(x) - mu)**2 / (2 * sigma**2))


def fit_lognormal_distribution(sizes, weights):
    """Fit a lognormal distribution to size distribution data."""
    try:
        # Initial parameter estimates
        size_log = np.log(sizes)
        initial_mu = np.average(size_log, weights=weights)
        initial_sigma = np.sqrt(np.average((size_log - initial_mu)**2, weights=weights))
        initial_amplitude = np.max(weights) * initial_sigma * np.sqrt(2 * np.pi) * np.exp(initial_mu)
        
        initial_params = [initial_mu, initial_sigma, initial_amplitude]
        
        # Perform curve fitting
        params, _ = curve_fit(
            lognormal_pdf, sizes, weights, p0=initial_params,
            bounds=([0, 0, 0], [np.inf, np.inf, np.inf]), maxfev=10000
        )
        
        # Generate fitted curve
        size_range = np.linspace(sizes.min(), sizes.max(), 200)
        fitted_curve = lognormal_pdf(size_range, *params)
        
        return True, (size_range, fitted_curve, params)
        
    except Exception as e:
        return False, f"Lognormal fit failed: {str(e)}"


def add_number_fit_curve_fixed(ax, plot_df, is_log_scale, fit_color='#F25C54'):
    """
    Add lognormal fit curve to number distribution plot.
    Calculates independent fits for each scale.
    """
    fit_legend_elements = []
    
    sizes = plot_df['size_nm'].values
    weights = plot_df['number_normalized_avg'].values
    
    # Remove any zero or negative weights
    valid_mask = (weights > 0) & (sizes > 0)
    if not np.any(valid_mask):
        return fit_legend_elements, None
    
    sizes_valid = sizes[valid_mask]
    weights_valid = weights[valid_mask]
    
    # Calculate fit independently for this scale
    success, result = fit_lognormal_distribution(sizes_valid, weights_valid)
    if not success:
        return fit_legend_elements, None
    
    size_range, fitted_curve, params = result
    
    # Plot the fit curve
    ax.plot(size_range, fitted_curve, '-', color=fit_color, linewidth=2.5, 
           alpha=0.9, label='Lognormal Fit', zorder=4)
    
    # Calculate geometric statistics
    geometric_mean = np.exp(params[0])
    geometric_std = np.exp(params[1])
    
    # Add fit info to legend
    fit_legend_elements.append(
        Line2D([0], [0], color=fit_color, linestyle='-', linewidth=2.5,
              label=f'Lognormal: μ={geometric_mean:.1f} nm, σ={geometric_std:.2f}')
    )
    
    return fit_legend_elements, ('lognormal', {
        'mu': params[0], 
        'sigma': params[1], 
        'amplitude': params[2],
        'geo_mean': geometric_mean, 
        'geo_std': geometric_std
    })


def add_d_value_lines_and_bands(ax, stats):
    """Add D-value lines and uncertainty bands to a subplot."""
    legend_elements = []
    
    if not stats or 'D10_avg' not in stats:
        return legend_elements
    
    d10_avg = stats['D10_avg']
    d10_lower = stats.get('D10_lower', d10_avg)
    d10_upper = stats.get('D10_upper', d10_avg)
    
    d50_avg = stats['D50_avg'] 
    d50_lower = stats.get('D50_lower', d50_avg)
    d50_upper = stats.get('D50_upper', d50_avg)
    
    d90_avg = stats['D90_avg']
    d90_lower = stats.get('D90_lower', d90_avg)
    d90_upper = stats.get('D90_upper', d90_avg)
    
    span = stats.get('span_avg', (d90_avg-d10_avg)/d50_avg if d50_avg > 0 else 0)
    
    # Add D-value lines and bands
    for d_val, d_lower, d_upper, style, width, alpha_band in [
        (d10_avg, d10_lower, d10_upper, '--', 1.5, 0.15),
        (d50_avg, d50_lower, d50_upper, '-', 2.5, 0.25), 
        (d90_avg, d90_lower, d90_upper, '--', 1.5, 0.15)
    ]:
        if not np.isnan(d_val):
            ax.axvline(x=d_val, color='gray', linestyle=style, alpha=0.8, linewidth=width, zorder=5)
            if not np.isnan(d_lower) and not np.isnan(d_upper) and (d_lower != d_val or d_upper != d_val):
                ax.axvspan(d_lower, d_upper, alpha=alpha_band, color='gray', zorder=1)
    
    # Create legend elements
    legend_elements.extend([
        Line2D([0], [0], color='gray', linestyle='--', linewidth=1.5, 
              label=f'D10: {d10_avg:.1f} nm ({d10_lower:.1f}-{d10_upper:.1f})'),
        Line2D([0], [0], color='gray', linestyle='-', linewidth=2.5, 
              label=f'D50: {d50_avg:.1f} nm ({d50_lower:.1f}-{d50_upper:.1f})'),
        Line2D([0], [0], color='gray', linestyle='--', linewidth=1.5, 
              label=f'D90: {d90_avg:.1f} nm ({d90_lower:.1f}-{d90_upper:.1f})'),
        Line2D([0], [0], color='white', linestyle='', 
              label=f'Span: {span:.3f}')
    ])
    
    return legend_elements


def create_number_plot(plot_df, is_log_scale, stats=None, uniqueID=None, metadata=None):
    """Create a two-subplot plot for number-weighted distribution."""
    
    scale_name = "Logarithmic" if is_log_scale else "Linear"
    xscale = 'log' if is_log_scale else 'linear'
    color = '#4C5B5C'  # Slate gray for number-weighted data
    
    # Sort by size
    plot_df = plot_df.sort_values('size_nm')
    
    # Create figure
    fig = plt.figure(figsize=(7, 9))
    gs = gridspec.GridSpec(2, 1, height_ratios=[0.6, 0.4], hspace=0.3, 
                          top=0.82, bottom=0.08)
    
    # =================================================================
    # TOP SUBPLOT: MAIN DISTRIBUTION WITH ERROR BARS AND FITS
    # =================================================================
    ax1 = fig.add_subplot(gs[0])
    
    # Plot main distribution with error bars
    if 'number_normalized_sd' in plot_df.columns:
        ax1.errorbar(plot_df['size_nm'], plot_df['number_normalized_avg'], 
                    yerr=plot_df['number_normalized_sd'],
                    fmt='o', color=color, ecolor=color, alpha=0.7,
                    capsize=3, capthick=1, markersize=6, linewidth=1.5,
                    label='Number Distribution')
    else:
        ax1.scatter(plot_df['size_nm'], plot_df['number_normalized_avg'], 
                   color=color, s=60, alpha=0.8, label='Number Distribution')
    
    # Add lognormal fit curve (calculated independently for this scale)
    fit_result = add_number_fit_curve_fixed(ax1, plot_df, is_log_scale)
    if isinstance(fit_result, tuple):
        fit_legend_elements, fit_results = fit_result
    else:
        fit_legend_elements = fit_result
        fit_results = None
    
    # Format top subplot
    ax1.set_ylabel('Normalized Number', color=color, fontsize=14, labelpad=10)
    ax1.tick_params(axis='y', labelcolor=color, labelsize=12)
    ax1.tick_params(axis='x', labelsize=12)
    ax1.spines['left'].set_color(color)
    
    # Set x-axis scale and limits
    ax1.set_xscale(xscale)
    if is_log_scale:
        # Log scale: start at 30 nm, focus on signal range
        min_size = max(30, plot_df['size_nm'].min())
        percentile_90 = np.percentile(plot_df['size_nm'], 90)
        max_size = min(percentile_90 * 1.3, 300)
        ax1.set_xlim([min_size, max_size])
    else:
        # Linear scale: start at 0, focus on main signal
        percentile_85 = np.percentile(plot_df['size_nm'], 85)
        max_size = min(percentile_85 * 1.2, 250)
        ax1.set_xlim([0, max_size])
    
    # Set y-axis to start from 0
    y_min, y_max = ax1.get_ylim()
    ax1.set_ylim([0, y_max])
    
    # Add D-value lines and bands
    d_legend_elements = add_d_value_lines_and_bands(ax1, stats)
    
    # Create comprehensive legend for top plot - PLACE OUTSIDE
    main_legend = [Line2D([0], [0], marker='o', color='w', markerfacecolor=color, 
                          markersize=8, label='Number Distribution')]
    
    all_legend_elements = main_legend + fit_legend_elements + d_legend_elements
    leg1 = ax1.legend(handles=all_legend_elements, fontsize=9, frameon=True, 
                     bbox_to_anchor=(1.05, 1), loc='upper left')
    leg1.get_frame().set_alpha(0.95)
    leg1.get_frame().set_edgecolor('lightgray')
    
    ax1.grid(True, linestyle='--', alpha=0.4)
    ax1.set_xlabel('')  # No x-label on top plot
    
    # =================================================================
    # BOTTOM SUBPLOT: CUMULATIVE DISTRIBUTION WITH UNCERTAINTY BANDS
    # =================================================================
    ax2 = fig.add_subplot(gs[1])
    
    if 'number_normalized_cumsum_avg' in plot_df.columns:
        # Plot cumulative curve as percentage
        cumsum_values = plot_df['number_normalized_cumsum_avg'] * 100
        ax2.plot(plot_df['size_nm'], cumsum_values, '-', 
                color=color, linewidth=3, alpha=0.9, label='Cumulative %')
        
        # Add uncertainty bands if available
        if 'number_normalized_cumsum_sd' in plot_df.columns:
            cumsum_sd = plot_df['number_normalized_cumsum_sd'] * 100
            ax2.fill_between(plot_df['size_nm'], 
                           cumsum_values - cumsum_sd,
                           cumsum_values + cumsum_sd,
                           color=color, alpha=0.25, zorder=1, label='± SD')
        
        ax2.set_ylim([0, 110])
        ax2.set_ylabel('Cumulative Percentage (%)', color=color, fontsize=14, labelpad=10)
        ax2.tick_params(axis='y', labelcolor=color, labelsize=12)
        ax2.spines['left'].set_color(color)
    
    # Format bottom subplot
    ax2.set_xlabel('Size (nm)', fontsize=14, labelpad=10)
    ax2.tick_params(axis='x', labelsize=12)
    ax2.set_xscale(xscale)
    ax2.set_xlim(ax1.get_xlim())  # Match top plot limits
    
    # Add D-value lines to bottom plot
    if 'number_normalized_cumsum_avg' in plot_df.columns:
        add_d_value_lines_and_bands(ax2, stats)
    
    # Legend for bottom plot - PLACE OUTSIDE
    cumulative_legend = [
        Line2D([0], [0], color=color, linewidth=3, label='Cumulative %'),
        Line2D([0], [0], color=color, alpha=0.25, linewidth=8, label='± SD')
    ]
    
    leg2 = ax2.legend(handles=cumulative_legend, fontsize=10, frameon=True, 
                     bbox_to_anchor=(1.05, 1), loc='upper left')
    leg2.get_frame().set_alpha(0.95)
    leg2.get_frame().set_edgecolor('lightgray')
    
    ax2.grid(True, linestyle='--', alpha=0.4)
    
    # =================================================================
    # TITLE AND METADATA
    # =================================================================
    
    # Extract replicate info
    replicate_info = ""
    if metadata and 'num_replicates' in metadata:
        num_reps = metadata['num_replicates']
        if num_reps and str(num_reps) != '1':
            replicate_info = f" (n={num_reps})"
    
    # Set main title
    main_title = f'{scale_name} Number-Weighted\nDistribution: {uniqueID}{replicate_info}'
    fig.suptitle(main_title, fontsize=14, fontweight='bold', y=0.94)
    
    # Add subtitle
    subtitle = f"Error bars/bands: ± SD | Fits: Lognormal"
    fig.text(0.5, 0.87, subtitle, ha='center', fontsize=11, style='italic')
    
    return fig, fit_results


def generate_number_plots(distribution_df, stats_dict=None, uniqueID=None, 
                         metadata=None, output_dir=None, config=None):
    """Generate number-weighted distribution plots for both linear and log scales."""
    
    if distribution_df is None or distribution_df.empty:
        return False, "No data available for plotting"
    
    plt.style.use('default')
    
    # Determine output directory
    if output_dir is None:
        if config is not None and "directory" in config:
            base_dir = config["directory"]
            output_dir = os.path.join(base_dir, "processed")
        else:
            output_dir = os.path.join(os.getcwd(), "processed")
    
    try:
        os.makedirs(output_dir, exist_ok=True)
    except Exception as e:
        return False, f"Failed to create output directory: {str(e)}"
    
    created_files = []
    
    # Generate linear and logarithmic plots (each with independent fits)
    for is_log_scale in [False, True]:
        scale_type = 'logarithmic' if is_log_scale else 'linear'
        scale_name = 'log' if is_log_scale else 'linear'
        
        print(f"Creating {scale_name} number-weighted plot...")
        
        # Filter data for this scale
        plot_df = distribution_df[distribution_df['scale'] == scale_type].copy()
        
        if plot_df.empty:
            print(f"  Warning: No {scale_type} scale data available")
            continue
        
        # Get statistics
        stats = None
        if stats_dict and scale_type in stats_dict and 'number' in stats_dict[scale_type]:
            stats = stats_dict[scale_type]['number']
        
        # Create the plot (fits calculated independently for each scale)
        fig, fit_results = create_number_plot(plot_df, is_log_scale, stats, uniqueID, metadata)
        
        if fig is None:
            print(f"  Failed to create plot")
            continue
        
        # Save fit results to comprehensive fits file
        if fit_results:
            fit_type, fit_data = fit_results
            try:
                
                # Load existing fits file or create new one
                fits_filename = f"Fits_{uniqueID}_all.json"
                fits_path = os.path.join(output_dir, fits_filename)
                
                if os.path.exists(fits_path):
                    with open(fits_path, 'r') as f:
                        all_fits = json.load(f)
                else:
                    all_fits = {'dataset': uniqueID, 'fits': {}}
                
                # Add this fit to the collection
                fit_key = f"number_{scale_name}"
                all_fits['fits'][fit_key] = {
                    'distribution_type': 'number',
                    'scale': scale_name,
                    'fit_type': fit_type,
                    'parameters': fit_data
                }
                
                # Save updated fits file
                with open(fits_path, 'w') as f:
                    json.dump(all_fits, f, indent=2, default=str)
                
                print(f"  ✓ Saved fit to: {fits_filename}")
                
            except Exception as e:
                print(f"  ⚠ Failed to save fit: {str(e)}")
        
        # Save the plot
        try:
            base_filename = f"Plot_{uniqueID}_{scale_name}_number"
            
            pdf_path = os.path.join(output_dir, f"{base_filename}.pdf")
            fig.savefig(pdf_path, bbox_inches='tight', dpi=300)
            
            png_path = os.path.join(output_dir, f"{base_filename}.png")
            fig.savefig(png_path, bbox_inches='tight', dpi=300)
            
            created_files.append(pdf_path)
            print(f"  ✓ Saved: {base_filename}.pdf/.png")
            
            plt.close(fig)
            
        except Exception as e:
            print(f"  ✗ Failed to save plot: {str(e)}")
            plt.close(fig)
            continue
    
    return True, created_files


# Execute if we have the required data
if 'current_distribution_df' in globals() and current_distribution_df is not None:
    print("=" * 80)
    print("GENERATING NUMBER-WEIGHTED DISTRIBUTION PLOTS (FIXED VERSION)")
    print("=" * 80)
    
    uniqueID = current_metadata.get('persistentID', 'unknown') if 'current_metadata' in globals() else 'unknown'
    stats = current_stats if 'current_stats' in globals() else None
    metadata = current_metadata if 'current_metadata' in globals() else None
    config = CONFIG if 'CONFIG' in globals() else None
    
    print(f"Creating number-weighted plots for: {uniqueID}")
    print("Includes: Linear + Logarithmic (Lognormal fits for both)")
    
    success, plot_files = generate_number_plots(
        current_distribution_df,
        stats_dict=stats,
        uniqueID=uniqueID,
        metadata=metadata,
        config=config
    )
    
    if not success:
        print(f"ERROR: {plot_files}")
    else:
        print(f"\n✓ Successfully created {len(plot_files)} number-weighted plots!")
        for filepath in plot_files:
            print(f"  - {os.path.basename(filepath)}")
        
        current_number_plots = plot_files

else:
    print("No data found. Run the complete workflow first.")


# ======================================================================
# CELL 10
# ======================================================================





def lognormal_pdf(x, mu, sigma, amplitude):
    """Calculate lognormal probability density function."""
    return amplitude * (1 / (x * sigma * np.sqrt(2 * np.pi))) * np.exp(-(np.log(x) - mu)**2 / (2 * sigma**2))



def fit_gaussian_mixture_direct(sizes, volumes, n_components_range=range(1, 4)):
    """
    Fit GMM directly to the distribution data without creating artificial datasets.
    
    Parameters:
    sizes (array): Size values (nm)
    volumes (array): Volume values (nm³/mL) 
    n_components_range: Range of component numbers to try
    
    Returns:
    success, (size_range, total_fit, individual_components, model_info)
    """
    try:
        # Remove zeros and sort by size
        valid_mask = (volumes > 0) & (sizes > 0)
        if not np.any(valid_mask):
            return False, "No valid data points"
        
        sizes_clean = sizes[valid_mask]
        volumes_clean = volumes[valid_mask]
        
        # Sort by size for consistency
        sort_idx = np.argsort(sizes_clean)
        sizes_clean = sizes_clean[sort_idx]
        volumes_clean = volumes_clean[sort_idx]
        
        if len(sizes_clean) < 6:  # Need minimum data points
            return False, "Insufficient data points for GMM"
        
        best_model = None
        best_aic = np.inf
        
        # Try different numbers of components
        for n_comp in n_components_range:
            try:
                # Initial parameter guess
                size_min, size_max = sizes_clean.min(), sizes_clean.max()
                initial_means = np.linspace(size_min + 0.1*(size_max-size_min), 
                                          size_max - 0.1*(size_max-size_min), n_comp)
                initial_stds = np.full(n_comp, (size_max - size_min) / (n_comp * 3))
                initial_weights = np.full(n_comp, 1.0 / n_comp)
                
                # Pack parameters: [means, stds, weights[:-1]]
                initial_params = np.concatenate([initial_means, initial_stds, initial_weights[:-1]])
                
                def unpack_params(params, n_comp):
                    means = params[:n_comp]
                    stds = params[n_comp:2*n_comp]
                    weights = np.zeros(n_comp)
                    weights[:-1] = params[2*n_comp:]
                    weights[-1] = 1.0 - np.sum(weights[:-1])  # Ensure sum = 1
                    return means, stds, weights
                
                def gmm_objective(params):
                    try:
                        means, stds, weights = unpack_params(params, n_comp)
                        
                        # Check constraints
                        if np.any(stds <= 0) or np.any(weights <= 0) or weights[-1] <= 0:
                            return 1e10
                        
                        # Calculate predicted volumes using GMM
                        predicted = np.zeros_like(sizes_clean)
                        for i in range(n_comp):
                            predicted += weights[i] * scipy_stats.norm.pdf(sizes_clean, means[i], stds[i])
                        
                        # Scale to match data magnitude
                        if np.max(predicted) > 0:
                            scale_factor = np.sum(volumes_clean) / np.sum(predicted)
                            predicted *= scale_factor
                        
                        # Calculate squared error
                        mse = np.mean((volumes_clean - predicted) ** 2)
                        return mse
                        
                    except:
                        return 1e10
                
                # Set bounds
                bounds = []
                # Bounds for means
                for _ in range(n_comp):
                    bounds.append((size_min, size_max))
                # Bounds for stds
                for _ in range(n_comp):
                    bounds.append((1.0, size_max - size_min))
                # Bounds for weights
                for _ in range(n_comp-1):
                    bounds.append((0.01, 0.98))
                
                # Optimize
                result = minimize(gmm_objective, initial_params, method='L-BFGS-B', 
                                bounds=bounds, options={'maxiter': 1000})
                
                if result.success:
                    final_means, final_stds, final_weights = unpack_params(result.x, n_comp)
                    
                    # Calculate final predicted curve
                    predicted_final = np.zeros_like(sizes_clean)
                    for i in range(n_comp):
                        predicted_final += final_weights[i] * scipy_stats.norm.pdf(sizes_clean, final_means[i], final_stds[i])
                    
                    scale_factor = np.sum(volumes_clean) / np.sum(predicted_final) if np.sum(predicted_final) > 0 else 1.0
                    
                    # Calculate AIC
                    mse = np.mean((volumes_clean - predicted_final * scale_factor) ** 2)
                    n_params = 3 * n_comp - 1
                    aic = len(sizes_clean) * np.log(mse) + 2 * n_params
                    
                    if aic < best_aic:
                        best_aic = aic
                        best_model = {
                            'n_components': n_comp,
                            'means': final_means,
                            'stds': final_stds,
                            'weights': final_weights,
                            'scale_factor': scale_factor,
                            'aic': aic,
                            'mse': mse
                        }
                        
            except Exception as e:
                print(f"    Failed to fit {n_comp} components: {e}")
                continue
        
        if best_model is None:
            return False, "No successful fits"
        
        # Generate smooth curves for plotting
        size_range = np.linspace(sizes_clean.min(), sizes_clean.max(), 200)
        
        # Calculate total fitted curve and individual components
        total_fit = np.zeros_like(size_range)
        individual_components = []
        
        for i in range(best_model['n_components']):
            # Individual component
            component = (best_model['weights'][i] * 
                        scipy_stats.norm.pdf(size_range, best_model['means'][i], best_model['stds'][i]) *
                        best_model['scale_factor'])
            individual_components.append(component)
            total_fit += component
        
        return True, (size_range, total_fit, individual_components, best_model)
        
    except Exception as e:
        return False, f"GMM fitting failed: {str(e)}"


def add_volume_fit_curve_fixed(ax, plot_df, is_log_scale, fit_color='#2C7F7F'):
    """Add lognormal fits for volume distributions. Calculates independent fits per scale."""
    fit_legend_elements = []
    
    sizes = plot_df['size_nm'].values
    volumes = plot_df['volume_nm^3_per_mL_avg'].values
    
    # Remove any zero or negative values
    valid_mask = (volumes > 0) & (sizes > 0)
    if not np.any(valid_mask):
        return fit_legend_elements, None
    
    sizes_valid = sizes[valid_mask]
    volumes_valid = volumes[valid_mask]
    
    # Calculate fit independently for this scale
    success, result = fit_lognormal_distribution(sizes_valid, volumes_valid)
    if not success:
        return fit_legend_elements, None
    
    size_range, fitted_curve, params = result
    
    # Plot the fit curve
    ax.plot(size_range, fitted_curve, '-', color=fit_color, linewidth=2.5, 
           alpha=0.9, label='Lognormal Fit', zorder=4)
    
    # Calculate geometric statistics
    geometric_mean = np.exp(params[0])
    geometric_std = np.exp(params[1])
    
    fit_legend_elements.append(
        Line2D([0], [0], color=fit_color, linestyle='-', linewidth=2.5,
              label=f'Lognormal: μ={geometric_mean:.1f} nm, σ={geometric_std:.2f}')
    )
    
    return fit_legend_elements, ('lognormal', {
        'mu': params[0], 
        'sigma': params[1], 
        'amplitude': params[2],
        'geo_mean': geometric_mean, 
        'geo_std': geometric_std
    })


def add_d_value_lines_and_bands(ax, stats):
    """Add D-value lines and uncertainty bands to a subplot."""
    legend_elements = []
    
    if not stats or 'D10_avg' not in stats:
        return legend_elements
    
    d10_avg = stats['D10_avg']
    d10_lower = stats.get('D10_lower', d10_avg)
    d10_upper = stats.get('D10_upper', d10_avg)
    
    d50_avg = stats['D50_avg'] 
    d50_lower = stats.get('D50_lower', d50_avg)
    d50_upper = stats.get('D50_upper', d50_avg)
    
    d90_avg = stats['D90_avg']
    d90_lower = stats.get('D90_lower', d90_avg)
    d90_upper = stats.get('D90_upper', d90_avg)
    
    span = stats.get('span_avg', (d90_avg-d10_avg)/d50_avg if d50_avg > 0 else 0)
    
    # Add D-value lines and bands
    for d_val, d_lower, d_upper, style, width, alpha_band in [
        (d10_avg, d10_lower, d10_upper, '--', 1.5, 0.15),
        (d50_avg, d50_lower, d50_upper, '-', 2.5, 0.25), 
        (d90_avg, d90_lower, d90_upper, '--', 1.5, 0.15)
    ]:
        if not np.isnan(d_val):
            ax.axvline(x=d_val, color='gray', linestyle=style, alpha=0.8, linewidth=width, zorder=5)
            if not np.isnan(d_lower) and not np.isnan(d_upper) and (d_lower != d_val or d_upper != d_val):
                ax.axvspan(d_lower, d_upper, alpha=alpha_band, color='gray', zorder=1)
    
    # Create legend elements
    legend_elements.extend([
        Line2D([0], [0], color='gray', linestyle='--', linewidth=1.5, 
              label=f'D10: {d10_avg:.1f} nm ({d10_lower:.1f}-{d10_upper:.1f})'),
        Line2D([0], [0], color='gray', linestyle='-', linewidth=2.5, 
              label=f'D50: {d50_avg:.1f} nm ({d50_lower:.1f}-{d50_upper:.1f})'),
        Line2D([0], [0], color='gray', linestyle='--', linewidth=1.5, 
              label=f'D90: {d90_avg:.1f} nm ({d90_lower:.1f}-{d90_upper:.1f})'),
        Line2D([0], [0], color='white', linestyle='', 
              label=f'Span: {span:.3f}')
    ])
    
    return legend_elements


def create_volume_plot(plot_df, is_log_scale, stats=None, uniqueID=None, metadata=None):
    """Create a two-subplot plot for volume-weighted distribution."""
    
    scale_name = "Logarithmic" if is_log_scale else "Linear"
    xscale = 'log' if is_log_scale else 'linear'
    color = '#2E7D32'  # Forest green for volume-weighted
    
    # Sort by size
    plot_df = plot_df.sort_values('size_nm')
    
    # Create figure
    fig = plt.figure(figsize=(7, 9))
    gs = gridspec.GridSpec(2, 1, height_ratios=[0.6, 0.4], hspace=0.3, 
                          top=0.82, bottom=0.08)
    
    # TOP SUBPLOT: MAIN DISTRIBUTION WITH ERROR BARS AND FITS
    ax1 = fig.add_subplot(gs[0])
    
    # Plot main distribution with error bars
    if 'volume_nm^3_per_mL_sd' in plot_df.columns:
        ax1.errorbar(plot_df['size_nm'], plot_df['volume_nm^3_per_mL_avg'], 
                    yerr=plot_df['volume_nm^3_per_mL_sd'],
                    fmt='o', color=color, ecolor=color, alpha=0.7,
                    capsize=3, capthick=1, markersize=6, linewidth=1.5,
                    label='Volume Distribution')
    else:
        ax1.scatter(plot_df['size_nm'], plot_df['volume_nm^3_per_mL_avg'], 
                   color=color, s=60, alpha=0.8, label='Volume Distribution')
    
    # Add fit curve and get fit results (calculated independently for this scale)
    fit_result = add_volume_fit_curve_fixed(ax1, plot_df, is_log_scale)
    if isinstance(fit_result, tuple):
        fit_legend_elements, fit_results = fit_result
    else:
        fit_legend_elements = fit_result
        fit_results = None
    
    # Format top subplot
    ax1.set_ylabel('Volume (nm³/mL)', color=color, fontsize=14, labelpad=10)
    ax1.tick_params(axis='y', labelcolor=color, labelsize=12)
    ax1.tick_params(axis='x', labelsize=12)
    ax1.spines['left'].set_color(color)
    
    # Set x-axis scale and add better tick labels for log scale
    ax1.set_xscale(xscale)
    if is_log_scale:
        # Add more detailed log scale ticks
        ax1.xaxis.set_major_locator(LogLocator(base=10, numticks=12))
        ax1.xaxis.set_minor_locator(LogLocator(base=10, subs=(0.2, 0.4, 0.6, 0.8), numticks=12))
        ax1.xaxis.set_major_formatter(LogFormatter(base=10, labelOnlyBase=False))
    
    # Smart volume-weighted range calculation (restored from original)
    weights_for_range = plot_df['volume_nm^3_per_mL_avg'].values
    sizes_for_range = plot_df['size_nm'].values
    
    # Find where 99% of the volume signal is contained
    cumsum_weights = np.cumsum(weights_for_range)
    total_weight = cumsum_weights[-1]
    
    # Find 1st and 99th percentiles of the volume-weighted distribution
    p1_idx = np.searchsorted(cumsum_weights, 0.01 * total_weight)
    p99_idx = np.searchsorted(cumsum_weights, 0.99 * total_weight)
    
    signal_min = sizes_for_range[max(0, p1_idx)]
    signal_max = sizes_for_range[min(len(sizes_for_range)-1, p99_idx)]
    data_max = plot_df['size_nm'].max()
    
    if is_log_scale:
        # Log scale: focus on the volume-weighted signal range with some padding
        min_size = max(signal_min * 0.7, 20)  # Don't go below 20 nm
        max_size = min(signal_max * 2.0, data_max * 1.2)  # Cap at reasonable range
        ax1.set_xlim([min_size, max_size])
    else:
        # Linear scale: tighten the range more, focus on main signal
        min_size = 0
        max_size = min(signal_max * 1.2, 400)  # Tighter range, cap at 400nm for most cases
        ax1.set_xlim([min_size, max_size])
    
    # Set y-axis to start from 0
    y_min, y_max = ax1.get_ylim()
    ax1.set_ylim([0, y_max])
    
    # Add D-value lines and bands
    d_legend_elements = add_d_value_lines_and_bands(ax1, stats)
    
    # Create comprehensive legend for top plot - PLACE OUTSIDE
    main_legend = [Line2D([0], [0], marker='o', color='w', markerfacecolor=color, 
                          markersize=8, label='Volume Distribution')]
    
    all_legend_elements = main_legend + fit_legend_elements + d_legend_elements
    leg1 = ax1.legend(handles=all_legend_elements, fontsize=9, frameon=True, 
                     bbox_to_anchor=(1.05, 1), loc='upper left')
    leg1.get_frame().set_alpha(0.95)
    leg1.get_frame().set_edgecolor('lightgray')
    
    ax1.grid(True, linestyle='--', alpha=0.4)
    ax1.set_xlabel('')
    
    # BOTTOM SUBPLOT: CUMULATIVE DISTRIBUTION
    ax2 = fig.add_subplot(gs[1])
    
    if 'volume_nm^3_per_mL_cumsum_avg' in plot_df.columns:
        cumsum_values = plot_df['volume_nm^3_per_mL_cumsum_avg']
        max_cumsum = np.max(cumsum_values)
        if max_cumsum > 0:
            cumsum_percentage = (cumsum_values / max_cumsum) * 100
            ax2.plot(plot_df['size_nm'], cumsum_percentage, '-', 
                    color=color, linewidth=3, alpha=0.9, label='Cumulative %')
            
            if 'volume_nm^3_per_mL_cumsum_sd' in plot_df.columns:
                cumsum_sd = plot_df['volume_nm^3_per_mL_cumsum_sd']
                cumsum_sd_percentage = (cumsum_sd / max_cumsum) * 100
                ax2.fill_between(plot_df['size_nm'], 
                               cumsum_percentage - cumsum_sd_percentage,
                               cumsum_percentage + cumsum_sd_percentage,
                               color=color, alpha=0.25, zorder=1, label='± SD')
        
        ax2.set_ylim([0, 110])
        ax2.set_ylabel('Cumulative Percentage (%)', color=color, fontsize=14, labelpad=10)
        ax2.tick_params(axis='y', labelcolor=color, labelsize=12)
        ax2.spines['left'].set_color(color)
    
    # Format bottom subplot
    ax2.set_xlabel('Size (nm)', fontsize=14, labelpad=10)
    ax2.tick_params(axis='x', labelsize=12)
    ax2.set_xscale(xscale)
    ax2.set_xlim(ax1.get_xlim())
    
    # Add D-value lines to bottom plot
    if 'volume_nm^3_per_mL_cumsum_avg' in plot_df.columns:
        add_d_value_lines_and_bands(ax2, stats)
    
    # Legend for bottom plot
    cumulative_legend = [
        Line2D([0], [0], color=color, linewidth=3, label='Cumulative %'),
        Line2D([0], [0], color=color, alpha=0.25, linewidth=8, label='± SD')
    ]
    
    leg2 = ax2.legend(handles=cumulative_legend, fontsize=10, frameon=True, 
                     bbox_to_anchor=(1.05, 1), loc='upper left')
    leg2.get_frame().set_alpha(0.95)
    leg2.get_frame().set_edgecolor('lightgray')
    
    ax2.grid(True, linestyle='--', alpha=0.4)
    
    # TITLE AND METADATA
    replicate_info = ""
    if metadata and 'num_replicates' in metadata:
        num_reps = metadata['num_replicates']
        if num_reps and str(num_reps) != '1':
            replicate_info = f" (n={num_reps})"
    
    main_title = f'{scale_name} Volume-Weighted\nDistribution: {uniqueID}{replicate_info}'
    fig.suptitle(main_title, fontsize=14, fontweight='bold', y=0.94)
    
    subtitle = f"Error bars/bands: ± SD | Fits: Lognormal"
    fig.text(0.5, 0.87, subtitle, ha='center', fontsize=11, style='italic')
    
    return fig, fit_results


def generate_volume_plots(distribution_df, stats_dict=None, uniqueID=None, 
                         metadata=None, output_dir=None, config=None):
    """Generate volume-weighted distribution plots for both linear and log scales."""
    
    if distribution_df is None or distribution_df.empty:
        return False, "No data available for plotting"
    
    plt.style.use('default')
    
    # Determine output directory
    if output_dir is None:
        if config is not None and "directory" in config:
            base_dir = config["directory"]
            output_dir = os.path.join(base_dir, "processed")
        else:
            output_dir = os.path.join(os.getcwd(), "processed")
    
    try:
        os.makedirs(output_dir, exist_ok=True)
    except Exception as e:
        return False, f"Failed to create output directory: {str(e)}"
    
    created_files = []
    
    # Generate linear and logarithmic plots (each with independent fits)
    for is_log_scale in [False, True]:
        scale_type = 'logarithmic' if is_log_scale else 'linear'
        scale_name = 'log' if is_log_scale else 'linear'
        
        print(f"Creating {scale_name} volume-weighted plot...")
        
        # Filter data for this scale
        plot_df = distribution_df[distribution_df['scale'] == scale_type].copy()
        
        if plot_df.empty:
            print(f"  Warning: No {scale_type} scale data available")
            continue
        
        # Get statistics
        stats = None
        if stats_dict and scale_type in stats_dict and 'volume' in stats_dict[scale_type]:
            stats = stats_dict[scale_type]['volume']
        
        # Create the plot (fits calculated independently for each scale)
        fig, fit_results = create_volume_plot(plot_df, is_log_scale, stats, uniqueID, metadata)
        
        if fig is None:
            print(f"  Failed to create plot")
            continue
        
        # Save fit results to comprehensive fits file
        if fit_results:
            fit_type, fit_data = fit_results
            try:
                
                # Load existing fits file or create new one
                fits_filename = f"Fits_{uniqueID}_all.json"
                fits_path = os.path.join(output_dir, fits_filename)
                
                if os.path.exists(fits_path):
                    with open(fits_path, 'r') as f:
                        all_fits = json.load(f)
                else:
                    all_fits = {'dataset': uniqueID, 'fits': {}}
                
                # Add this fit to the collection
                fit_key = f"volume_{scale_name}"
                all_fits['fits'][fit_key] = {
                    'distribution_type': 'volume',
                    'scale': scale_name,
                    'fit_type': fit_type,
                    'parameters': fit_data
                }
                
                # Save updated fits file
                with open(fits_path, 'w') as f:
                    json.dump(all_fits, f, indent=2, default=str)
                
                print(f"  ✓ Saved fit to: {fits_filename}")
                
            except Exception as e:
                print(f"  ⚠ Failed to save fit: {str(e)}")
                
            except Exception as e:
                print(f"  ⚠ Failed to save fit: {str(e)}")
        
        # Save the plot
        try:
            base_filename = f"Plot_{uniqueID}_{scale_name}_volume"
            
            pdf_path = os.path.join(output_dir, f"{base_filename}.pdf")
            fig.savefig(pdf_path, bbox_inches='tight', dpi=300)
            
            png_path = os.path.join(output_dir, f"{base_filename}.png")
            fig.savefig(png_path, bbox_inches='tight', dpi=300)
            
            created_files.append(pdf_path)
            print(f"  ✓ Saved: {base_filename}.pdf/.png")
            
            plt.close(fig)
            
        except Exception as e:
            print(f"  ✗ Failed to save plot: {str(e)}")
            plt.close(fig)
            continue
    
    return True, created_files


# Execute if we have the required data
if 'current_distribution_df' in globals() and current_distribution_df is not None:
    print("=" * 80)
    print("GENERATING VOLUME-WEIGHTED DISTRIBUTION PLOTS (FIXED VERSION)")
    print("=" * 80)
    
    uniqueID = current_metadata.get('persistentID', 'unknown') if 'current_metadata' in globals() else 'unknown'
    stats = current_stats if 'current_stats' in globals() else None
    metadata = current_metadata if 'current_metadata' in globals() else None
    config = CONFIG if 'CONFIG' in globals() else None
    
    print(f"Creating volume-weighted plots for: {uniqueID}")
    print("Includes: Linear + Logarithmic (Lognormal fits for both)")
    
    success, plot_files = generate_volume_plots(
        current_distribution_df,
        stats_dict=stats,
        uniqueID=uniqueID,
        metadata=metadata,
        config=config
    )
    
    if not success:
        print(f"ERROR: {plot_files}")
    else:
        print(f"\n✓ Successfully created {len(plot_files)} volume-weighted plots!")
        for filepath in plot_files:
            print(f"  - {os.path.basename(filepath)}")
        
        current_volume_plots = plot_files

else:
    print("No data found. Run the complete workflow first.")


# ======================================================================
# CELL 11
# ======================================================================





def lognormal_pdf(x, mu, sigma, amplitude):
    """Calculate lognormal probability density function."""
    return amplitude * (1 / (x * sigma * np.sqrt(2 * np.pi))) * np.exp(-(np.log(x) - mu)**2 / (2 * sigma**2))



def add_surface_area_fit_curve(ax, plot_df, is_log_scale, fit_color='#C45B5B'):
    """Add lognormal fits for surface area distributions. Calculates independent fits per scale."""
    fit_legend_elements = []
    
    sizes = plot_df['size_nm'].values
    surface_areas = plot_df['area_nm^2_per_mL_avg'].values
    
    # Remove any zero or negative values
    valid_mask = (surface_areas > 0) & (sizes > 0)
    if not np.any(valid_mask):
        return fit_legend_elements, None
    
    sizes_valid = sizes[valid_mask]
    surface_areas_valid = surface_areas[valid_mask]
    
    # Calculate fit independently for this scale
    success, result = fit_lognormal_distribution(sizes_valid, surface_areas_valid)
    if not success:
        return fit_legend_elements, None
    
    size_range, fitted_curve, params = result
    
    # Plot the fit curve
    ax.plot(size_range, fitted_curve, '-', color=fit_color, linewidth=2.5, 
           alpha=0.9, label='Lognormal Fit', zorder=4)
    
    # Calculate geometric statistics
    geometric_mean = np.exp(params[0])
    geometric_std = np.exp(params[1])
    
    # Add fit info to legend
    fit_legend_elements.append(
        Line2D([0], [0], color=fit_color, linestyle='-', linewidth=2.5,
              label=f'Lognormal: μ={geometric_mean:.1f} nm, σ={geometric_std:.2f}')
    )
    
    return fit_legend_elements, ('lognormal', {
        'mu': params[0], 
        'sigma': params[1], 
        'amplitude': params[2],
        'geo_mean': geometric_mean, 
        'geo_std': geometric_std
    })


def add_d_value_lines_and_bands(ax, stats):
    """Add D-value lines and uncertainty bands to a subplot."""
    legend_elements = []
    
    if not stats or 'D10_avg' not in stats:
        return legend_elements
    
    d10_avg = stats['D10_avg']
    d10_lower = stats.get('D10_lower', d10_avg)
    d10_upper = stats.get('D10_upper', d10_avg)
    
    d50_avg = stats['D50_avg'] 
    d50_lower = stats.get('D50_lower', d50_avg)
    d50_upper = stats.get('D50_upper', d50_avg)
    
    d90_avg = stats['D90_avg']
    d90_lower = stats.get('D90_lower', d90_avg)
    d90_upper = stats.get('D90_upper', d90_avg)
    
    span = stats.get('span_avg', (d90_avg-d10_avg)/d50_avg if d50_avg > 0 else 0)
    
    # Add D-value lines and bands
    for d_val, d_lower, d_upper, style, width, alpha_band in [
        (d10_avg, d10_lower, d10_upper, '--', 1.5, 0.15),
        (d50_avg, d50_lower, d50_upper, '-', 2.5, 0.25), 
        (d90_avg, d90_lower, d90_upper, '--', 1.5, 0.15)
    ]:
        if not np.isnan(d_val):
            ax.axvline(x=d_val, color='gray', linestyle=style, alpha=0.8, linewidth=width, zorder=5)
            if not np.isnan(d_lower) and not np.isnan(d_upper) and (d_lower != d_val or d_upper != d_val):
                ax.axvspan(d_lower, d_upper, alpha=alpha_band, color='gray', zorder=1)
    
    # Create legend elements
    legend_elements.extend([
        Line2D([0], [0], color='gray', linestyle='--', linewidth=1.5, 
              label=f'D10: {d10_avg:.1f} nm ({d10_lower:.1f}-{d10_upper:.1f})'),
        Line2D([0], [0], color='gray', linestyle='-', linewidth=2.5, 
              label=f'D50: {d50_avg:.1f} nm ({d50_lower:.1f}-{d50_upper:.1f})'),
        Line2D([0], [0], color='gray', linestyle='--', linewidth=1.5, 
              label=f'D90: {d90_avg:.1f} nm ({d90_lower:.1f}-{d90_upper:.1f})'),
        Line2D([0], [0], color='white', linestyle='', 
              label=f'Span: {span:.3f}')
    ])
    
    return legend_elements


def create_surface_area_plot(plot_df, is_log_scale, stats=None, uniqueID=None, metadata=None):
    """Create a two-subplot plot for surface area-weighted distribution."""
    
    scale_name = "Logarithmic" if is_log_scale else "Linear"
    xscale = 'log' if is_log_scale else 'linear'
    color = '#4059AD'  # Indigo blue for surface area-weighted
    
    # Sort by size
    plot_df = plot_df.sort_values('size_nm')
    
    # Create figure
    fig = plt.figure(figsize=(7, 9))
    gs = gridspec.GridSpec(2, 1, height_ratios=[0.6, 0.4], hspace=0.3, 
                          top=0.82, bottom=0.08)
    
    # =================================================================
    # TOP SUBPLOT: MAIN DISTRIBUTION WITH ERROR BARS AND FITS
    # =================================================================
    ax1 = fig.add_subplot(gs[0])
    
    # Plot main distribution with error bars
    if 'area_nm^2_per_mL_sd' in plot_df.columns:
        ax1.errorbar(plot_df['size_nm'], plot_df['area_nm^2_per_mL_avg'], 
                    yerr=plot_df['area_nm^2_per_mL_sd'],
                    fmt='o', color=color, ecolor=color, alpha=0.7,
                    capsize=3, capthick=1, markersize=6, linewidth=1.5,
                    label='Surface Area Distribution')
    else:
        ax1.scatter(plot_df['size_nm'], plot_df['area_nm^2_per_mL_avg'], 
                   color=color, s=60, alpha=0.8, label='Surface Area Distribution')
    
    # Add lognormal fit curve (calculated independently for this scale)
    fit_result = add_surface_area_fit_curve(ax1, plot_df, is_log_scale)
    if isinstance(fit_result, tuple):
        fit_legend_elements, fit_results = fit_result
    else:
        fit_legend_elements = fit_result
        fit_results = None
    
    # Format top subplot
    ax1.set_ylabel('Surface Area (nm²/mL)', color=color, fontsize=14, labelpad=10)
    ax1.tick_params(axis='y', labelcolor=color, labelsize=12)
    ax1.tick_params(axis='x', labelsize=12)
    ax1.spines['left'].set_color(color)
    
    # Set x-axis scale and add better tick labels for log scale
    ax1.set_xscale(xscale)
    if is_log_scale:
        ax1.xaxis.set_major_locator(LogLocator(base=10, numticks=12))
        ax1.xaxis.set_minor_locator(LogLocator(base=10, subs=(0.2, 0.4, 0.6, 0.8), numticks=12))
        ax1.xaxis.set_major_formatter(LogFormatter(base=10, labelOnlyBase=False))
    
    # Smart surface area-weighted range calculation
    weights_for_range = plot_df['area_nm^2_per_mL_avg'].values
    sizes_for_range = plot_df['size_nm'].values
    
    # Find where 99% of the surface area signal is contained
    cumsum_weights = np.cumsum(weights_for_range)
    total_weight = cumsum_weights[-1]
    
    # Find 1st and 99th percentiles of the surface area-weighted distribution
    p1_idx = np.searchsorted(cumsum_weights, 0.01 * total_weight)
    p99_idx = np.searchsorted(cumsum_weights, 0.99 * total_weight)
    
    signal_min = sizes_for_range[max(0, p1_idx)]
    signal_max = sizes_for_range[min(len(sizes_for_range)-1, p99_idx)]
    data_max = plot_df['size_nm'].max()
    
    if is_log_scale:
        # Log scale: focus on the surface area-weighted signal range with some padding
        min_size = max(signal_min * 0.7, 20)  # Don't go below 20 nm
        max_size = min(signal_max * 1.8, data_max * 1.2)  # Surface area less skewed than volume
        ax1.set_xlim([min_size, max_size])
    else:
        # Linear scale: tighten the range more, focus on main signal
        min_size = 0
        max_size = min(signal_max * 1.15, 350)  # Tighter range for surface area
        ax1.set_xlim([min_size, max_size])
    
    # Set y-axis to start from 0
    y_min, y_max = ax1.get_ylim()
    ax1.set_ylim([0, y_max])
    
    # Add D-value lines and bands
    d_legend_elements = add_d_value_lines_and_bands(ax1, stats)
    
    # Create comprehensive legend for top plot - PLACE OUTSIDE
    main_legend = [Line2D([0], [0], marker='o', color='w', markerfacecolor=color, 
                          markersize=8, label='Surface Area Distribution')]
    
    all_legend_elements = main_legend + fit_legend_elements + d_legend_elements
    leg1 = ax1.legend(handles=all_legend_elements, fontsize=9, frameon=True, 
                     bbox_to_anchor=(1.05, 1), loc='upper left')
    leg1.get_frame().set_alpha(0.95)
    leg1.get_frame().set_edgecolor('lightgray')
    
    ax1.grid(True, linestyle='--', alpha=0.4)
    ax1.set_xlabel('')  # No x-label on top plot
    
    # =================================================================
    # BOTTOM SUBPLOT: CUMULATIVE DISTRIBUTION WITH UNCERTAINTY BANDS
    # =================================================================
    ax2 = fig.add_subplot(gs[1])
    
    if 'area_nm^2_per_mL_cumsum_avg' in plot_df.columns:
        # For surface area cumulative, normalize to percentage
        cumsum_values = plot_df['area_nm^2_per_mL_cumsum_avg']
        max_cumsum = np.max(cumsum_values)
        if max_cumsum > 0:
            cumsum_percentage = (cumsum_values / max_cumsum) * 100
            ax2.plot(plot_df['size_nm'], cumsum_percentage, '-', 
                    color=color, linewidth=3, alpha=0.9, label='Cumulative %')
            
            # Add uncertainty bands if available
            if 'area_nm^2_per_mL_cumsum_sd' in plot_df.columns:
                cumsum_sd = plot_df['area_nm^2_per_mL_cumsum_sd']
                cumsum_sd_percentage = (cumsum_sd / max_cumsum) * 100
                ax2.fill_between(plot_df['size_nm'], 
                               cumsum_percentage - cumsum_sd_percentage,
                               cumsum_percentage + cumsum_sd_percentage,
                               color=color, alpha=0.25, zorder=1, label='± SD')
        
        ax2.set_ylim([0, 110])
        ax2.set_ylabel('Cumulative Percentage (%)', color=color, fontsize=14, labelpad=10)
        ax2.tick_params(axis='y', labelcolor=color, labelsize=12)
        ax2.spines['left'].set_color(color)
    
    # Format bottom subplot
    ax2.set_xlabel('Size (nm)', fontsize=14, labelpad=10)
    ax2.tick_params(axis='x', labelsize=12)
    ax2.set_xscale(xscale)
    ax2.set_xlim(ax1.get_xlim())  # Match top plot limits
    
    if is_log_scale:
        # Add same detailed log scale ticks to bottom plot
        ax2.xaxis.set_major_locator(LogLocator(base=10, numticks=12))
        ax2.xaxis.set_minor_locator(LogLocator(base=10, subs=(0.2, 0.4, 0.6, 0.8), numticks=12))
        ax2.xaxis.set_major_formatter(LogFormatter(base=10, labelOnlyBase=False))
    
    # Add D-value lines to bottom plot
    if 'area_nm^2_per_mL_cumsum_avg' in plot_df.columns:
        add_d_value_lines_and_bands(ax2, stats)
    
    # Legend for bottom plot - PLACE OUTSIDE
    cumulative_legend = [
        Line2D([0], [0], color=color, linewidth=3, label='Cumulative %'),
        Line2D([0], [0], color=color, alpha=0.25, linewidth=8, label='± SD')
    ]
    
    leg2 = ax2.legend(handles=cumulative_legend, fontsize=10, frameon=True, 
                     bbox_to_anchor=(1.05, 1), loc='upper left')
    leg2.get_frame().set_alpha(0.95)
    leg2.get_frame().set_edgecolor('lightgray')
    
    ax2.grid(True, linestyle='--', alpha=0.4)
    
    # =================================================================
    # TITLE AND METADATA
    # =================================================================
    
    # Extract replicate info
    replicate_info = ""
    if metadata and 'num_replicates' in metadata:
        num_reps = metadata['num_replicates']
        if num_reps and str(num_reps) != '1':
            replicate_info = f" (n={num_reps})"
    
    # Set main title
    main_title = f'{scale_name} Surface Area-Weighted\nDistribution: {uniqueID}{replicate_info}'
    fig.suptitle(main_title, fontsize=14, fontweight='bold', y=0.94)
    
    # Add subtitle
    subtitle = f"Error bars/bands: ± SD | Fits: Lognormal"
    fig.text(0.5, 0.87, subtitle, ha='center', fontsize=11, style='italic')
    
    return fig, fit_results


def generate_surface_area_plots(distribution_df, stats_dict=None, uniqueID=None, 
                               metadata=None, output_dir=None, config=None):
    """Generate surface area-weighted distribution plots for both linear and log scales."""
    
    if distribution_df is None or distribution_df.empty:
        return False, "No data available for plotting"
    
    plt.style.use('default')
    
    # Determine output directory
    if output_dir is None:
        if config is not None and "directory" in config:
            base_dir = config["directory"]
            output_dir = os.path.join(base_dir, "processed")
        else:
            output_dir = os.path.join(os.getcwd(), "processed")
    
    try:
        os.makedirs(output_dir, exist_ok=True)
    except Exception as e:
        return False, f"Failed to create output directory: {str(e)}"
    
    created_files = []
    
    # Generate linear and logarithmic plots (each with independent fits)
    for is_log_scale in [False, True]:
        scale_type = 'logarithmic' if is_log_scale else 'linear'
        scale_name = 'log' if is_log_scale else 'linear'
        
        print(f"Creating {scale_name} surface area-weighted plot...")
        
        # Filter data for this scale
        plot_df = distribution_df[distribution_df['scale'] == scale_type].copy()
        
        if plot_df.empty:
            print(f"  Warning: No {scale_type} scale data available")
            continue
        
        # Get statistics
        stats = None
        if stats_dict and scale_type in stats_dict and 'surface_area' in stats_dict[scale_type]:
            stats = stats_dict[scale_type]['surface_area']
        
        # Create the plot (fits calculated independently for each scale)
        fig, fit_results = create_surface_area_plot(plot_df, is_log_scale, stats, uniqueID, metadata)
        
        if fig is None:
            print(f"  Failed to create plot")
            continue
        
        # Save fit results to comprehensive fits file
        if fit_results:
            fit_type, fit_data = fit_results
            try:
                
                # Load existing fits file or create new one
                fits_filename = f"Fits_{uniqueID}_all.json"
                fits_path = os.path.join(output_dir, fits_filename)
                
                if os.path.exists(fits_path):
                    with open(fits_path, 'r') as f:
                        all_fits = json.load(f)
                else:
                    all_fits = {'dataset': uniqueID, 'fits': {}}
                
                # Add this fit to the collection
                fit_key = f"surface_area_{scale_name}"
                all_fits['fits'][fit_key] = {
                    'distribution_type': 'surface_area',
                    'scale': scale_name,
                    'fit_type': fit_type,
                    'parameters': fit_data
                }
                
                # Save updated fits file
                with open(fits_path, 'w') as f:
                    json.dump(all_fits, f, indent=2, default=str)
                
                print(f"  ✓ Saved fit to: {fits_filename}")
                
            except Exception as e:
                print(f"  ⚠ Failed to save fit: {str(e)}")
        
        # Save the plot
        try:
            base_filename = f"Plot_{uniqueID}_{scale_name}_surface_area"
            
            pdf_path = os.path.join(output_dir, f"{base_filename}.pdf")
            fig.savefig(pdf_path, bbox_inches='tight', dpi=300)
            
            png_path = os.path.join(output_dir, f"{base_filename}.png")
            fig.savefig(png_path, bbox_inches='tight', dpi=300)
            
            created_files.append(pdf_path)
            print(f"  ✓ Saved: {base_filename}.pdf/.png")
            
            plt.close(fig)
            
        except Exception as e:
            print(f"  ✗ Failed to save plot: {str(e)}")
            plt.close(fig)
            continue
    
    return True, created_files


# Execute if we have the required data
if 'current_distribution_df' in globals() and current_distribution_df is not None:
    print("=" * 80)
    print("GENERATING SURFACE AREA-WEIGHTED DISTRIBUTION PLOTS")
    print("=" * 80)
    
    uniqueID = current_metadata.get('persistentID', 'unknown') if 'current_metadata' in globals() else 'unknown'
    stats = current_stats if 'current_stats' in globals() else None
    metadata = current_metadata if 'current_metadata' in globals() else None
    config = CONFIG if 'CONFIG' in globals() else None
    
    print(f"Creating surface area-weighted plots for: {uniqueID}")
    print("Includes: Linear + Logarithmic (Lognormal fits for both)")
    
    success, plot_files = generate_surface_area_plots(
        current_distribution_df,
        stats_dict=stats,
        uniqueID=uniqueID,
        metadata=metadata,
        config=config
    )
    
    if not success:
        print(f"ERROR: {plot_files}")
    else:
        print(f"\n✓ Successfully created {len(plot_files)} surface area-weighted plots!")
        for filepath in plot_files:
            print(f"  - {os.path.basename(filepath)}")
        
        current_surface_area_plots = plot_files

else:
    print("No data found. Run the complete workflow first.")


# ======================================================================
# CELL 12
# ======================================================================





def plot_raw_counts_with_settings(distribution_df, uniqueID=None, metadata=None, config=None):
    """
    Create visualizations of raw particle counts vs. size for averaged data, with instrument settings.
    
    Parameters:
    distribution_df (DataFrame): Processed NTA data with averaged values
    uniqueID (str): Unique identifier for the dataset
    metadata (dict): Metadata dictionary
    config (dict): Configuration dictionary
    
    Returns:
    tuple: (success_flag, filepath)
    """
    # Validate input
    if distribution_df is None or distribution_df.empty:
        print("No data available for plotting")
        return False, "No data available for plotting"
    
    # Extract sample ID
    id_text = uniqueID or (distribution_df['uniqueID'].iloc[0] if 'uniqueID' in distribution_df.columns else 'Unknown')
    
    # Create figure
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    
    # Define color scheme - same as other plots
    color = '#4C5B5C'  # Slate gray for raw counts
    
    # Determine which columns to use for raw counts (before dilution correction)
    raw_count_col = None
    raw_count_sd_col = None
    
    # Look for the averaged number data (before normalization)
    if 'number_avg' in distribution_df.columns:
        raw_count_col = 'number_avg'
        raw_count_sd_col = 'number_sd' if 'number_sd' in distribution_df.columns else None
    else:
        print("Raw count data column not found in the dataset")
        return False, "Raw count data column not found in the dataset"
    
    # Plot 1: Linear scale raw counts
    lin_df = distribution_df[distribution_df['scale'] == 'linear'].sort_values('size_nm')
    if not lin_df.empty:
        # Plot raw counts with error bars if available
        if raw_count_sd_col and raw_count_sd_col in lin_df.columns:
            ax1.errorbar(lin_df['size_nm'], lin_df[raw_count_col], 
                        yerr=lin_df[raw_count_sd_col],
                        fmt='o', color=color, ecolor=color, alpha=0.7,
                        capsize=3, capthick=1, markersize=4, linewidth=1.5)
        else:
            ax1.scatter(lin_df['size_nm'], lin_df[raw_count_col], color=color, s=30, alpha=0.8)
        
        ax1.set_title('Raw Particle Counts - Linear Scale', fontsize=14)
        ax1.set_xlabel('Size (nm)', fontsize=12)
        ax1.set_ylabel('Raw Counts (particles)', fontsize=12, color=color)
        ax1.tick_params(axis='y', labelcolor=color)
        ax1.spines['left'].set_color(color)
        ax1.grid(True, alpha=0.3)
        
        # Set reasonable x-axis limits
        percentile_95 = np.percentile(lin_df['size_nm'], 95)
        max_size = min(percentile_95 * 1.2, 500)  # Cap at 500 nm or use 95th percentile + 20%
        ax1.set_xlim([0, max_size])
        
        # Set y-axis to start from 0
        y_min, y_max = ax1.get_ylim()
        ax1.set_ylim([0, y_max])
    
    # Plot 2: Log scale raw counts 
    log_df = distribution_df[distribution_df['scale'] == 'logarithmic'].sort_values('size_nm')
    if not log_df.empty:
        # Plot raw counts with error bars if available
        if raw_count_sd_col and raw_count_sd_col in log_df.columns:
            ax2.errorbar(log_df['size_nm'], log_df[raw_count_col], 
                        yerr=log_df[raw_count_sd_col],
                        fmt='o', color=color, ecolor=color, alpha=0.7,
                        capsize=3, capthick=1, markersize=4, linewidth=1.5)
        else:
            ax2.scatter(log_df['size_nm'], log_df[raw_count_col], color=color, s=30, alpha=0.8)
        
        ax2.set_title('Raw Particle Counts - Log Scale', fontsize=14)
        ax2.set_xlabel('Size (nm)', fontsize=12)
        ax2.set_ylabel('Raw Counts (particles)', fontsize=12, color=color)
        ax2.tick_params(axis='y', labelcolor=color)
        ax2.spines['left'].set_color(color)
        ax2.set_xscale('log')
        ax2.grid(True, alpha=0.3)
        
        # Set reasonable x-axis limits for log scale
        min_size = max(30, log_df['size_nm'].min())  # Don't go below 30 nm
        percentile_90 = np.percentile(log_df['size_nm'], 90)
        max_size = min(percentile_90 * 1.3, 300)
        ax2.set_xlim([min_size, max_size])
        
        # Set y-axis to start from 0
        y_min, y_max = ax2.get_ylim()
        ax2.set_ylim([0, y_max])
    
    # Determine replicate info for title
    replicate_info = ""
    if metadata and 'num_replicates' in metadata:
        num_reps = metadata['num_replicates']
        if num_reps and str(num_reps) != '1':
            replicate_info = f" (n={num_reps} replicates)"
    
    # Add overall title
    fig.suptitle(f'Raw Particle Counts vs Size: {id_text}{replicate_info}', fontsize=16, y=0.98)
    
    # Calculate total raw count (use linear data preferentially)
    if not lin_df.empty:
        total_raw_count = lin_df[raw_count_col].sum()
        if raw_count_sd_col and raw_count_sd_col in lin_df.columns:
            # Error propagation for sum: sqrt(sum of variances)
            total_raw_count_sd = np.sqrt((lin_df[raw_count_sd_col] ** 2).sum())
        else:
            total_raw_count_sd = 0
    elif not log_df.empty:
        total_raw_count = log_df[raw_count_col].sum()
        if raw_count_sd_col and raw_count_sd_col in log_df.columns:
            total_raw_count_sd = np.sqrt((log_df[raw_count_sd_col] ** 2).sum())
        else:
            total_raw_count_sd = 0
    else:
        total_raw_count = 0
        total_raw_count_sd = 0
    
    # Extract and format instrument settings from metadata
    instrument_settings = []
    if metadata:
        # Extract key NTA parameters
        cycles = metadata.get('nta_cycles', metadata.get('cycles', 'Unknown'))
        if cycles != 'Unknown':
            instrument_settings.append(f"Cycles: {cycles}")
        
        fps = metadata.get('nta_fps', metadata.get('fps', 'Unknown'))
        if fps != 'Unknown':
            instrument_settings.append(f"FPS: {fps}")
        
        # Extract total number of traces (summed across replicates)
        traces = metadata.get('nta_number_of_traces_sum', metadata.get('nta_number_of_traces', 'Unknown'))
        if traces != 'Unknown':
            instrument_settings.append(f"Total traces: {traces}")
        
        # Extract video file info
        avi_size = metadata.get('nta_avi_filesize', 'Unknown')
        if avi_size != 'Unknown':
            instrument_settings.append(f"Video size: {avi_size}")
        
        # Extract temperature
        temp = metadata.get('nta_temperature', metadata.get('temperature', 'Unknown'))
        if temp != 'Unknown':
            instrument_settings.append(f"Temp: {temp}")
    
    # Add explanation text about what raw counts represent
    explanation_parts = [
        "Raw counts represent the actual number of particles detected at each size bin,",
        "before normalization or concentration calculation"
    ]
    
    if replicate_info:
        explanation_parts.append("Error bars show ± SD across replicates")
    
    explanation = " ".join(explanation_parts)
    fig.text(0.5, 0.01, explanation, ha='center', fontsize=10, style='italic')
    
    # Add data statistics and instrument settings
    if total_raw_count_sd > 0:
        stats_text = f"Total raw particles detected: {total_raw_count:.0f} ± {total_raw_count_sd:.0f}"
    else:
        stats_text = f"Total raw particles detected: {total_raw_count:.0f}"
    
    # Create text block for statistics and settings
    text_block = [stats_text]
    if instrument_settings:
        text_block.append("Instrument Settings: " + " | ".join(instrument_settings))
    
    # Add the text block to the figure
    fig.text(0.5, 0.04, "\n".join(text_block), ha='center', fontsize=11)
    
    # Adjust layout to make room for the text
    bottom_margin = 0.06 + 0.02 * len(text_block)
    plt.tight_layout(rect=[0, bottom_margin, 1, 0.95])
    
    # Save figure
    filepath = None
    
    # Determine the output directory
    if config is not None and "directory" in config:
        output_dir = os.path.join(config["directory"], "processed")
    else:
        output_dir = os.path.join(os.getcwd(), "processed")
    
    # Ensure the output directory exists
    try:
        os.makedirs(output_dir, exist_ok=True)
        
        # Create descriptive filename
        base_path = os.path.join(output_dir, f"Plot_{id_text}_raw_counts")
        
        # Save in both PDF and PNG formats
        pdf_path = f"{base_path}.pdf"
        plt.savefig(pdf_path, bbox_inches='tight', dpi=300)
        
        png_path = f"{base_path}.png"
        plt.savefig(png_path, bbox_inches='tight', dpi=300)
        
        print(f"✓ Saved raw count plots: {os.path.basename(pdf_path)}/.png")
        filepath = pdf_path
    except Exception as e:
        print(f"⚠ Failed to save plot: {str(e)}")
        return False, str(e)
    
    # Update metadata with total raw particle count
    if metadata and total_raw_count > 0:
        try:
            # Get output directory for metadata
            if config is not None and "directory" in config:
                metadata_dir = os.path.join(config["directory"], "metadata")
            else:
                metadata_dir = os.path.join(os.getcwd(), "metadata")
            
            # Ensure metadata directory exists
            os.makedirs(metadata_dir, exist_ok=True)
            
            # Construct metadata filepath
            unique_id = metadata.get('persistentID', id_text)
            metadata_path = os.path.join(metadata_dir, f"Data_{unique_id}_metadata.txt")
            
            # Read existing metadata to preserve all fields
            existing_data = {}
            if os.path.exists(metadata_path):
                with open(metadata_path, 'r') as f:
                    for line in f:
                        parts = line.strip().split('\t')
                        if len(parts) >= 2:
                            existing_data[parts[0]] = parts[1]
            
            # Add total raw particle count with uncertainty
            if total_raw_count_sd > 0:
                existing_data['nta_total_raw_particles'] = f"{total_raw_count:.0f} ± {total_raw_count_sd:.0f}"
            else:
                existing_data['nta_total_raw_particles'] = f"{total_raw_count:.0f}"
            
            # Write updated metadata
            with open(metadata_path, 'w') as f:
                for key, value in existing_data.items():
                    f.write(f"{key}\t{value}\t\n")
            
            print(f"✓ Updated metadata with total raw particle count")
            
            # Update global metadata variable if it exists
            if 'current_metadata' in globals():
                if total_raw_count_sd > 0:
                    globals()['current_metadata']['nta_total_raw_particles'] = f"{total_raw_count:.0f} ± {total_raw_count_sd:.0f}"
                else:
                    globals()['current_metadata']['nta_total_raw_particles'] = f"{total_raw_count:.0f}"
                
        except Exception as e:
            print(f"⚠ Could not update metadata: {str(e)}")
    
    plt.show()
    return True, filepath


# Execute if we have the required data
if 'current_distribution_df' in globals() and current_distribution_df is not None:
    print("=" * 80)
    print("GENERATING RAW PARTICLE COUNT VISUALIZATION")
    print("=" * 80)
    
    uniqueID = current_metadata.get('persistentID', 'unknown') if 'current_metadata' in globals() else 'unknown'
    metadata = current_metadata if 'current_metadata' in globals() else None
    config = CONFIG if 'CONFIG' in globals() else None
    
    print(f"Creating raw particle count plots for: {uniqueID}")
    
    success, filepath = plot_raw_counts_with_settings(
        current_distribution_df, 
        uniqueID, 
        metadata,
        config
    )
    
    if success:
        current_raw_count_plot = filepath
        print(f"\n✓ Successfully created raw count visualization!")
    else:
        print(f"ERROR: {filepath}")

else:
    print("No data found. Run the complete workflow first.")


# ======================================================================
# CELL 13
# ======================================================================





def lognormal_pdf(x, mu, sigma, amplitude):
    """Calculate lognormal probability density function with numerical stability."""
    # Add small epsilon to avoid log(0)
    x_safe = np.maximum(x, 1e-10)
    
    # Calculate the lognormal PDF
    log_x = np.log(x_safe)
    exponent = -((log_x - mu) ** 2) / (2 * sigma ** 2)
    
    # Avoid numerical overflow/underflow
    exponent = np.clip(exponent, -50, 50)
    
    pdf = amplitude * np.exp(exponent) / (x_safe * sigma * np.sqrt(2 * np.pi))
    
    # Ensure no NaN or inf values
    pdf = np.nan_to_num(pdf, nan=0.0, posinf=0.0, neginf=0.0)
    
    return pdf



def add_d_value_lines_and_bands_surface_area(ax, stats):
    """Add D-value lines and uncertainty bands to a subplot, converted to surface area in µm²."""
    legend_elements = []
    
    if not stats or 'D10_avg' not in stats:
        return legend_elements
    
    # Convert diameter D-values to surface area D-values using A = π × d² and convert to µm²
    d10_avg_size = stats['D10_avg']
    d10_lower_size = stats.get('D10_lower', d10_avg_size)
    d10_upper_size = stats.get('D10_upper', d10_avg_size)
    
    d50_avg_size = stats['D50_avg'] 
    d50_lower_size = stats.get('D50_lower', d50_avg_size)
    d50_upper_size = stats.get('D50_upper', d50_avg_size)
    
    d90_avg_size = stats['D90_avg']
    d90_lower_size = stats.get('D90_lower', d90_avg_size)
    d90_upper_size = stats.get('D90_upper', d90_avg_size)
    
    # Convert to surface areas in µm² (divide by 1e6 to convert from nm² to µm²)
    d10_avg_sa = (np.pi * (d10_avg_size ** 2)) / 1e6
    d10_lower_sa = (np.pi * (d10_lower_size ** 2)) / 1e6
    d10_upper_sa = (np.pi * (d10_upper_size ** 2)) / 1e6
    
    d50_avg_sa = (np.pi * (d50_avg_size ** 2)) / 1e6
    d50_lower_sa = (np.pi * (d50_lower_size ** 2)) / 1e6
    d50_upper_sa = (np.pi * (d50_upper_size ** 2)) / 1e6
    
    d90_avg_sa = (np.pi * (d90_avg_size ** 2)) / 1e6
    d90_lower_sa = (np.pi * (d90_lower_size ** 2)) / 1e6
    d90_upper_sa = (np.pi * (d90_upper_size ** 2)) / 1e6
    
    span = stats.get('span_avg', (d90_avg_size-d10_avg_size)/d50_avg_size if d50_avg_size > 0 else 0)
    
    # Add D-value lines and bands using surface area values in µm²
    for d_val_sa, d_lower_sa, d_upper_sa, style, width, alpha_band in [
        (d10_avg_sa, d10_lower_sa, d10_upper_sa, '--', 1.5, 0.15),
        (d50_avg_sa, d50_lower_sa, d50_upper_sa, '-', 2.5, 0.25), 
        (d90_avg_sa, d90_lower_sa, d90_upper_sa, '--', 1.5, 0.15)
    ]:
        if not np.isnan(d_val_sa):
            ax.axvline(x=d_val_sa, color='gray', linestyle=style, alpha=0.8, linewidth=width, zorder=5)
            if not np.isnan(d_lower_sa) and not np.isnan(d_upper_sa) and (d_lower_sa != d_val_sa or d_upper_sa != d_val_sa):
                ax.axvspan(d_lower_sa, d_upper_sa, alpha=alpha_band, color='gray', zorder=1)
    
    # Create legend elements showing diameter in nm and surface area in µm²
    legend_elements.extend([
        Line2D([0], [0], color='gray', linestyle='--', linewidth=1.5, 
              label=f'D10: {d10_avg_size:.1f} nm → {d10_avg_sa:.3f} µm² ({d10_lower_size:.1f}-{d10_upper_size:.1f} nm)'),
        Line2D([0], [0], color='gray', linestyle='-', linewidth=2.5, 
              label=f'D50: {d50_avg_size:.1f} nm → {d50_avg_sa:.3f} µm² ({d50_lower_size:.1f}-{d50_upper_size:.1f} nm)'),
        Line2D([0], [0], color='gray', linestyle='--', linewidth=1.5, 
              label=f'D90: {d90_avg_size:.1f} nm → {d90_avg_sa:.3f} µm² ({d90_lower_size:.1f}-{d90_upper_size:.1f} nm)'),
        Line2D([0], [0], color='white', linestyle='', 
              label=f'Span: {span:.3f}')
    ])
    
    return legend_elements


def add_count_fit_curve(ax, plot_df, is_log_scale, fit_color='#F25C54'):
    """Add lognormal fits for count vs surface area distributions in µm²."""
    fit_legend_elements = []
    
    surface_areas = plot_df['surface_area_um2'].values
    counts = plot_df['number_avg'].values
    
    # Remove any zero or negative values
    valid_mask = (counts > 0) & (surface_areas > 0) & np.isfinite(counts) & np.isfinite(surface_areas)
    if not np.any(valid_mask):
        print(f"      No valid data points for fitting")
        return fit_legend_elements, None
    
    surface_areas_clean = surface_areas[valid_mask]
    counts_clean = counts[valid_mask]
    
    print(f"      Fitting with {len(surface_areas_clean)} valid data points")
    print(f"      Surface area range: {surface_areas_clean.min():.4f} - {surface_areas_clean.max():.4f} µm²")
    print(f"      Count range: {counts_clean.min():.1f} - {counts_clean.max():.1f}")
    
    # Use lognormal fit for both linear and log scales
    success, result = fit_lognormal_distribution(surface_areas_clean, counts_clean)
    if success:
        sa_range, fitted_curve, params = result
        ax.plot(sa_range, fitted_curve, '-', color=fit_color, linewidth=2.5, 
               alpha=0.9, label='Lognormal Fit', zorder=4)
        
        geometric_mean_sa = np.exp(params[0])
        geometric_std = np.exp(params[1])
        
        print(f"      Fit parameters: geo_mean={geometric_mean_sa:.4f} µm², geo_std={geometric_std:.3f}")
        
        # Add fit info to legend with µm² units
        fit_legend_elements.append(
            Line2D([0], [0], color=fit_color, linestyle='-', linewidth=2.5,
                  label=f'Lognormal: geo_mean={geometric_mean_sa:.3f} µm², geo_std={geometric_std:.2f}')
        )
        
        return fit_legend_elements, ('lognormal', {'mu': params[0], 'sigma': params[1], 'amplitude': params[2]})
    else:
        print(f"      Lognormal fit failed: {result}")
        return fit_legend_elements, None


def create_count_vs_surface_area_plot(plot_df, is_log_scale, stats=None, uniqueID=None, metadata=None):
    """Create a two-subplot plot for raw counts vs theoretical surface area."""
    
    scale_name = "Logarithmic" if is_log_scale else "Linear"
    xscale = 'log' if is_log_scale else 'linear'
    color = '#4059AD'  # Indigo blue (same as cell_10c surface area plots)
    
    # Calculate theoretical surface area for spherical particles: A = π * d² and convert to µm²
    plot_df = plot_df.copy()
    plot_df['surface_area_um2'] = (np.pi * (plot_df['size_nm'] ** 2)) / 1e6  # Convert nm² to µm²
    
    # Sort by surface area
    plot_df = plot_df.sort_values('surface_area_um2')
    
    # Create figure
    fig = plt.figure(figsize=(7, 9))
    gs = gridspec.GridSpec(2, 1, height_ratios=[0.6, 0.4], hspace=0.3, 
                          top=0.82, bottom=0.08)
    
    # =================================================================
    # TOP SUBPLOT: MAIN DISTRIBUTION WITH ERROR BARS AND FITS
    # =================================================================
    ax1 = fig.add_subplot(gs[0])
    
    # Plot main distribution with error bars
    if 'number_sd' in plot_df.columns:
        ax1.errorbar(plot_df['surface_area_um2'], plot_df['number_avg'], 
                    yerr=plot_df['number_sd'],
                    fmt='o', color=color, ecolor=color, alpha=0.7,
                    capsize=3, capthick=1, markersize=6, linewidth=1.5,
                    label='Raw Counts')
    else:
        ax1.scatter(plot_df['surface_area_um2'], plot_df['number_avg'], 
                   color=color, s=60, alpha=0.8, label='Raw Counts')
    
    # Add lognormal fit curve
    fit_result = add_count_fit_curve(ax1, plot_df, is_log_scale)
    if isinstance(fit_result, tuple):
        fit_legend_elements, fit_results = fit_result
    else:
        fit_legend_elements = fit_result
        fit_results = None
    
    # Debug: Print fit status
    if fit_results:
        print(f"    ✓ Lognormal fit successful for {scale_name} scale")
    else:
        print(f"    ⚠ Lognormal fit failed for {scale_name} scale")
    
    # Format top subplot
    ax1.set_ylabel('Raw Particle Counts', color=color, fontsize=14, labelpad=10)
    ax1.tick_params(axis='y', labelcolor=color, labelsize=12)
    ax1.tick_params(axis='x', labelsize=12)
    ax1.spines['left'].set_color(color)
    
    # Set x-axis scale and add better tick labels for log scale
    ax1.set_xscale(xscale)
    if is_log_scale:
        ax1.xaxis.set_major_locator(LogLocator(base=10, numticks=12))
        ax1.xaxis.set_minor_locator(LogLocator(base=10, subs=(0.2, 0.4, 0.6, 0.8), numticks=12))
        ax1.xaxis.set_major_formatter(LogFormatter(base=10, labelOnlyBase=False))
    
    # Smart range calculation based on where the count signal is
    weights_for_range = plot_df['number_avg'].values
    sa_for_range = plot_df['surface_area_um2'].values
    
    # Find where 99% of the count signal is contained
    cumsum_weights = np.cumsum(weights_for_range)
    total_weight = cumsum_weights[-1]
    
    if total_weight > 0:
        # Find 1st and 99th percentiles of the count-weighted distribution
        p1_idx = np.searchsorted(cumsum_weights, 0.01 * total_weight)
        p99_idx = np.searchsorted(cumsum_weights, 0.99 * total_weight)
        
        signal_min_sa = sa_for_range[max(0, p1_idx)]
        signal_max_sa = sa_for_range[min(len(sa_for_range)-1, p99_idx)]
        data_max_sa = plot_df['surface_area_um2'].max()
        
        if is_log_scale:
            # Log scale: focus on the count signal range with some padding
            min_sa = max(signal_min_sa * 0.7, (np.pi * (20**2)) / 1e6)  # Don't go below 20nm diameter equivalent
            max_sa = min(signal_max_sa * 1.8, data_max_sa * 1.2)
            ax1.set_xlim([min_sa, max_sa])
        else:
            # Linear scale: tighten the range, focus on main signal
            min_sa = 0
            max_sa = min(signal_max_sa * 1.15, (np.pi * (350**2)) / 1e6)  # Cap at ~350nm diameter equivalent
            ax1.set_xlim([min_sa, max_sa])
    else:
        # Fallback if no signal
        if is_log_scale:
            ax1.set_xlim([(np.pi * (30**2)) / 1e6, (np.pi * (300**2)) / 1e6])
        else:
            ax1.set_xlim([0, (np.pi * (250**2)) / 1e6])
    
    # Set y-axis to start from 0
    y_min, y_max = ax1.get_ylim()
    ax1.set_ylim([0, y_max])
    
    # Add D-value lines and bands (converted to surface area)
    d_legend_elements = add_d_value_lines_and_bands_surface_area(ax1, stats)
    
    # Create comprehensive legend for top plot - PLACE OUTSIDE
    main_legend = [Line2D([0], [0], marker='o', color='w', markerfacecolor=color, 
                          markersize=8, label='Raw Counts')]
    
    all_legend_elements = main_legend + fit_legend_elements + d_legend_elements
    leg1 = ax1.legend(handles=all_legend_elements, fontsize=9, frameon=True, 
                     bbox_to_anchor=(1.05, 1), loc='upper left')
    leg1.get_frame().set_alpha(0.95)
    leg1.get_frame().set_edgecolor('lightgray')
    
    ax1.grid(True, linestyle='--', alpha=0.4)
    ax1.set_xlabel('')  # No x-label on top plot
    
    # =================================================================
    # BOTTOM SUBPLOT: CUMULATIVE DISTRIBUTION WITH UNCERTAINTY BANDS
    # =================================================================
    ax2 = fig.add_subplot(gs[1])
    
    # Calculate cumulative counts as percentage
    cumsum_counts = np.cumsum(plot_df['number_avg'])
    max_cumsum = cumsum_counts.iloc[-1] if len(cumsum_counts) > 0 else 1
    
    if max_cumsum > 0:
        cumsum_percentage = (cumsum_counts / max_cumsum) * 100
        ax2.plot(plot_df['surface_area_um2'], cumsum_percentage, '-', 
                color=color, linewidth=3, alpha=0.9, label='Cumulative %')
        
        # Add uncertainty bands if available
        if 'number_sd' in plot_df.columns:
            # Calculate cumulative uncertainty
            cumsum_sd = np.sqrt(np.cumsum(plot_df['number_sd'] ** 2))
            cumsum_sd_percentage = (cumsum_sd / max_cumsum) * 100
            ax2.fill_between(plot_df['surface_area_um2'], 
                           cumsum_percentage - cumsum_sd_percentage,
                           cumsum_percentage + cumsum_sd_percentage,
                           color=color, alpha=0.25, zorder=1, label='± SD')
    
    ax2.set_ylim([0, 110])
    ax2.set_ylabel('Cumulative Percentage (%)', color=color, fontsize=14, labelpad=10)
    ax2.tick_params(axis='y', labelcolor=color, labelsize=12)
    ax2.spines['left'].set_color(color)
    
    # Format bottom subplot
    ax2.set_xlabel('Theoretical Surface Area (µm²)', fontsize=14, labelpad=10)
    ax2.tick_params(axis='x', labelsize=12)
    ax2.set_xscale(xscale)
    ax2.set_xlim(ax1.get_xlim())  # Match top plot limits
    
    if is_log_scale:
        # Add same detailed log scale ticks to bottom plot
        ax2.xaxis.set_major_locator(LogLocator(base=10, numticks=12))
        ax2.xaxis.set_minor_locator(LogLocator(base=10, subs=(0.2, 0.4, 0.6, 0.8), numticks=12))
        ax2.xaxis.set_major_formatter(LogFormatter(base=10, labelOnlyBase=False))
    
    # Add D-value lines to bottom plot (converted to surface area)
    add_d_value_lines_and_bands_surface_area(ax2, stats)
    
    # Legend for bottom plot - PLACE OUTSIDE
    cumulative_legend = [
        Line2D([0], [0], color=color, linewidth=3, label='Cumulative %'),
        Line2D([0], [0], color=color, alpha=0.25, linewidth=8, label='± SD')
    ]
    
    leg2 = ax2.legend(handles=cumulative_legend, fontsize=10, frameon=True, 
                     bbox_to_anchor=(1.05, 1), loc='upper left')
    leg2.get_frame().set_alpha(0.95)
    leg2.get_frame().set_edgecolor('lightgray')
    
    ax2.grid(True, linestyle='--', alpha=0.4)
    
    # =================================================================
    # TITLE AND METADATA
    # =================================================================
    
    # Extract replicate info
    replicate_info = ""
    if metadata and 'num_replicates' in metadata:
        num_reps = metadata['num_replicates']
        if num_reps and str(num_reps) != '1':
            replicate_info = f" (n={num_reps})"
    
    # Set main title
    main_title = f'{scale_name} Raw Counts vs\nTheoretical Surface Area: {uniqueID}{replicate_info}'
    fig.suptitle(main_title, fontsize=14, fontweight='bold', y=0.94)
    
    # Add subtitle
    subtitle = f"Error bars/bands: ± SD | Surface Area = π × diameter² | Fits: Lognormal"
    fig.text(0.5, 0.87, subtitle, ha='center', fontsize=11, style='italic')
    
    return fig, fit_results


def generate_count_vs_surface_area_plots(distribution_df, stats_dict=None, uniqueID=None, 
                                        metadata=None, output_dir=None, config=None):
    """Generate raw count vs surface area plots for both linear and log scales."""
    
    if distribution_df is None or distribution_df.empty:
        return False, "No data available for plotting"
    
    # Check if we have the required columns
    if 'number_avg' not in distribution_df.columns:
        return False, "Raw count data (number_avg) not found in dataset"
    
    plt.style.use('default')
    
    # Determine output directory
    if output_dir is None:
        if config is not None and "directory" in config:
            base_dir = config["directory"]
            output_dir = os.path.join(base_dir, "processed")
        else:
            output_dir = os.path.join(os.getcwd(), "processed")
    
    try:
        os.makedirs(output_dir, exist_ok=True)
    except Exception as e:
        return False, f"Failed to create output directory: {str(e)}"
    
    created_files = []
    
    # Generate linear and logarithmic plots
    for is_log_scale in [False, True]:
        scale_type = 'logarithmic' if is_log_scale else 'linear'
        scale_name = 'log' if is_log_scale else 'linear'
        
        print(f"Creating {scale_name} count vs surface area plot...")
        
        # Filter data for this scale
        plot_df = distribution_df[distribution_df['scale'] == scale_type].copy()
        
        if plot_df.empty:
            print(f"  Warning: No {scale_type} scale data available")
            continue
        
        # Get statistics (use number-weighted stats from the appropriate scale)
        stats = None
        if stats_dict and scale_type in stats_dict and 'number' in stats_dict[scale_type]:
            stats = stats_dict[scale_type]['number']
        
        # Create the plot (matching function signature)
        fig, fit_results = create_count_vs_surface_area_plot(plot_df, is_log_scale, stats, uniqueID, metadata)
        
        if fig is None:
            print(f"  Failed to create plot")
            continue
        
        # Save fit results to comprehensive fits file
        if fit_results:
            fit_type, fit_data = fit_results
            try:
                
                # Load existing fits file or create new one
                fits_filename = f"Fits_{uniqueID}_all.json"
                fits_path = os.path.join(output_dir, fits_filename)
                
                if os.path.exists(fits_path):
                    with open(fits_path, 'r') as f:
                        all_fits = json.load(f)
                else:
                    all_fits = {'dataset': uniqueID, 'fits': {}}
                
                # Add this fit to the collection
                fit_key = f"count_vs_surface_area_{scale_name}"
                all_fits['fits'][fit_key] = {
                    'distribution_type': 'count_vs_surface_area',
                    'scale': scale_name,
                    'fit_type': fit_type,
                    'parameters': fit_data
                }
                
                # Save updated fits file
                with open(fits_path, 'w') as f:
                    json.dump(all_fits, f, indent=2, default=str)
                
                print(f"  ✓ Saved fit to: {fits_filename}")
                
            except Exception as e:
                print(f"  ⚠ Failed to save fit: {str(e)}")
        
        # Save the plot
        try:
            base_filename = f"Plot_{uniqueID}_{scale_name}_count_vs_surface_area"
            
            pdf_path = os.path.join(output_dir, f"{base_filename}.pdf")
            fig.savefig(pdf_path, bbox_inches='tight', dpi=300)
            
            png_path = os.path.join(output_dir, f"{base_filename}.png")
            fig.savefig(png_path, bbox_inches='tight', dpi=300)
            
            created_files.append(pdf_path)
            print(f"  ✓ Saved: {base_filename}.pdf/.png")
            
            plt.close(fig)
            
        except Exception as e:
            print(f"  ✗ Failed to save plot: {str(e)}")
            plt.close(fig)
            continue
    
    return True, created_files


# Execute if we have the required data
if 'current_distribution_df' in globals() and current_distribution_df is not None:
    print("=" * 80)
    print("GENERATING RAW COUNT vs THEORETICAL SURFACE AREA PLOTS WITH D-VALUES")
    print("=" * 80)
    
    uniqueID = current_metadata.get('persistentID', 'unknown') if 'current_metadata' in globals() else 'unknown'
    stats = current_stats if 'current_stats' in globals() else None
    metadata = current_metadata if 'current_metadata' in globals() else None
    config = CONFIG if 'CONFIG' in globals() else None
    
    print(f"Creating count vs surface area plots for: {uniqueID}")
    print("Surface Area = π × diameter² (spherical particles) in µm²")
    print("Includes: Linear + Logarithmic (Lognormal fits + D-values for both)")
    
    success, plot_files = generate_count_vs_surface_area_plots(
        current_distribution_df,
        stats_dict=stats,
        uniqueID=uniqueID,
        metadata=metadata,
        config=config
    )
    
    if not success:
        print(f"ERROR: {plot_files}")
    else:
        print(f"\n✓ Successfully created {len(plot_files)} count vs surface area plots!")
        for filepath in plot_files:
            print(f"  - {os.path.basename(filepath)}")
        
        current_count_vs_surface_area_plots = plot_files

else:
    print("No data found. Run the complete workflow first.")


# ======================================================================
# CELL 14
# ======================================================================





def lognormal_pdf(x, mu, sigma, amplitude):
    """Calculate lognormal probability density function with numerical stability."""
    # Add small epsilon to avoid log(0)
    x_safe = np.maximum(x, 1e-10)
    
    # Calculate the lognormal PDF
    log_x = np.log(x_safe)
    exponent = -((log_x - mu) ** 2) / (2 * sigma ** 2)
    
    # Avoid numerical overflow/underflow
    exponent = np.clip(exponent, -50, 50)
    
    pdf = amplitude * np.exp(exponent) / (x_safe * sigma * np.sqrt(2 * np.pi))
    
    # Ensure no NaN or inf values
    pdf = np.nan_to_num(pdf, nan=0.0, posinf=0.0, neginf=0.0)
    
    return pdf



def add_d_value_lines_and_bands_volume(ax, stats):
    """Add D-value lines and uncertainty bands to a subplot, converted to volume in µm³."""
    legend_elements = []
    
    if not stats or 'D10_avg' not in stats:
        return legend_elements
    
    # Convert diameter D-values to volume D-values using V = (π/6) × d³ and convert to µm³
    d10_avg_size = stats['D10_avg']
    d10_lower_size = stats.get('D10_lower', d10_avg_size)
    d10_upper_size = stats.get('D10_upper', d10_avg_size)
    
    d50_avg_size = stats['D50_avg'] 
    d50_lower_size = stats.get('D50_lower', d50_avg_size)
    d50_upper_size = stats.get('D50_upper', d50_avg_size)
    
    d90_avg_size = stats['D90_avg']
    d90_lower_size = stats.get('D90_lower', d90_avg_size)
    d90_upper_size = stats.get('D90_upper', d90_avg_size)
    
    # Convert to volumes in µm³ (divide by 1e9 to convert from nm³ to µm³)
    d10_avg_vol = (np.pi / 6) * (d10_avg_size ** 3) / 1e9
    d10_lower_vol = (np.pi / 6) * (d10_lower_size ** 3) / 1e9
    d10_upper_vol = (np.pi / 6) * (d10_upper_size ** 3) / 1e9
    
    d50_avg_vol = (np.pi / 6) * (d50_avg_size ** 3) / 1e9
    d50_lower_vol = (np.pi / 6) * (d50_lower_size ** 3) / 1e9
    d50_upper_vol = (np.pi / 6) * (d50_upper_size ** 3) / 1e9
    
    d90_avg_vol = (np.pi / 6) * (d90_avg_size ** 3) / 1e9
    d90_lower_vol = (np.pi / 6) * (d90_lower_size ** 3) / 1e9
    d90_upper_vol = (np.pi / 6) * (d90_upper_size ** 3) / 1e9
    
    span = stats.get('span_avg', (d90_avg_size-d10_avg_size)/d50_avg_size if d50_avg_size > 0 else 0)
    
    # Add D-value lines and bands using volume values in µm³
    for d_val_vol, d_lower_vol, d_upper_vol, style, width, alpha_band in [
        (d10_avg_vol, d10_lower_vol, d10_upper_vol, '--', 1.5, 0.15),
        (d50_avg_vol, d50_lower_vol, d50_upper_vol, '-', 2.5, 0.25), 
        (d90_avg_vol, d90_lower_vol, d90_upper_vol, '--', 1.5, 0.15)
    ]:
        if not np.isnan(d_val_vol):
            ax.axvline(x=d_val_vol, color='gray', linestyle=style, alpha=0.8, linewidth=width, zorder=5)
            if not np.isnan(d_lower_vol) and not np.isnan(d_upper_vol) and (d_lower_vol != d_val_vol or d_upper_vol != d_val_vol):
                ax.axvspan(d_lower_vol, d_upper_vol, alpha=alpha_band, color='gray', zorder=1)
    
    # Create legend elements showing diameter in nm and volume in µm³
    legend_elements.extend([
        Line2D([0], [0], color='gray', linestyle='--', linewidth=1.5, 
              label=f'D10: {d10_avg_size:.1f} nm → {d10_avg_vol:.4f} µm³ ({d10_lower_size:.1f}-{d10_upper_size:.1f} nm)'),
        Line2D([0], [0], color='gray', linestyle='-', linewidth=2.5, 
              label=f'D50: {d50_avg_size:.1f} nm → {d50_avg_vol:.4f} µm³ ({d50_lower_size:.1f}-{d50_upper_size:.1f} nm)'),
        Line2D([0], [0], color='gray', linestyle='--', linewidth=1.5, 
              label=f'D90: {d90_avg_size:.1f} nm → {d90_avg_vol:.4f} µm³ ({d90_lower_size:.1f}-{d90_upper_size:.1f} nm)'),
        Line2D([0], [0], color='white', linestyle='', 
              label=f'Span: {span:.3f}')
    ])
    
    return legend_elements


def add_count_fit_curve_volume(ax, plot_df, is_log_scale, fit_color='#F25C54'):
    """No fitting for volume distributions - they're too skewed for standard fits."""
    # Return empty fit elements - no fitting needed
    return [], None


def create_count_vs_volume_plot(plot_df, is_log_scale, stats=None, uniqueID=None, metadata=None):
    """Create a two-subplot plot for raw counts vs theoretical volume."""
    
    scale_name = "Logarithmic" if is_log_scale else "Linear"
    xscale = 'log' if is_log_scale else 'linear'
    color = '#2E7D32'  # Forest green (same as cell_10b volume plots)
    
    # Calculate theoretical volume for spherical particles: V = (π/6) * d³ and convert to µm³
    plot_df = plot_df.copy()
    plot_df['volume_um3'] = (np.pi / 6) * (plot_df['size_nm'] ** 3) / 1e9  # Convert nm³ to µm³
    
    # Sort by volume
    plot_df = plot_df.sort_values('volume_um3')
    
    # Create figure
    fig = plt.figure(figsize=(7, 9))
    gs = gridspec.GridSpec(2, 1, height_ratios=[0.6, 0.4], hspace=0.3, 
                          top=0.82, bottom=0.08)
    
    # =================================================================
    # TOP SUBPLOT: MAIN DISTRIBUTION WITH ERROR BARS AND FITS
    # =================================================================
    ax1 = fig.add_subplot(gs[0])
    
    # Plot main distribution with error bars
    if 'number_sd' in plot_df.columns:
        ax1.errorbar(plot_df['volume_um3'], plot_df['number_avg'], 
                    yerr=plot_df['number_sd'],
                    fmt='o', color=color, ecolor=color, alpha=0.7,
                    capsize=3, capthick=1, markersize=6, linewidth=1.5,
                    label='Raw Counts')
    else:
        ax1.scatter(plot_df['volume_um3'], plot_df['number_avg'], 
                   color=color, s=60, alpha=0.8, label='Raw Counts')
    
    # Add lognormal fit curve (disabled for volume - too skewed)
    fit_result = add_count_fit_curve_volume(ax1, plot_df, is_log_scale)
    if isinstance(fit_result, tuple):
        fit_legend_elements, fit_results = fit_result
    else:
        fit_legend_elements = fit_result
        fit_results = None
    
    # No fitting for volume distributions - they're inherently too skewed
    
    # Format top subplot
    ax1.set_ylabel('Raw Particle Counts', color=color, fontsize=14, labelpad=10)
    ax1.tick_params(axis='y', labelcolor=color, labelsize=12)
    ax1.tick_params(axis='x', labelsize=12)
    ax1.spines['left'].set_color(color)
    
    # Set x-axis scale and add better tick labels for log scale
    ax1.set_xscale(xscale)
    if is_log_scale:
        ax1.xaxis.set_major_locator(LogLocator(base=10, numticks=12))
        ax1.xaxis.set_minor_locator(LogLocator(base=10, subs=(0.2, 0.4, 0.6, 0.8), numticks=12))
        ax1.xaxis.set_major_formatter(LogFormatter(base=10, labelOnlyBase=False))
    
    # Smart range calculation based on where the count signal is
    weights_for_range = plot_df['number_avg'].values
    vol_for_range = plot_df['volume_um3'].values
    
    # Find where 99% of the count signal is contained
    cumsum_weights = np.cumsum(weights_for_range)
    total_weight = cumsum_weights[-1]
    
    if total_weight > 0:
        # Find 1st and 99th percentiles of the count-weighted distribution
        p1_idx = np.searchsorted(cumsum_weights, 0.01 * total_weight)
        p99_idx = np.searchsorted(cumsum_weights, 0.99 * total_weight)
        
        signal_min_vol = vol_for_range[max(0, p1_idx)]
        signal_max_vol = vol_for_range[min(len(vol_for_range)-1, p99_idx)]
        data_max_vol = plot_df['volume_um3'].max()
        
        if is_log_scale:
            # Log scale: focus on the count signal range with some padding
            min_vol = max(signal_min_vol * 0.7, (np.pi / 6) * (20**3) / 1e9)  # Don't go below 20nm diameter equivalent
            max_vol = min(signal_max_vol * 1.8, data_max_vol * 1.2)
            ax1.set_xlim([min_vol, max_vol])
        else:
            # Linear scale: tighten the range, focus on main signal
            min_vol = 0
            max_vol = min(signal_max_vol * 1.15, (np.pi / 6) * (350**3) / 1e9)  # Cap at ~350nm diameter equivalent
            ax1.set_xlim([min_vol, max_vol])
    else:
        # Fallback if no signal
        if is_log_scale:
            ax1.set_xlim([(np.pi / 6) * (30**3) / 1e9, (np.pi / 6) * (300**3) / 1e9])
        else:
            ax1.set_xlim([0, (np.pi / 6) * (250**3) / 1e9])
    
    # Set y-axis to start from 0
    y_min, y_max = ax1.get_ylim()
    ax1.set_ylim([0, y_max])
    
    # Add D-value lines and bands (converted to volume)
    d_legend_elements = add_d_value_lines_and_bands_volume(ax1, stats)
    
    # Create comprehensive legend for top plot - PLACE OUTSIDE
    main_legend = [Line2D([0], [0], marker='o', color='w', markerfacecolor=color, 
                          markersize=8, label='Raw Counts')]
    
    all_legend_elements = main_legend + fit_legend_elements + d_legend_elements
    leg1 = ax1.legend(handles=all_legend_elements, fontsize=9, frameon=True, 
                     bbox_to_anchor=(1.05, 1), loc='upper left')
    leg1.get_frame().set_alpha(0.95)
    leg1.get_frame().set_edgecolor('lightgray')
    
    ax1.grid(True, linestyle='--', alpha=0.4)
    ax1.set_xlabel('')  # No x-label on top plot
    
    # =================================================================
    # BOTTOM SUBPLOT: CUMULATIVE DISTRIBUTION WITH UNCERTAINTY BANDS
    # =================================================================
    ax2 = fig.add_subplot(gs[1])
    
    # Calculate cumulative counts as percentage
    cumsum_counts = np.cumsum(plot_df['number_avg'])
    max_cumsum = cumsum_counts.iloc[-1] if len(cumsum_counts) > 0 else 1
    
    if max_cumsum > 0:
        cumsum_percentage = (cumsum_counts / max_cumsum) * 100
        ax2.plot(plot_df['volume_um3'], cumsum_percentage, '-', 
                color=color, linewidth=3, alpha=0.9, label='Cumulative %')
        
        # Add uncertainty bands if available
        if 'number_sd' in plot_df.columns:
            # Calculate cumulative uncertainty
            cumsum_sd = np.sqrt(np.cumsum(plot_df['number_sd'] ** 2))
            cumsum_sd_percentage = (cumsum_sd / max_cumsum) * 100
            ax2.fill_between(plot_df['volume_um3'], 
                           cumsum_percentage - cumsum_sd_percentage,
                           cumsum_percentage + cumsum_sd_percentage,
                           color=color, alpha=0.25, zorder=1, label='± SD')
    
    ax2.set_ylim([0, 110])
    ax2.set_ylabel('Cumulative Percentage (%)', color=color, fontsize=14, labelpad=10)
    ax2.tick_params(axis='y', labelcolor=color, labelsize=12)
    ax2.spines['left'].set_color(color)
    
    # Format bottom subplot
    ax2.set_xlabel('Theoretical Volume (µm³)', fontsize=14, labelpad=10)
    ax2.tick_params(axis='x', labelsize=12)
    ax2.set_xscale(xscale)
    ax2.set_xlim(ax1.get_xlim())  # Match top plot limits
    
    if is_log_scale:
        # Add same detailed log scale ticks to bottom plot
        ax2.xaxis.set_major_locator(LogLocator(base=10, numticks=12))
        ax2.xaxis.set_minor_locator(LogLocator(base=10, subs=(0.2, 0.4, 0.6, 0.8), numticks=12))
        ax2.xaxis.set_major_formatter(LogFormatter(base=10, labelOnlyBase=False))
    
    # Add D-value lines to bottom plot (converted to volume)
    add_d_value_lines_and_bands_volume(ax2, stats)
    
    # Legend for bottom plot - PLACE OUTSIDE
    cumulative_legend = [
        Line2D([0], [0], color=color, linewidth=3, label='Cumulative %'),
        Line2D([0], [0], color=color, alpha=0.25, linewidth=8, label='± SD')
    ]
    
    leg2 = ax2.legend(handles=cumulative_legend, fontsize=10, frameon=True, 
                     bbox_to_anchor=(1.05, 1), loc='upper left')
    leg2.get_frame().set_alpha(0.95)
    leg2.get_frame().set_edgecolor('lightgray')
    
    ax2.grid(True, linestyle='--', alpha=0.4)
    
    # =================================================================
    # TITLE AND METADATA
    # =================================================================
    
    # Extract replicate info
    replicate_info = ""
    if metadata and 'num_replicates' in metadata:
        num_reps = metadata['num_replicates']
        if num_reps and str(num_reps) != '1':
            replicate_info = f" (n={num_reps})"
    
    # Set main title
    main_title = f'{scale_name} Raw Counts vs\nTheoretical Volume: {uniqueID}{replicate_info}'
    fig.suptitle(main_title, fontsize=14, fontweight='bold', y=0.94)
    
    # Add subtitle
    subtitle = f"Error bars/bands: ± SD | Volume = (π/6) × diameter³ | No fits (volume too skewed)"
    fig.text(0.5, 0.87, subtitle, ha='center', fontsize=11, style='italic')
    
    return fig, fit_results


def generate_count_vs_volume_plots(distribution_df, stats_dict=None, uniqueID=None, 
                                  metadata=None, output_dir=None, config=None):
    """Generate raw count vs volume plots for both linear and log scales."""
    
    if distribution_df is None or distribution_df.empty:
        return False, "No data available for plotting"
    
    # Check if we have the required columns
    if 'number_avg' not in distribution_df.columns:
        return False, "Raw count data (number_avg) not found in dataset"
    
    plt.style.use('default')
    
    # Determine output directory
    if output_dir is None:
        if config is not None and "directory" in config:
            base_dir = config["directory"]
            output_dir = os.path.join(base_dir, "processed")
        else:
            output_dir = os.path.join(os.getcwd(), "processed")
    
    try:
        os.makedirs(output_dir, exist_ok=True)
    except Exception as e:
        return False, f"Failed to create output directory: {str(e)}"
    
    created_files = []
    
    # Generate linear and logarithmic plots
    for is_log_scale in [False, True]:
        scale_type = 'logarithmic' if is_log_scale else 'linear'
        scale_name = 'log' if is_log_scale else 'linear'
        
        print(f"Creating {scale_name} count vs volume plot...")
        
        # Filter data for this scale
        plot_df = distribution_df[distribution_df['scale'] == scale_type].copy()
        
        if plot_df.empty:
            print(f"  Warning: No {scale_type} scale data available")
            continue
        
        # Get statistics (use number-weighted stats from the appropriate scale)
        stats = None
        if stats_dict and scale_type in stats_dict and 'number' in stats_dict[scale_type]:
            stats = stats_dict[scale_type]['number']
        
        # Create the plot (matching function signature)
        fig, fit_results = create_count_vs_volume_plot(plot_df, is_log_scale, stats, uniqueID, metadata)
        
        if fig is None:
            print(f"  Failed to create plot")
            continue
        
        # Save fit results to comprehensive fits file
        if fit_results:
            fit_type, fit_data = fit_results
            try:
                
                # Load existing fits file or create new one
                fits_filename = f"Fits_{uniqueID}_all.json"
                fits_path = os.path.join(output_dir, fits_filename)
                
                if os.path.exists(fits_path):
                    with open(fits_path, 'r') as f:
                        all_fits = json.load(f)
                else:
                    all_fits = {'dataset': uniqueID, 'fits': {}}
                
                # Add this fit to the collection
                fit_key = f"count_vs_volume_{scale_name}"
                all_fits['fits'][fit_key] = {
                    'distribution_type': 'count_vs_volume',
                    'scale': scale_name,
                    'fit_type': fit_type,
                    'parameters': fit_data
                }
                
                # Save updated fits file
                with open(fits_path, 'w') as f:
                    json.dump(all_fits, f, indent=2, default=str)
                
                print(f"  ✓ Saved fit to: {fits_filename}")
                
            except Exception as e:
                print(f"  ⚠ Failed to save fit: {str(e)}")
        
        # Save the plot
        try:
            base_filename = f"Plot_{uniqueID}_{scale_name}_count_vs_volume"
            
            pdf_path = os.path.join(output_dir, f"{base_filename}.pdf")
            fig.savefig(pdf_path, bbox_inches='tight', dpi=300)
            
            png_path = os.path.join(output_dir, f"{base_filename}.png")
            fig.savefig(png_path, bbox_inches='tight', dpi=300)
            
            created_files.append(pdf_path)
            print(f"  ✓ Saved: {base_filename}.pdf/.png")
            
            plt.close(fig)
            
        except Exception as e:
            print(f"  ✗ Failed to save plot: {str(e)}")
            plt.close(fig)
            continue
    
    return True, created_files


# Execute if we have the required data
if 'current_distribution_df' in globals() and current_distribution_df is not None:
    print("=" * 80)
    print("GENERATING RAW COUNT vs THEORETICAL VOLUME PLOTS WITH D-VALUES")
    print("=" * 80)
    
    uniqueID = current_metadata.get('persistentID', 'unknown') if 'current_metadata' in globals() else 'unknown'
    stats = current_stats if 'current_stats' in globals() else None
    metadata = current_metadata if 'current_metadata' in globals() else None
    config = CONFIG if 'CONFIG' in globals() else None
    
    print(f"Creating count vs volume plots for: {uniqueID}")
    print("Volume = (π/6) × diameter³ (spherical particles) in µm³")
    print("Includes: Linear + Logarithmic (Lognormal fits + D-values for both)")
    
    success, plot_files = generate_count_vs_volume_plots(
        current_distribution_df,
        stats_dict=stats,
        uniqueID=uniqueID,
        metadata=metadata,
        config=config
    )
    
    if not success:
        print(f"ERROR: {plot_files}")
    else:
        print(f"\n✓ Successfully created {len(plot_files)} count vs volume plots!")
        for filepath in plot_files:
            print(f"  - {os.path.basename(filepath)}")
        
        current_count_vs_volume_plots = plot_files

else:
    print("No data found. Run the complete workflow first.")
