"""
Download Manager for NTA Analysis
Handles file organization, ZIP creation, and fit data export
"""

import os
import json
import pandas as pd
import numpy as np
from pathlib import Path
import zipfile
import io


def create_fit_data_tsv(distribution_df, stats_dict, scale_type='linear', dist_type='number'):
    """
    Create a simple TSV with distribution data + fit curve for a specific plot.
    
    Parameters:
    - distribution_df: Full distribution dataframe
    - stats_dict: Statistics dictionary with fit info
    - scale_type: 'linear' or 'logarithmic'
    - dist_type: 'number', 'volume', or 'surface_area'
    
    Returns:
    - DataFrame with: size_nm, measured_count, fitted_curve, fit_parameters_as_comment
    """
    
    # Filter to this scale
    plot_df = distribution_df[distribution_df['scale'] == scale_type].copy()
    
    if plot_df.empty:
        return None
    
    # Get the appropriate column
    if dist_type == 'number':
        measured_col = 'number_normalized_avg'
        measured_sd_col = 'number_normalized_sd'
    elif dist_type == 'volume':
        measured_col = 'volume_nm^3_per_mL_avg'
        measured_sd_col = 'volume_nm^3_per_mL_sd'
    elif dist_type == 'surface_area':
        measured_col = 'area_nm^2_per_mL_avg'
        measured_sd_col = 'area_nm^2_per_mL_sd'
    else:
        return None
    
    if measured_col not in plot_df.columns:
        return None
    
    # Sort by size
    plot_df = plot_df.sort_values('size_nm').reset_index(drop=True)
    
    # Create output dataframe
    result_df = pd.DataFrame({
        'size_nm': plot_df['size_nm'].values,
        f'{dist_type}_measured': plot_df[measured_col].values,
        f'{dist_type}_sd': plot_df.get(measured_sd_col, pd.Series([0]*len(plot_df))).values,
    })
    
    # Try to add fit curve if available
    try:
        from scipy.optimize import curve_fit
        
        sizes = plot_df['size_nm'].values
        weights = plot_df[measured_col].values
        
        # Simple lognormal fit
        def lognormal_pdf(x, mu, sigma, amplitude):
            return amplitude * (1 / (x * sigma * np.sqrt(2 * np.pi))) * \
                   np.exp(-(np.log(x) - mu)**2 / (2 * sigma**2))
        
        # Fit
        valid_mask = (weights > 0) & (sizes > 0)
        if np.sum(valid_mask) > 3:
            sizes_valid = sizes[valid_mask]
            weights_valid = weights[valid_mask]
            
            size_log = np.log(sizes_valid)
            mu_init = np.average(size_log, weights=weights_valid)
            sigma_init = np.sqrt(np.average((size_log - mu_init)**2, weights=weights_valid))
            amplitude_init = np.max(weights_valid) * sigma_init * np.sqrt(2 * np.pi) * np.exp(mu_init)
            
            params, _ = curve_fit(
                lognormal_pdf, sizes_valid, weights_valid,
                p0=[mu_init, sigma_init, amplitude_init],
                bounds=([0, 0, 0], [np.inf, np.inf, np.inf]),
                maxfev=10000
            )
            
            # Generate fitted curve
            fit_curve = lognormal_pdf(sizes, *params)
            result_df[f'{dist_type}_fitted'] = fit_curve
            
            # Store parameters as metadata
            result_df.attrs['fit_type'] = 'lognormal'
            result_df.attrs['fit_mu'] = float(params[0])
            result_df.attrs['fit_sigma'] = float(params[1])
            result_df.attrs['fit_amplitude'] = float(params[2])
            result_df.attrs['geo_mean'] = float(np.exp(params[0]))
            result_df.attrs['geo_std'] = float(np.exp(params[1]))
    except:
        pass
    
    return result_df


def export_fit_data_tsv(distribution_df, stats_dict, scale_type='linear', dist_type='number'):
    """
    Create TSV text with fit data and parameters as comments.
    
    Returns string suitable for download.
    """
    fit_df = create_fit_data_tsv(distribution_df, stats_dict, scale_type, dist_type)
    
    if fit_df is None:
        return None
    
    # Create header with parameters
    lines = []
    lines.append(f"# Fit Data for {dist_type.replace('_', ' ').title()} Distribution ({scale_type.title()} Scale)")
    lines.append("#")
    
    if 'fit_type' in fit_df.attrs:
        lines.append(f"# Fit Type: {fit_df.attrs['fit_type']}")
        lines.append(f"# Geometric Mean: {fit_df.attrs.get('geo_mean', 'N/A'):.2f} nm")
        lines.append(f"# Geometric Std Dev: {fit_df.attrs.get('geo_std', 'N/A'):.2f}")
        lines.append(f"# Parameters: mu={fit_df.attrs.get('fit_mu', 'N/A'):.4f}, "
                    f"sigma={fit_df.attrs.get('fit_sigma', 'N/A'):.4f}, "
                    f"amplitude={fit_df.attrs.get('fit_amplitude', 'N/A'):.4e}")
        lines.append("#")
    
    # Get TSV content
    tsv_content = fit_df.to_csv(sep='\t', index=False)
    
    return '\n'.join(lines) + '\n' + tsv_content


def create_download_zip(output_dir, unique_id, num_replicates, metadata_dict=None, 
                       distribution_df=None, statistics_dict=None):
    """
    Create a ZIP file with all analysis results.
    
    Includes:
    - Distribution data (CSV)
    - Metadata (TSV)
    - Statistics (TSV)
    - Fit data (TSV)
    - Plot PDFs (skips PNGs)
    
    Returns: bytes for download
    """
    zip_buffer = io.BytesIO()
    
    with zipfile.ZipFile(zip_buffer, 'w', zipfile.ZIP_DEFLATED) as zip_file:
        output_path = Path(output_dir)
        
        # Add data files section
        # 1. Distribution data
        if distribution_df is not None:
            try:
                dist_csv = distribution_df.to_csv(index=False)
                zip_file.writestr(
                    f"data/distribution_{unique_id}_avg{num_replicates}.csv",
                    dist_csv
                )
            except Exception as e:
                print(f"  ⚠ Could not add distribution data: {str(e)}")
        
        # 2. Metadata
        if metadata_dict is not None:
            try:
                metadata_df = pd.DataFrame(
                    [(k, v) for k, v in metadata_dict.items()],
                    columns=['Field', 'Value']
                )
                metadata_tsv = metadata_df.to_csv(sep='\t', index=False)
                zip_file.writestr(
                    f"data/metadata_{unique_id}_avg{num_replicates}.txt",
                    metadata_tsv
                )
            except Exception as e:
                print(f"  ⚠ Could not add metadata: {str(e)}")
        
        # 3. Statistics
        if statistics_dict is not None:
            try:
                stats_list = []
                for scale in statistics_dict:
                    if isinstance(statistics_dict[scale], dict):
                        for dist_type in statistics_dict[scale]:
                            stat_dict = statistics_dict[scale][dist_type]
                            row = {
                                'Scale': scale,
                                'Distribution': dist_type,
                                **stat_dict
                            }
                            stats_list.append(row)
                
                if stats_list:
                    stats_df = pd.DataFrame(stats_list)
                    stats_tsv = stats_df.to_csv(sep='\t', index=False)
                    zip_file.writestr(
                        f"data/statistics_{unique_id}_avg{num_replicates}.txt",
                        stats_tsv
                    )
            except Exception as e:
                print(f"  ⚠ Could not add statistics: {str(e)}")
        
        # Add plot PDFs (skip PNGs)
        if output_path.exists():
            pdf_files = sorted(output_path.glob('*Plot_*.pdf'))
            if pdf_files:
                for pdf_file in pdf_files:
                    try:
                        zip_file.write(pdf_file, arcname=f"plots/{pdf_file.name}")
                    except Exception as e:
                        print(f"  ⚠ Could not add plot: {str(e)}")
            
            # Add fit data TSV files
            fit_files = sorted(output_path.glob('FitData_*.txt'))
            if fit_files:
                for fit_file in fit_files:
                    try:
                        zip_file.write(fit_file, arcname=f"fit_data/{fit_file.name}")
                    except Exception as e:
                        print(f"  ⚠ Could not add fit data: {str(e)}")
    
    zip_buffer.seek(0)
    return zip_buffer.getvalue()


def get_all_output_files(output_dir):
    """
    Get all generated files organized by type.
    
    Returns dict with lists of files by category.
    """
    output_path = Path(output_dir)
    
    files = {
        'data': [],
        'plots_pdf': [],
        'plots_png': [],
        'fits': [],
        'metadata': []
    }
    
    if not output_path.exists():
        return files
    
    for file_path in output_path.glob('*'):
        if not file_path.is_file():
            continue
        
        filename = file_path.name
        
        if 'distribution' in filename and filename.endswith('.txt'):
            files['data'].append((filename, str(file_path)))
        elif 'Plot_' in filename and filename.endswith('.pdf'):
            files['plots_pdf'].append((filename, str(file_path)))
        elif 'Plot_' in filename and filename.endswith('.png'):
            files['plots_png'].append((filename, str(file_path)))
        elif 'statistics' in filename and filename.endswith('.txt'):
            files['data'].append((filename, str(file_path)))
        elif 'metadata' in filename and filename.endswith('.txt'):
            files['metadata'].append((filename, str(file_path)))
        elif 'Fits_' in filename or '_fits_' in filename:
            files['fits'].append((filename, str(file_path)))
    
    return files


def format_file_list(files_dict):
    """Format file list for display."""
    display = {}
    for category, files in files_dict.items():
        if files:
            display[category.replace('_', ' ').title()] = len(files)
    return display


class DownloadPackage:
    """
    Manage a complete download package with all files.
    """
    
    def __init__(self, unique_id, num_replicates, output_dir):
        self.unique_id = unique_id
        self.num_replicates = num_replicates
        self.output_dir = output_dir
        self.files = get_all_output_files(output_dir)
    
    def get_zip_bytes(self):
        """Get ZIP file as bytes for download."""
        return create_download_zip(self.output_dir, self.unique_id, self.num_replicates)
    
    def get_file_count(self):
        """Get total count of files."""
        count = 0
        for files_list in self.files.values():
            count += len(files_list)
        return count
    
    def list_files(self):
        """Get organized list of files."""
        return self.files
    
    def get_file_summary(self):
        """Get summary of available files."""
        summary = {}
        summary['data_files'] = len(self.files['data'])
        summary['plot_pdfs'] = len(self.files['plots_pdf'])
        summary['plot_pngs'] = len(self.files['plots_png'])
        summary['fit_files'] = len(self.files['fits'])
        summary['metadata_files'] = len(self.files['metadata'])
        summary['total'] = self.get_file_count()
        return summary
