"""
NTA Data Analysis - Simplified Standalone Script

Adapted for actual ZetaView output format.
Handles single or multiple file analysis with averaging and uncertainty.
"""

import os
import re
import json
from datetime import date

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import LogFormatter

from scipy.optimize import curve_fit


# ============================================================================
# FILE READING & PARSING
# ============================================================================

def read_nta_file(filepath):
    """Read NTA file with proper encoding."""
    try:
        with open(filepath, 'r', encoding='latin1') as f:
            return f.read()
    except Exception as e:
        raise Exception(f"Error reading file: {str(e)}")


def extract_metadata(content, filename):
    """Extract metadata from file header."""
    metadata = {
        'filename': filename,
        'extraction_date': str(date.today()),
    }
    
    # Parse header lines
    lines = content.split('\n')
    for line in lines[:100]:
        if ':' in line:
            parts = line.split(':', 1)
            if len(parts) == 2:
                key = parts[0].strip()
                value = parts[1].strip()
                
                # Extract important fields
                if key == 'Sample':
                    metadata['sample_name'] = value
                elif key == 'Date':
                    metadata['measurement_date'] = value
                elif key == 'Dilution::' or key == 'Dilution':
                    try:
                        metadata['dilution_factor'] = float(value)
                    except:
                        metadata['dilution_factor'] = 1.0
                elif key == 'Temperature':
                    try:
                        metadata['temperature'] = float(value.split()[0])
                    except:
                        pass
    
    # Extract persistent ID from filename
    id_match = re.search(r'([A-Za-z0-9_]+)_\d+_rawdata', filename)
    metadata['persistentID'] = id_match.group(1) if id_match else filename.replace('.txt', '')
    
    return metadata


def extract_distribution_data(content):
    """Extract particle size distribution from file content."""
    
    # Find the data section
    lines = content.split('\n')
    
    data_start = None
    header_line = None
    
    for i, line in enumerate(lines):
        if 'Size / nm' in line and 'Number' in line:
            data_start = i
            header_line = line
            break
    
    if data_start is None:
        raise Exception("Could not find data section in file")
    
    # Parse header
    headers = [h.strip() for h in header_line.split('\t')]
    
    # Parse data rows
    data_rows = []
    for line in lines[data_start + 1:]:
        if not line.strip() or not any(c.isdigit() for c in line):
            break
        
        try:
            values = line.split('\t')
            if len(values) >= len(headers):
                row = {}
                for header, value in zip(headers, values[:len(headers)]):
                    try:
                        row[header.strip()] = float(value.strip())
                    except:
                        row[header.strip()] = None
                
                # Only keep rows with positive size and some data
                size = row.get('Size / nm')
                if size is not None and size > 0 and any(v is not None and v > 0 for v in row.values() if v != size):
                    data_rows.append(row)
        except:
            continue
    
    if not data_rows:
        raise Exception("No data rows found in file")
    
    df = pd.DataFrame(data_rows)
    
    # Standardize column names
    df.columns = [col.strip() for col in df.columns]
    
    return df


# ============================================================================
# DATA PROCESSING
# ============================================================================

def process_files(filepaths, dilution_override=None):
    """Process one or more NTA files with averaging."""
    
    if isinstance(filepaths, str):
        filepaths = [filepaths]
    
    # Read all files
    distributions = []
    metadatas = []
    filenames = []
    
    for filepath in filepaths:
        content = read_nta_file(filepath)
        filename = os.path.basename(filepath)
        
        metadata = extract_metadata(content, filename)
        dist = extract_distribution_data(content)
        
        distributions.append(dist)
        metadatas.append(metadata)
        filenames.append(filename)
    
    # Apply dilution override if provided
    if dilution_override is not None:
        for meta in metadatas:
            meta['dilution_factor'] = dilution_override
    
    # Average if multiple files
    if len(distributions) > 1:
        # Sort all by size and align them
        for dist in distributions:
            dist.sort_values('Size / nm', inplace=True)
        
        # Use first distribution as baseline
        avg_dist = distributions[0].copy()
        avg_dist = avg_dist.reset_index(drop=True)
        
        # Average numeric columns
        for col in avg_dist.columns:
            if col not in ['Size / nm', 'size_nm'] and pd.api.types.is_numeric_dtype(avg_dist[col]):
                col_values = []
                
                for dist in distributions:
                    if col in dist.columns:
                        col_values.append(dist[col].values)
                
                if col_values:
                    # Handle different lengths by averaging aligned values
                    try:
                        col_array = np.array(col_values)
                        avg_dist[col] = np.nanmean(col_array, axis=0)
                        
                        # Store std dev
                        if len(distributions) > 1:
                            avg_dist[f'{col}_std'] = np.nanstd(col_array, axis=0, ddof=1)
                    except:
                        # If arrays are different lengths, use first as base
                        avg_dist[col] = distributions[0][col].values
                        if len(distributions) > 1:
                            avg_dist[f'{col}_std'] = 0.0
        
        avg_dist['num_replicates'] = len(distributions)
    else:
        avg_dist = distributions[0].copy()
        avg_dist['num_replicates'] = 1
        
        # Add zero std for consistency
        for col in avg_dist.columns:
            if col not in ['Size / nm', 'size_nm', 'num_replicates'] and pd.api.types.is_numeric_dtype(avg_dist[col]):
                avg_dist[f'{col}_std'] = 0.0
    
    # Use first metadata as base
    metadata = metadatas[0].copy()
    metadata['num_files'] = len(distributions)
    
    return avg_dist, metadata, filenames


def apply_dilution_correction(df, dilution_factor=1.0):
    """Apply dilution factor to concentration columns."""
    df = df.copy()
    
    for col in df.columns:
        if 'Concentration' in col:
            df[col] = df[col] * dilution_factor
            # Also scale uncertainty if it exists
            std_col = f'{col}_std'
            if std_col in df.columns:
                df[std_col] = df[std_col] * dilution_factor
    
    return df


def normalize_distribution(df):
    """Normalize to unit area (for number-weighted distribution)."""
    df = df.copy()
    
    if 'Size / nm' in df.columns and 'Number' in df.columns:
        sizes = df['Size / nm'].values
        numbers = df['Number'].values
        
        # Trapezoidal integration
        from scipy.integrate import trapezoid
        integral = trapezoid(numbers, sizes)
        
        if integral > 0:
            df['Number_normalized'] = numbers / integral
            
            if 'Number_std' in df.columns:
                df['Number_normalized_std'] = df['Number_std'] / integral
    
    return df


def calculate_cumulative(df):
    """Calculate cumulative distribution."""
    df = df.copy()
    
    # Sort by size
    df = df.sort_values('Size / nm').reset_index(drop=True)
    
    if 'Number_normalized' in df.columns:
        df['Number_cumulative'] = df['Number_normalized'].cumsum()
        
        if 'Number_normalized_std' in df.columns:
            df['Number_cumulative_std'] = np.sqrt((df['Number_normalized_std'] ** 2).cumsum())
    
    return df


# ============================================================================
# PLOTTING
# ============================================================================

def lognormal_pdf(x, mu, sigma, amplitude):
    """Lognormal probability density function."""
    return (amplitude / (x * sigma * np.sqrt(2 * np.pi))) * np.exp(-(np.log(x) - mu) ** 2 / (2 * sigma ** 2))


def fit_lognormal(sizes, numbers):
    """Fit lognormal to number distribution."""
    try:
        mu_init = np.mean(np.log(sizes))
        sigma_init = np.std(np.log(sizes))
        amplitude_init = np.max(numbers)
        
        # Filter out zero values for fitting
        mask = (numbers > 0) & (sizes > 0)
        sizes_fit = sizes[mask]
        numbers_fit = numbers[mask]
        
        if len(sizes_fit) < 3:
            return None
        
        popt, _ = curve_fit(
            lognormal_pdf, sizes_fit, numbers_fit,
            p0=[mu_init, sigma_init, amplitude_init],
            maxfev=10000
        )
        return popt
    except:
        return None


def create_linear_number_plot(df, figsize=(12, 8)):
    """Create linear scale number-weighted distribution plot."""
    
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=figsize, sharex=True)
    
    if 'Size / nm' not in df.columns or 'Number_normalized' not in df.columns:
        raise Exception("Required columns not found in distribution data")
    
    sizes = df['Size / nm'].values
    numbers_norm = df['Number_normalized'].values
    
    # Plot 1: Distribution
    ax1.bar(sizes, numbers_norm, width=np.diff(sizes).mean(), alpha=0.6, label='Number distribution')
    
    # Try lognormal fit
    fit_popt = fit_lognormal(sizes, numbers_norm)
    if fit_popt is not None:
        sizes_fit = np.linspace(sizes.min(), sizes.max(), 500)
        fit_curve = lognormal_pdf(sizes_fit, *fit_popt)
        ax1.plot(sizes_fit, fit_curve, 'r-', linewidth=2, label='Lognormal fit')
    
    if 'Number_normalized_std' in df.columns:
        stds = df['Number_normalized_std'].values
        ax1.fill_between(sizes, numbers_norm - stds, numbers_norm + stds, alpha=0.3)
    
    ax1.set_ylabel('Normalized Number Distribution', fontsize=11)
    ax1.set_title('Linear Scale - Number-Weighted Particle Size Distribution', fontsize=12, fontweight='bold')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # Plot 2: Cumulative
    if 'Number_cumulative' in df.columns:
        cum = df['Number_cumulative'].values
        ax2.plot(sizes, cum, 'b-', linewidth=2, label='Cumulative')
        
        if 'Number_cumulative_std' in df.columns:
            cum_std = df['Number_cumulative_std'].values
            ax2.fill_between(sizes, cum - cum_std, cum + cum_std, alpha=0.3, label='±1 SD')
        
        ax2.axhline(y=0.5, color='r', linestyle='--', alpha=0.5, label='D50')
        ax2.set_ylabel('Cumulative Distribution', fontsize=11)
        ax2.set_xlabel('Size (nm)', fontsize=11)
        ax2.legend()
        ax2.grid(True, alpha=0.3)
        ax2.set_ylim([0, 1.05])
    
    plt.tight_layout()
    return fig


# ============================================================================
# MAIN ANALYZER CLASS
# ============================================================================

class NTAAnalyzer:
    """Main NTA analysis class."""
    
    def __init__(self, config=None):
        self.config = config or {}
        self.results = {}
    
    def process(self, filepaths, metadata_overrides=None):
        """Process NTA files and calculate everything."""
        
        # Parse files
        dist, metadata, filenames = process_files(
            filepaths,
            dilution_override=metadata_overrides.get('dilution_factor') if metadata_overrides else None
        )
        
        # Update metadata with user inputs
        if metadata_overrides:
            metadata.update(metadata_overrides)
        
        # Apply corrections
        dist = apply_dilution_correction(dist, metadata.get('dilution_factor', 1.0))
        dist = normalize_distribution(dist)
        dist = calculate_cumulative(dist)
        
        self.results = {
            'distribution': dist,
            'metadata': metadata,
            'filenames': filenames
        }
        
        return self.results
    
    def get_plot(self):
        """Generate the linear number plot."""
        if 'distribution' not in self.results:
            raise Exception("No results. Run process() first.")
        
        return create_linear_number_plot(self.results['distribution'])
    
    def save_outputs(self, output_dir):
        """Save all outputs."""
        os.makedirs(output_dir, exist_ok=True)
        
        unique_id = self.results['metadata'].get('persistentID', 'analysis')
        created_files = []
        
        # Save metadata
        meta_path = os.path.join(output_dir, f'Data_{unique_id}_metadata.txt')
        with open(meta_path, 'w') as f:
            for key, value in self.results['metadata'].items():
                f.write(f"{key}: {value}\n")
        created_files.append(meta_path)
        
        # Save distribution
        dist_path = os.path.join(output_dir, f'Data_{unique_id}_PSD.txt')
        self.results['distribution'].to_csv(dist_path, sep='\t', index=False)
        created_files.append(dist_path)
        
        # Save plot
        plot_path = os.path.join(output_dir, f'Plot_{unique_id}_linear_number.png')
        fig = self.get_plot()
        fig.savefig(plot_path, dpi=150, bbox_inches='tight')
        plt.close(fig)
        created_files.append(plot_path)
        
        # Save basic stats
        stats_path = os.path.join(output_dir, f'Stats_{unique_id}.txt')
        with open(stats_path, 'w') as f:
            dist = self.results['distribution']
            if 'Number_cumulative' in dist.columns:
                # Find D50 (median)
                cum = dist['Number_cumulative'].values
                idx = np.argmin(np.abs(cum - 0.5))
                d50 = dist['Size / nm'].values[idx]
                f.write(f"D50 (nm): {d50:.2f}\n")
        created_files.append(stats_path)
        
        return created_files
