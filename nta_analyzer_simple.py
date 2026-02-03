"""
NTA Data Analysis - Standalone Script (Updated with your actual code)

Uses critical functions directly from your original notebook:
- Cell 02: File I/O (read, identify sections, validate)
- Cell 03: Data extraction & averaging
- Cell 04: Metadata extraction

With professional plotting and analysis.
"""

import os
import re
import json
from datetime import date

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from scipy.optimize import curve_fit
from scipy.integrate import trapezoid

# ============================================================================
# IMPORTS & CONFIG
# ============================================================================

DEFAULT_CONFIG = {
    "nta_concentration_calibration_factor": 4.61E+5,
    "project_metadata": {
        "experimenter": "Your_Initials",
        "location": "Your_Lab_Location",
        "project": "Your_Project_Name",
        "meta_version": "v03",
        "pi": "Principal_Investigator_Initials",
        "funding": "Funding_Source",
        "data_collection_method": "NTA",
    }
}

# ============================================================================
# CELL 02 - FILE I/O FUNCTIONS (from your notebook)
# ============================================================================

def read_nta_file(filepath):
    """
    Read an NTA data file with appropriate encoding.
    
    Parameters:
    filepath (str): Path to the NTA file
    
    Returns:
    tuple: (success_flag, result)
        - If successful: (True, file_content)
        - If failed: (False, error_message)
    """
    try:
        with open(filepath, 'r', encoding='latin1') as file:
            content = file.read()
        
        if not content or len(content) < 100:
            return False, f"File appears to be empty or too small: {filepath}"
        
        return True, content
    except Exception as e:
        return False, f"Error reading file: {str(e)}"


def identify_sections(content):
    """
    Identify key data sections in the file content.
    
    Parameters:
    content (str): File content to analyze
    
    Returns:
    tuple: (success_flag, result)
        - If successful: (True, dict_of_sections)
        - If failed: (False, error_message)
    """
    sections = {}
    
    # Find linear size distribution section
    size_lin_start = content.find("Size Distribution")
    if size_lin_start == -1:
        return False, "Could not find 'Size Distribution' section for linear data"
    sections['linear_start'] = size_lin_start
    
    # Find logarithmic data section (starts with -1.000E+0 separator)
    size_log_start = content.find("-1.000E+0")
    if size_log_start == -1:
        # Try alternative approach - look for second "Size / nm" header
        second_header = content.find("Size / nm", size_lin_start + 100)
        if second_header == -1:
            return False, "Could not find logarithmic data section"
        sections['logarithmic_start'] = second_header
    else:
        sections['logarithmic_start'] = size_log_start
    
    # Validate section order
    if sections['linear_start'] >= sections['logarithmic_start']:
        return False, "Invalid file structure: linear section should come before logarithmic section"
    
    return True, sections


def validate_file_structure(content, sections):
    """
    Validate file structure to ensure it can be processed.
    
    Parameters:
    content (str): File content
    sections (dict): Section positions from identify_sections
    
    Returns:
    tuple: (success_flag, message)
    """
    # Check linear section has expected header pattern
    linear_section = content[sections['linear_start']:sections['logarithmic_start']]
    if not re.search(r'Size / nm\s+Number\s+Concentration', linear_section):
        return False, "Missing expected header pattern in linear section"
    
    # Check logarithmic section has data in expected format
    log_section = content[sections['logarithmic_start']:]
    if not re.search(r'[\d.-]+E[\+\-]\d+\s+[\d.-]+E[\+\-]\d+', log_section):
        return False, "Could not find data rows in logarithmic section"
    
    return True, "File structure is valid"


# ============================================================================
# CELL 03 - DATA EXTRACTION FUNCTIONS (from your notebook)
# ============================================================================

def extract_data_section(content, start_pos, end_pos=None, is_log_section=False):
    """
    Extract tabular data from a section of the file content.
    
    Parameters:
    content (str): File content
    start_pos (int): Starting position of the section
    end_pos (int, optional): Ending position of the section
    is_log_section (bool): Whether this is the logarithmic section
    
    Returns:
    tuple: (success_flag, result)
        - If successful: (True, (header, data_lines))
        - If failed: (False, error_message)
    """
    # Extract the section
    section = content[start_pos:end_pos]
    lines = section.split('\n')
    
    # Find header line
    header_line = None
    data_lines = []
    
    for i, line in enumerate(lines):
        # Look for line containing "Size / nm"
        if 'Size / nm' in line:
            header_line = line
            # Data starts from next line
            data_lines = lines[i+1:]
            break
    
    if header_line is None:
        return False, "Could not find header line in section"
    
    if not data_lines:
        return False, "No data lines found after header"
    
    return True, (header_line, data_lines)


def parse_data_lines(header_line, data_lines, scale_type):
    """
    Parse data lines into a structured DataFrame.
    
    Parameters:
    header_line (str): Header line with column names
    data_lines (list): List of data lines
    scale_type (str): 'linear' or 'logarithmic'
    
    Returns:
    pd.DataFrame: Parsed data, or None if parsing fails
    """
    # Parse header
    headers = [h.strip() for h in header_line.split('\t')]
    
    # Parse data rows
    data = []
    for line in data_lines:
        if not line.strip():
            continue
        
        try:
            values = line.split('\t')
            if len(values) >= len(headers):
                row = {}
                for header, value in zip(headers, values[:len(headers)]):
                    try:
                        row[header.strip()] = float(value.strip())
                    except (ValueError, TypeError):
                        pass
                
                # Only keep rows with valid size data
                if 'Size / nm' in row and row['Size / nm'] > 0:
                    data.append(row)
        except Exception as e:
            continue
    
    if not data:
        return None
    
    df = pd.DataFrame(data)
    df['scale'] = scale_type
    return df


def extract_single_file_distribution(content, sections, filename):
    """
    Extract complete particle distribution from a single file.
    
    Extracts both linear and logarithmic scale data.
    
    Parameters:
    content (str): File content
    sections (dict): Section positions from identify_sections
    filename (str): Filename for reference
    
    Returns:
    pd.DataFrame: Combined distribution data, or None if fails
    """
    dfs = []
    
    # Extract linear scale (linear_start to logarithmic_start)
    success, result = extract_data_section(content, sections['linear_start'], 
                                          sections['logarithmic_start'], False)
    if success:
        header, data_lines = result
        df_linear = parse_data_lines(header, data_lines, 'linear')
        if df_linear is not None:
            dfs.append(df_linear)
    
    # Extract logarithmic scale (logarithmic_start to end)
    success, result = extract_data_section(content, sections['logarithmic_start'], 
                                          None, True)
    if success:
        header, data_lines = result
        df_log = parse_data_lines(header, data_lines, 'logarithmic')
        if df_log is not None:
            dfs.append(df_log)
    
    if dfs:
        result_df = pd.concat(dfs, ignore_index=True)
        result_df['filename'] = filename
        return result_df
    
    return None


def average_replicate_data(dataframes_list, filenames_list):
    """
    Average particle distribution data across replicates.
    
    Performs bin-by-bin averaging with standard deviation calculation:
    - Mean: x̄ᵢ = (1/n)Σⱼ xᵢⱼ
    - Std Dev: sᵢ = √[(1/(n-1))Σⱼ(xᵢⱼ - x̄ᵢ)²]
    
    Parameters:
    dataframes_list (list): List of DataFrames from different files
    filenames_list (list): List of filenames for reference
    
    Returns:
    pd.DataFrame: Averaged distribution with uncertainty
    """
    if not dataframes_list:
        return None
    
    # If only one file, return as-is with zero uncertainties
    if len(dataframes_list) == 1:
        df = dataframes_list[0].copy()
        for col in df.columns:
            if col not in ['Size / nm', 'scale', 'filename'] and '_std' not in col:
                if pd.api.types.is_numeric_dtype(df[col]):
                    df[f'{col}_std'] = 0.0
        return df
    
    # Multiple files: average by size bin
    all_data = pd.concat(dataframes_list, ignore_index=True)
    
    # Group by size
    grouped = all_data.groupby('Size / nm')
    
    result_data = []
    for size, group in grouped:
        row = {'Size / nm': size}
        
        # Average numeric columns
        for col in all_data.columns:
            if col not in ['Size / nm', 'scale', 'filename'] and pd.api.types.is_numeric_dtype(all_data[col]):
                values = group[col].dropna()
                if len(values) > 0:
                    row[col] = values.mean()
                    # Calculate standard deviation
                    if len(values) > 1:
                        row[f'{col}_std'] = values.std(ddof=1)
                    else:
                        row[f'{col}_std'] = 0.0
        
        result_data.append(row)
    
    result_df = pd.DataFrame(result_data)
    if 'scale' in all_data.columns:
        result_df['scale'] = all_data['scale'].iloc[0]
    
    return result_df


# ============================================================================
# CELL 04 - METADATA EXTRACTION (from your notebook)
# ============================================================================

def extract_all_metadata_fields(content, filename):
    """
    Extract ALL possible metadata fields from file content.
    
    Uses comprehensive regex patterns to pull measurement parameters,
    settings, and instrument info.
    
    Parameters:
    content (str): File content
    filename (str): Filename
    
    Returns:
    dict: Extracted metadata
    """
    metadata = {
        'filename': filename,
        'extraction_date': str(date.today()),
    }
    
    # Split into lines for parsing
    lines = content.split('\n')
    
    # Parse header lines (key: value pairs)
    for line in lines[:150]:
        if ':' not in line:
            continue
        
        parts = line.split(':', 1)
        if len(parts) != 2:
            continue
        
        key = parts[0].strip()
        value = parts[1].strip()
        
        # Map common fields
        field_mapping = {
            'Sample': 'sample_name',
            'Date': 'measurement_date',
            'Time': 'measurement_time',
            'Temperature': 'temperature',
            'Dilution::': 'dilution_factor',
            'Dilution': 'dilution_factor',
            'pH': 'pH',
            'Conductivity': 'conductivity',
            'Viscosity': 'viscosity',
            'ZetaView S/N': 'zetaview_sn',
            'Cell S/N': 'cell_sn',
            'Software': 'software_version',
        }
        
        for key_pattern, meta_key in field_mapping.items():
            if key_pattern in key:
                # Try to extract numeric value if applicable
                try:
                    if 'factor' in meta_key or 'temperature' in meta_key.lower():
                        metadata[meta_key] = float(value.split()[0])
                    else:
                        metadata[meta_key] = value
                except (ValueError, IndexError):
                    metadata[meta_key] = value
                break
    
    # Extract persistent ID from filename
    id_match = re.search(r'([A-Za-z0-9_]+)_\d+_rawdata', filename)
    if id_match:
        metadata['persistentID'] = id_match.group(1)
    else:
        metadata['persistentID'] = filename.replace('.txt', '')
    
    return metadata


# ============================================================================
# PLOTTING FUNCTIONS
# ============================================================================

def lognormal_pdf(x, mu, sigma, amplitude):
    """Lognormal probability density function."""
    return (amplitude / (x * sigma * np.sqrt(2 * np.pi))) * np.exp(-(np.log(x) - mu) ** 2 / (2 * sigma ** 2))


def fit_lognormal(sizes, numbers):
    """Fit lognormal to number distribution."""
    try:
        mask = (numbers > 0) & (sizes > 0)
        sizes_fit = sizes[mask]
        numbers_fit = numbers[mask]
        
        if len(sizes_fit) < 5:
            return None
        
        mu_init = np.mean(np.log(sizes_fit))
        sigma_init = np.std(np.log(sizes_fit))
        amplitude_init = np.max(numbers_fit)
        
        popt, _ = curve_fit(
            lognormal_pdf, sizes_fit, numbers_fit,
            p0=[mu_init, sigma_init, amplitude_init],
            maxfev=10000
        )
        return popt
    except:
        return None


def calculate_d_values(df):
    """Calculate D10, D50, D90 from cumulative distribution."""
    if 'Number_cumulative' not in df.columns:
        return {}, {}
    
    sizes = df['Size / nm'].values
    cum = df['Number_cumulative'].values
    
    if cum[-1] > 0:
        cum_norm = cum / cum[-1]
    else:
        return {}, {}
    
    cum_std = np.zeros_like(cum)
    if 'Number_cumulative_std' in df.columns:
        cum_std = df['Number_cumulative_std'].values / cum[-1]
    
    d_values = {}
    
    for percentile, target in [(10, 0.1), (50, 0.5), (90, 0.9)]:
        idx = np.argmin(np.abs(cum_norm - target))
        d_values[f'D{percentile}'] = sizes[idx]
    
    return d_values, {}


def create_linear_number_plot(df, figsize=(14, 10)):
    """Create professional linear number-weighted distribution plot."""
    
    if 'Size / nm' not in df.columns or 'Number_normalized' not in df.columns:
        raise Exception("Required columns not found")
    
    d_values, _ = calculate_d_values(df)
    
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=figsize, sharex=True)
    
    sizes = df['Size / nm'].values
    numbers_norm = df['Number_normalized'].values
    numbers_std = df['Number_normalized_std'].values if 'Number_normalized_std' in df.columns else np.zeros_like(numbers_norm)
    
    # Subplot 1: Distribution with fit
    width = np.diff(sizes).mean()
    ax1.bar(sizes, numbers_norm, width=width, alpha=0.65, color='#0173B2', 
            edgecolor='#0173B2', linewidth=0.5, label='Number Distribution')
    
    ax1.fill_between(sizes, numbers_norm - numbers_std, numbers_norm + numbers_std, 
                     alpha=0.25, color='#0173B2', label='±1 SD')
    
    fit_popt = fit_lognormal(sizes, numbers_norm)
    if fit_popt is not None:
        sizes_smooth = np.linspace(sizes.min(), sizes.max(), 500)
        fit_curve = lognormal_pdf(sizes_smooth, *fit_popt)
        ax1.plot(sizes_smooth, fit_curve, color='#DE8F05', linewidth=2.5, 
                label='Lognormal Fit', zorder=5)
    
    ax1.set_ylabel('Normalized Number Distribution', fontsize=12, fontweight='bold')
    ax1.set_title('Linear Scale - Number-Weighted Particle Size Distribution', 
                 fontsize=13, fontweight='bold', pad=15)
    ax1.grid(True, alpha=0.3, linestyle='--', linewidth=0.5)
    ax1.legend(loc='upper right', fontsize=10, framealpha=0.95)
    ax1.set_xlim([sizes.min(), sizes.max()])
    
    # Subplot 2: Cumulative
    if 'Number_cumulative' in df.columns:
        cum = df['Number_cumulative'].values
        cum_max = cum[-1]
        
        if cum_max > 0:
            cum_norm = cum / cum_max
            cum_std = np.zeros_like(cum)
            
            if 'Number_cumulative_std' in df.columns:
                cum_std = df['Number_cumulative_std'].values / cum_max
            
            ax2.plot(sizes, cum_norm, color='#029E73', linewidth=2.5, label='Cumulative Distribution')
            ax2.fill_between(sizes, np.maximum(cum_norm - cum_std, 0), 
                           np.minimum(cum_norm + cum_std, 1), 
                           alpha=0.2, color='#029E73', label='±1 SD')
            
            colors = {'D10': '#CC78BC', 'D50': '#CA0020', 'D90': '#009ACD'}
            
            for d_name, d_val in d_values.items():
                if d_val > 0:
                    target = float(d_name[1:]) / 100.0
                    ax2.axvline(d_val, color=colors[d_name], linestyle='--', 
                               linewidth=1.5, alpha=0.7, zorder=3)
                    ax2.axhline(target, color=colors[d_name], linestyle=':', 
                               linewidth=1, alpha=0.4)
                    ax2.text(d_val, target + 0.05, f'{d_name}\n{d_val:.1f} nm', 
                            ha='center', fontsize=9, color=colors[d_name],
                            fontweight='bold', bbox=dict(boxstyle='round,pad=0.3', 
                            facecolor='white', edgecolor=colors[d_name], alpha=0.8))
        
        ax2.set_ylabel('Cumulative Distribution', fontsize=12, fontweight='bold')
        ax2.set_xlabel('Size (nm)', fontsize=12, fontweight='bold')
        ax2.set_ylim([0, 1.05])
        ax2.grid(True, alpha=0.3, linestyle='--', linewidth=0.5)
        ax2.legend(loc='lower right', fontsize=10, framealpha=0.95)
    
    for ax in [ax1, ax2]:
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.tick_params(labelsize=10)
    
    plt.tight_layout()
    return fig, d_values


# ============================================================================
# PROCESSING FUNCTIONS
# ============================================================================

def apply_dilution_correction(df, dilution_factor=1.0):
    """Apply dilution factor to concentration columns."""
    df = df.copy()
    
    for col in df.columns:
        if 'Concentration' in col and '_std' not in col:
            df[col] = df[col] * dilution_factor
            if f'{col}_std' in df.columns:
                df[f'{col}_std'] = df[f'{col}_std'] * dilution_factor
    
    return df


def normalize_distribution(df):
    """Normalize to unit area."""
    df = df.copy()
    
    if 'Size / nm' in df.columns and 'Number' in df.columns:
        sizes = df['Size / nm'].values
        numbers = df['Number'].values
        
        integral = trapezoid(numbers, sizes)
        
        if integral > 0:
            df['Number_normalized'] = numbers / integral
            
            if 'Number_std' in df.columns:
                df['Number_normalized_std'] = df['Number_std'] / integral
    
    return df


def calculate_cumulative(df):
    """Calculate cumulative distribution."""
    df = df.copy()
    df = df.sort_values('Size / nm').reset_index(drop=True)
    
    if 'Number_normalized' in df.columns:
        df['Number_cumulative'] = df['Number_normalized'].cumsum()
        
        if 'Number_normalized_std' in df.columns:
            df['Number_cumulative_std'] = np.sqrt((df['Number_normalized_std'] ** 2).cumsum())
    
    return df


# ============================================================================
# MAIN ANALYZER CLASS
# ============================================================================

class NTAAnalyzer:
    """Main NTA analysis class using your original code."""
    
    def __init__(self, config=None):
        self.config = config or DEFAULT_CONFIG
        self.results = {}
    
    def process(self, filepaths, metadata_overrides=None):
        """Process NTA files."""
        
        if isinstance(filepaths, str):
            filepaths = [filepaths]
        
        # Read and process all files
        distributions = []
        metadatas = []
        
        for filepath in filepaths:
            success, content = read_nta_file(filepath)
            if not success:
                print(f"⚠ Skipped {filepath}: {content}")
                continue
            
            success, sections = identify_sections(content)
            if not success:
                print(f"⚠ Skipped {filepath}: {sections}")
                continue
            
            success, msg = validate_file_structure(content, sections)
            if not success:
                print(f"⚠ Skipped {filepath}: {msg}")
                continue
            
            filename = os.path.basename(filepath)
            
            # Extract distribution
            dist = extract_single_file_distribution(content, sections, filename)
            if dist is None:
                print(f"⚠ Could not extract distribution from {filepath}")
                continue
            
            distributions.append(dist)
            
            # Extract metadata
            metadata = extract_all_metadata_fields(content, filename)
            metadatas.append(metadata)
        
        if not distributions:
            raise Exception("No files could be processed successfully")
        
        # Average distributions
        avg_dist = average_replicate_data(distributions, [m['filename'] for m in metadatas])
        
        # Use first metadata
        metadata = metadatas[0].copy()
        if metadata_overrides:
            metadata.update(metadata_overrides)
        
        metadata['num_files'] = len(distributions)
        
        # Apply corrections
        avg_dist = apply_dilution_correction(avg_dist, metadata.get('dilution_factor', 1.0))
        avg_dist = normalize_distribution(avg_dist)
        avg_dist = calculate_cumulative(avg_dist)
        
        self.results = {
            'distribution': avg_dist,
            'metadata': metadata,
        }
        
        return self.results
    
    def get_plot(self):
        """Generate the plot."""
        if 'distribution' not in self.results:
            raise Exception("No results. Run process() first.")
        
        fig, d_values = create_linear_number_plot(self.results['distribution'])
        self.results['d_values'] = d_values
        return fig
    
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
        
        # Save stats
        stats_path = os.path.join(output_dir, f'Stats_{unique_id}.txt')
        with open(stats_path, 'w') as f:
            f.write("=== PARTICLE SIZE DISTRIBUTION STATISTICS ===\n\n")
            
            if 'd_values' in self.results:
                f.write("D-VALUE STATISTICS (nm):\n")
                for d_name in ['D10', 'D50', 'D90']:
                    if d_name in self.results['d_values']:
                        d_val = self.results['d_values'][d_name]
                        f.write(f"  {d_name}: {d_val:.2f}\n")
                f.write("\n")
            
            f.write("DISTRIBUTION SUMMARY:\n")
            dist = self.results['distribution']
            f.write(f"  Size range: {dist['Size / nm'].min():.2f} - {dist['Size / nm'].max():.2f} nm\n")
            f.write(f"  Number of files: {self.results['metadata'].get('num_files', 1)}\n")
            f.write(f"  Dilution factor: {self.results['metadata'].get('dilution_factor', 1):.1f}\n")
        
        created_files.append(stats_path)
        
        return created_files
