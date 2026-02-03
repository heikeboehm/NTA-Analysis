"""
NTA Data Analysis - Complete Integrated Script

Fully integrates all functions from your original notebook:
- Cell 02: File I/O, section identification, validation (6 functions)
- Cell 03: Data extraction, averaging, replication handling (4 functions)
- Cell 04: Comprehensive metadata extraction & analysis (6 functions)
- Plus: Professional plotting and complete workflow

This is your actual production code from the notebook, adapted for standalone use.
"""

import os
import re
import json
import ntpath
from datetime import date

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from scipy.optimize import curve_fit
from scipy.integrate import trapezoid

# ============================================================================
# CELL 02 - FILE I/O MODULE (6 functions)
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


def process_single_file_content(filepath):
    """
    Complete file processing workflow for a single file.
    
    Parameters:
    filepath (str): Path to file
    
    Returns:
    tuple: (success_flag, result)
        - If successful: (True, (content, sections))
        - If failed: (False, error_message)
    """
    success, content = read_nta_file(filepath)
    if not success:
        return False, content
    
    success, sections = identify_sections(content)
    if not success:
        return False, sections
    
    success, msg = validate_file_structure(content, sections)
    if not success:
        return False, msg
    
    return True, (content, sections)


def process_multiple_file_contents(filepaths):
    """
    Process multiple NTA files, skipping any that fail with error messages.
    
    Parameters:
    filepaths (list): List of file paths
    
    Returns:
    tuple: (successful_files, failed_files)
        - successful_files: List of (filepath, content, sections) tuples
        - failed_files: List of (filepath, error_message) tuples
    """
    successful_files = []
    failed_files = []
    
    for filepath in filepaths:
        success, result = process_single_file_content(filepath)
        
        if success:
            content, sections = result
            successful_files.append((filepath, content, sections))
        else:
            failed_files.append((filepath, result))
    
    return successful_files, failed_files


def preview_file_content(content, max_chars=500):
    """
    Display a preview of the file content to verify it was read correctly.
    
    Parameters:
    content (str): File content
    max_chars (int): Maximum characters to preview
    
    Returns:
    str: Preview string
    """
    if len(content) > max_chars:
        preview = content[:max_chars] + f"\n... [truncated, {len(content)} total chars]"
    else:
        preview = content
    
    return preview


# ============================================================================
# CELL 03 - DATA EXTRACTION & AVERAGING (4 functions)
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
    
    Extracts both linear and logarithmic scale data with proper section handling.
    
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
# CELL 04 - METADATA EXTRACTION & ANALYSIS (6 functions)
# ============================================================================

def extract_all_metadata_fields(content, filename):
    """
    Extract ALL possible metadata fields from file content using comprehensive regex patterns.
    
    Parameters:
    content (str): File content to extract metadata from
    filename (str): Original filename for reference
    
    Returns:
    dict: Dictionary with all found metadata fields
    """
    # Comprehensive regex patterns for ALL possible NTA metadata fields
    metadata_patterns = [
        ('original_file', r'Original File:\s+(.+?)(?:\s+Section:|$)'),
        ('section', r'Section:\s+(.+)'),
        ('operator', r'Operator:\s+(.+)'),
        ('experiment', r'Experiment:\s+(.+)'),
        ('zetaview_sn', r'ZetaView S/N:\s+(.+)'),
        ('cell_sn', r'Cell S/N:\s+(.+)'),
        ('software', r'Software:\s+(.+?)(?:\s+Analyze:|$)'),
        ('analyze', r'Analyze:\s+(.+)'),
        ('sop', r'SOP:\s+(.+)'),
        ('sample', r'Sample:\s+(.+)'),
        ('electrolyte', r'Electrolyte:(?:\s*(.*?))?(?:\r?\n|$)'),
        ('ph', r'pH:\s+(.+?)(?:\s+entered|$)'),
        ('conductivity', r'Conductivity:\s+(.+?)(?:\s+sensed|$)'),
        ('temp_control', r'TempControl:\s+(.+)'),
        ('set_temperature', r'SetTemperature:\s+(.+)'),
        ('temperature', r'Temperature:\s+(.+?)(?:\s+sensed|$)'),
        ('viscosity', r'Viscosity:\s+(.+)'),
        ('date', r'Date:\s+(.+)'),
        ('time', r'Time:\s+(.+)'),
        ('general_remarks', r'General Remarks:\s+(.+)'),
        ('remarks', r'Remarks:\s+(.+)'),
        ('sample_info_1', r'Sample Info 1:\s+(.+)'),
        ('sample_info_2', r'Sample Info 2:\s+(.+)'),
        ('sample_info_3', r'Sample Info 3:\s*(.*)'),
        ('scattering_intensity', r'Scattering Intensity:\s+(.+)'),
        ('detected_particles', r'Detected Particles:\s+(.+)'),
        ('particle_drift_checked', r'Particle Drift Checked:\s+(.+)'),
        ('particle_drift_check_result', r'Particle Drift Check Result:\s+(.+)'),
        ('cell_check_date', r'Cell Checked:\s+(\d{4}-\d{2}-\d{2})'),
        ('cell_check_result', r'Cell Check Result:\s+(.+)'),
        ('type_of_measurement', r'Type of Measurement:\s+(.+)'),
        ('positions', r'Positions:\s+(.+)'),
        ('microscope_position', r'Microscope Position:\s+(.+)'),
        ('number_of_traces', r'Number of Traces:\s+(\d+)'),
        ('average_number_of_particles', r'Average Number of Particles:\s+(\d+\.\d+)'),
        ('dilution', r'Dilution::\s+(\d+\.\d+)'),
        ('concentration_correction_factor', r'Concentration Correction Factor:\s+(.+)'),
        ('laser_wavelength', r'Laser Wavelength nm:\s+(\d+\.\d+)'),
        ('median_number_d50', r'Median Number \(D50\):\s+(.+)'),
        ('median_concentration_d50', r'Median Concentration \(D50\):\s+(.+)'),
        ('median_volume_d50', r'Median Volume \(D50\):\s+(.+)'),
        ('minimum_brightness', r'Minimum Brightness:\s+(\d+)'),
        ('minimum_area', r'Minimum Area:\s+(\d+)'),
        ('maximum_area', r'Maximum Area:\s+(\d+)'),
        ('maximum_brightness', r'Maximum Brightness:\s+(\d+)'),
        ('tracking_radius2', r'Tracking Radius2:\s+(\d+)'),
        ('minimum_tracelength', r'Minimum Tracelength:\s+(\d+)'),
        ('fps', r'Camera:\s*FpSec\s+(\d+)\s+#Cycles'),
        ('cycles', r'#Cycles\s+(\d+)'),
        ('camera_settings', r'Camera:\s+(.+)'),
        ('frame_rate', r'FRate\s+(\d+\.\d+)'),
        ('auto_settings', r'Auto:\s+(.+)'),
        ('sensitivity', r'Sensitivity:\s+(.+)'),
        ('shutter', r'Shutter:\s+(.+)'),
        ('gain', r'Gain:\s+(.+)'),
    ]
    
    # Extract metadata
    metadata = {}
    
    for key, pattern in metadata_patterns:
        match = re.search(pattern, content)
        if match:
            value = match.group(1).strip() if match.group(1) else ''
            # Clean up common issues
            if value and not value.lower() in ['none', 'null', '']:
                metadata[key] = value
    
    # Add filename and derived fields
    metadata['filename'] = filename
    
    # Extract unique ID from filename
    base_name = os.path.splitext(filename)[0]
    if base_name.endswith("_rawdata"):
        base_name = base_name[:-8]
    if base_name.startswith("Data_"):
        base_name = base_name[5:]
    metadata['uniqueID'] = base_name
    
    # Try to find AVI file information
    original_file = metadata.get('original_file', '')
    if original_file:
        avi_filename = ntpath.basename(original_file)
        metadata['avi_filename'] = avi_filename
    
    return metadata


def extract_metadata_from_all_files(files_data):
    """
    Extract metadata from all files and organize by filename.
    
    Parameters:
    files_data (list): List of (filepath, content, sections) tuples
    
    Returns:
    tuple: (success_flag, all_files_metadata_dict)
    """
    if not files_data:
        return False, "No files provided for metadata extraction"
    
    all_files_metadata = {}
    
    for filepath, content, sections in files_data:
        filename = os.path.basename(filepath)
        metadata = extract_all_metadata_fields(content, filename)
        all_files_metadata[filename] = metadata
    
    return True, all_files_metadata


def analyze_field_differences(all_files_metadata):
    """
    Analyze which fields are identical vs. different across files.
    
    Parameters:
    all_files_metadata (dict): Dictionary of {filename: metadata_dict}
    
    Returns:
    tuple: (identical_fields, different_fields, field_analysis)
    """
    if not all_files_metadata:
        return {}, {}, {}
    
    # Get all unique field names across all files
    all_field_names = set()
    for metadata in all_files_metadata.values():
        all_field_names.update(metadata.keys())
    
    identical_fields = {}
    different_fields = {}
    field_analysis = {}
    
    for field_name in sorted(all_field_names):
        # Collect values for this field from all files
        values = []
        files_with_field = []
        
        for filename, metadata in all_files_metadata.items():
            if field_name in metadata:
                values.append(metadata[field_name])
                files_with_field.append(filename)
        
        # Analyze this field
        field_info = {
            'values': values,
            'files_with_field': files_with_field,
            'present_in_files': len(files_with_field),
            'total_files': len(all_files_metadata)
        }
        
        if len(set(values)) == 1:
            # All values are identical
            identical_fields[field_name] = values[0]
            field_info['status'] = 'identical'
        else:
            # Values differ
            different_fields[field_name] = values
            field_info['status'] = 'different'
        
        field_analysis[field_name] = field_info
    
    return identical_fields, different_fields, field_analysis


def smart_format_field(field_name, values):
    """
    Apply smart formatting rules based on field type and content.
    
    Parameters:
    field_name (str): Name of the metadata field
    values (list): List of values from different files
    
    Returns:
    tuple: (formatted_value, notes)
    """
    # File-specific info - keep as arrays
    file_specific_fields = [
        'filename', 'avi_filename', 'experiment', 'original_file', 'uniqueID', 
        'time', 'particle_drift_checked'
    ]
    
    # Text fields - use first value only
    text_first_fields = [
        'sample_info_1', 'sample_info_2', 'sample_info_3', 'remarks', 'general_remarks'
    ]
    
    # Sum fields - add them up
    sum_fields = [
        'number_of_traces', 'detected_particles'
    ]
    
    # Instrument-determined fields - keep as arrays
    instrument_determined_fields = [
        'median_number_d50', 'median_concentration_d50', 'median_volume_d50'
    ]
    
    notes = ""
    
    if field_name in file_specific_fields:
        return json.dumps(values), "file_specific"
    
    elif field_name in text_first_fields:
        return values[0], f"using_first_of_{len(values)}"
    
    elif field_name in sum_fields:
        try:
            numeric_values = [float(v) for v in values]
            total = sum(numeric_values)
            return f"{total:.0f}", f"sum_of_{len(values)}_measurements"
        except ValueError:
            return json.dumps(values), "non_numeric_sum_field"
    
    elif field_name in instrument_determined_fields:
        return json.dumps(values), "instrument_determined_per_measurement"
    
    else:
        # Default: use first value if all same, otherwise keep array
        unique_values = list(set(values))
        if len(unique_values) == 1:
            return values[0], "consistent"
        else:
            return json.dumps(values), f"varies_{len(unique_values)}_values"


def create_automated_metadata(all_files_metadata, identical_fields, different_fields, config=None):
    """
    Create standardized metadata for multi-file analysis.
    
    Combines identical fields and smartly formats different fields.
    
    Parameters:
    all_files_metadata (dict): Metadata from all files
    identical_fields (dict): Fields that are identical across files
    different_fields (dict): Fields that differ across files
    config (dict): Optional configuration
    
    Returns:
    dict: Standardized metadata
    """
    metadata = identical_fields.copy()
    
    # Add info about multi-file analysis
    metadata['num_files'] = len(all_files_metadata)
    metadata['extraction_date'] = str(date.today())
    metadata['file_list'] = json.dumps(list(all_files_metadata.keys()))
    
    # Smart format different fields
    for field_name, values in different_fields.items():
        formatted_value, notes = smart_format_field(field_name, values)
        metadata[field_name] = formatted_value
        if notes:
            metadata[f'{field_name}_note'] = notes
    
    # If we have a persistent ID, use the first one
    if 'uniqueID' not in metadata and all_files_metadata:
        first_file_metadata = list(all_files_metadata.values())[0]
        metadata['persistentID'] = first_file_metadata.get('uniqueID', 'multi_file_analysis')
    
    return metadata


def save_metadata_file(metadata, output_dir=None, config=None):
    """
    Save metadata to a file.
    
    Parameters:
    metadata (dict): Metadata to save
    output_dir (str): Output directory
    config (dict): Configuration
    
    Returns:
    tuple: (success, filepath_or_error)
    """
    if output_dir is None:
        output_dir = os.getcwd()
    
    os.makedirs(output_dir, exist_ok=True)
    
    unique_id = metadata.get('persistentID', 'metadata')
    filepath = os.path.join(output_dir, f'Data_{unique_id}_metadata.txt')
    
    try:
        with open(filepath, 'w') as f:
            for key, value in sorted(metadata.items()):
                f.write(f"{key}: {value}\n")
        return True, filepath
    except Exception as e:
        return False, str(e)


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
    """
    Main NTA analysis class integrating all cells 02-04.
    
    Full workflow:
    1. Read multiple files (Cell 02)
    2. Extract and average distributions (Cell 03)
    3. Extract and analyze metadata from all files (Cell 04)
    4. Apply corrections, normalize, calculate statistics
    5. Generate plots
    """
    
    def __init__(self, config=None):
        self.config = config or {}
        self.results = {}
    
    def process(self, filepaths, metadata_overrides=None):
        """Process NTA files with full workflow."""
        
        if isinstance(filepaths, str):
            filepaths = [filepaths]
        
        # ===== CELL 02: File I/O =====
        successful_files, failed_files = process_multiple_file_contents(filepaths)
        
        if not successful_files:
            raise Exception("No files could be processed successfully")
        
        # ===== CELL 03: Data Extraction =====
        distributions = []
        for filepath, content, sections in successful_files:
            filename = os.path.basename(filepath)
            dist = extract_single_file_distribution(content, sections, filename)
            if dist is not None:
                distributions.append(dist)
        
        if not distributions:
            raise Exception("Could not extract distribution data from any files")
        
        # Average distributions
        filenames = [os.path.basename(f[0]) for f in successful_files]
        avg_dist = average_replicate_data(distributions, filenames)
        
        # ===== CELL 04: Metadata Analysis =====
        success, all_metadata = extract_metadata_from_all_files(successful_files)
        
        # Analyze differences
        identical, different, analysis = analyze_field_differences(all_metadata)
        
        # Create standardized metadata
        metadata = create_automated_metadata(all_metadata, identical, different)
        
        # Apply overrides
        if metadata_overrides:
            metadata.update(metadata_overrides)
        
        metadata['num_files'] = len(successful_files)
        
        # Apply corrections
        avg_dist = apply_dilution_correction(avg_dist, metadata.get('dilution_factor', 1.0))
        avg_dist = normalize_distribution(avg_dist)
        avg_dist = calculate_cumulative(avg_dist)
        
        self.results = {
            'distribution': avg_dist,
            'metadata': metadata,
            'identical_fields': identical,
            'different_fields': different,
            'field_analysis': analysis,
            'all_file_metadata': all_metadata,
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
        
        # Save standardized metadata
        meta_path = os.path.join(output_dir, f'Data_{unique_id}_metadata.txt')
        with open(meta_path, 'w') as f:
            for key, value in sorted(self.results['metadata'].items()):
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
            
            # File-specific metadata summary
            if self.results.get('identical_fields'):
                f.write(f"\n=== IDENTICAL METADATA ACROSS FILES ===\n")
                f.write(f"Fields identical in all files: {len(self.results['identical_fields'])}\n")
            
            if self.results.get('different_fields'):
                f.write(f"\nFIELDS WITH DIFFERENCES:\n")
                for field, values in sorted(self.results['different_fields'].items()):
                    f.write(f"  {field}: {values}\n")
        
        created_files.append(stats_path)
        
        return created_files
