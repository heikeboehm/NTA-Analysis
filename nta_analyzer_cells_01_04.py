"""
NTA Data Analysis - Complete Integrated Script (Cells 01-06)

This script combines all functions from your original notebook:
- Cell 01: Configuration & utilities
- Cell 02: File I/O (read, identify sections, validate)
- Cell 03: Data extraction & averaging
- Cell 04: Metadata extraction & analysis
- Cell 05: Dilution correction, normalization, cumulative distributions, total metrics
- Cell 06: Statistics with bounds-based uncertainty (D-values, percentiles, span)

All code is extracted from your EXACT notebook, adapted for Streamlit.
"""

import os
import re
import json
import ntpath
from datetime import date

import pandas as pd
import numpy as np
from scipy import integrate


# ============================================================================
# CELL 01 - CONFIGURATION & UTILITIES
# ============================================================================

if 'CONFIG' not in globals():
    CONFIG = {
        "directory": ".",
        "file_identifier": ".txt",
        "output_subdirs": ["metadata", "processed"],
        "nta_concentration_calibration_factor": 4.61E+5,
        "project_metadata": {
            "experimenter": "Your_Initials",
            "location": "Your_Lab_Location",
            "project": "Your_Project_Name",
            "meta_version": "v03",
            "pi": "Principal_Investigator_Initials",
            "funding": "none",
            "data_collection_method": "NTA",
            "unit_of_analysis": '["nm", "nm^2", "nm^3"]',
            "keywords": "particle_size_distribution",
            "publications": "None",
        }
    }


def extract_sample_info(filename):
    """
    Extract sample information from filename for more readable display.
    Also extracts the persistentID (everything up to and including the date YYYYMMDD).
    (From Cell 01)
    """
    base_name = filename.replace("_rawdata.txt", "").replace("_rawdata", "")
    
    if base_name.startswith("Data_"):
        base_name = base_name[5:]
    
    # Extract persistentID: everything up to and including YYYYMMDD date
    import re
    date_match = re.search(r'(.+?)_(\d{8})(?:_|$)', base_name)
    persistent_id = ""
    date_str = ""
    
    if date_match:
        persistent_id = date_match.group(1) + "_" + date_match.group(2)  # e.g., CBP_LEAF_5455892_2_20250710
        date = date_match.group(2)
        date_str = f" (Date: {date[0:4]}-{date[4:6]}-{date[6:8]})"
    else:
        persistent_id = base_name
    
    return persistent_id, f"{persistent_id}{date_str}"


# ============================================================================
# CELL 02 - FILE I/O MODULE
# ============================================================================

def read_nta_file(filepath):
    """
    Read an NTA data file with appropriate encoding.
    (From Cell 02)
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
    (From Cell 02)
    """
    sections = {}
    
    size_lin_start = content.find("Size Distribution")
    if size_lin_start == -1:
        return False, "Could not find 'Size Distribution' section for linear data"
    sections['linear_start'] = size_lin_start
    
    size_log_start = content.find("-1.000E+0")
    if size_log_start == -1:
        second_header = content.find("Size / nm", size_lin_start + 100)
        if second_header == -1:
            return False, "Could not find logarithmic data section"
        sections['logarithmic_start'] = second_header
    else:
        sections['logarithmic_start'] = size_log_start
    
    if sections['linear_start'] >= sections['logarithmic_start']:
        return False, "Invalid file structure: linear section should come before logarithmic section"
    
    return True, sections


def validate_file_structure(content, sections):
    """
    Validate file structure to ensure it can be processed.
    (From Cell 02)
    """
    linear_section = content[sections['linear_start']:sections['logarithmic_start']]
    if not re.search(r'Size / nm\s+Number\s+Concentration', linear_section):
        return False, "Missing expected header pattern in linear section"
    
    log_section = content[sections['logarithmic_start']:]
    if not re.search(r'[\d.-]+E[\+\-]\d+\s+[\d.-]+E[\+\-]\d+', log_section):
        return False, "Could not find data rows in logarithmic section"
    
    return True, "File structure is valid"


def process_single_file_content(filepath):
    """
    Complete file processing workflow for a single file.
    (From Cell 02)
    """
    success, result = read_nta_file(filepath)
    if not success:
        return False, result
    content = result
    
    success, result = identify_sections(content)
    if not success:
        return False, result
    sections = result
    
    success, message = validate_file_structure(content, sections)
    if not success:
        return False, message
    
    return True, (content, sections)


def process_multiple_file_contents(filepaths):
    """
    Process multiple NTA files, skipping any that fail with error messages.
    (From Cell 02)
    """
    successful_files = []
    failed_files = {}
    
    for filepath in filepaths:
        filename = os.path.basename(filepath)
        
        success, result = process_single_file_content(filepath)
        
        if success:
            content, sections = result
            successful_files.append((filename, content, sections))
        else:
            failed_files[filename] = result
    
    return successful_files, failed_files


def preview_file_content(content, max_chars=500):
    """
    Display a preview of the file content to verify it was read correctly.
    (From Cell 02)
    """
    preview = content[:max_chars]
    return "-" * 50 + "\n" + preview + "\n" + "-" * 50


# ============================================================================
# CELL 03 - DATA EXTRACTION & AVERAGING
# ============================================================================

def extract_data_section(content, start_pos, end_pos=None, is_log_section=False):
    """
    Extract tabular data from a section of the file content.
    (From Cell 03)
    """
    section = content[start_pos:end_pos]
    
    if is_log_section:
        sep_pos = section.find("-1.000E+0")
        if sep_pos != -1:
            sep_line_end = section.find('\n', sep_pos)
            if sep_line_end != -1:
                section = section[sep_line_end + 1:]
            else:
                section = section[sep_pos + 20:]
        
        second_header = section.find("Size Distribution")
        if second_header != -1:
            section = section[second_header:]
    
    header_match = re.search(r'Size / nm\s+Number\s+Concentration.+', section)
    
    if not header_match and is_log_section:
        header_match = re.search(r'Size.*Number', section)
        if not header_match:
            header_line = "Size / nm\tNumber\tConcentration / cm-3\tVolume / nm^3\tArea / nm^2"
            data_lines = []
            for line in section.split('\n'):
                if "-1.000E+0" in line:
                    continue
                if re.match(r'^\s*[\d.-]+E[\+\-]\d+\s+[\d.-]+E[\+\-]\d+', line):
                    data_lines.append(line)
            
            if not data_lines:
                return False, "Could not find any valid data lines in logarithmic section"
            return True, (header_line, data_lines)
        else:
            header_line = header_match.group(0)
            data_start = section.find(header_line) + len(header_line)
            data_section = section[data_start:]
    else:
        if header_match:
            header_line = header_match.group(0)
            data_start = section.find(header_line) + len(header_line)
            data_section = section[data_start:]
        else:
            return False, "Could not find header line in data section"
    
    all_lines = data_section.split('\n')
    data_lines = []
    
    for line in all_lines:
        if "-1.000E+0" in line:
            continue
        if re.match(r'^\s*[\d.-]+E[\+\-]\d+\s+[\d.-]+E[\+\-]\d+', line):
            data_lines.append(line)
        elif len(data_lines) > 0 and (line.strip() == '' or line.strip().startswith('-')):
            break
    
    if not data_lines:
        return False, "No data lines found in section"
    
    return True, (header_line, data_lines)


def parse_data_lines(header_line, data_lines, scale_type):
    """
    Parse data lines into a structured DataFrame.
    (From Cell 03)
    """
    parsed_data = []
    
    for line in data_lines:
        if "-1.000E+0" in line:
            continue
            
        values = re.findall(r'[\d.-]+E[\+\-]\d+', line)
        
        if len(values) >= 5:
            try:
                if values[0] == "-1.000E+0":
                    continue
                
                row = {
                    'size_nm': float(values[0]),
                    'number': float(values[1]),
                    'concentration_cm-3': float(values[2]),
                    'volume_nm^3': float(values[3]),
                    'area_nm^2': float(values[4]),
                    'scale': scale_type
                }
                parsed_data.append(row)
            except ValueError:
                continue
    
    if not parsed_data:
        return False, f"Failed to parse any data lines for {scale_type} scale"
    
    return True, pd.DataFrame(parsed_data)


def extract_single_file_distribution(content, sections, filename):
    """
    Extract distribution data from a single file.
    (From Cell 03)
    """
    lin_start = sections['linear_start']
    log_start = sections['logarithmic_start']
    
    success, lin_result = extract_data_section(content, lin_start, log_start, is_log_section=False)
    if not success:
        return False, f"Failed to extract linear data from {filename}: {lin_result}"
    
    lin_header, lin_data_lines = lin_result
    success, lin_df = parse_data_lines(lin_header, lin_data_lines, 'linear')
    if not success:
        return False, f"Failed to parse linear data from {filename}: {lin_df}"
    
    success, log_result = extract_data_section(content, log_start, is_log_section=True)
    if not success:
        combined_df = lin_df.copy()
    else:
        log_header, log_data_lines = log_result
        success, log_df = parse_data_lines(log_header, log_data_lines, 'logarithmic')
        if not success:
            combined_df = lin_df.copy()
        else:
            combined_df = pd.concat([lin_df, log_df], ignore_index=True)
    
    combined_df['source_file'] = filename
    
    return True, combined_df


def average_replicate_data(dataframes_list, filenames_list):
    """
    Average distribution data across multiple replicates bin-by-bin.
    (From Cell 03)
    """
    if not dataframes_list:
        return False, "No dataframes to average"
    
    if len(dataframes_list) == 1:
        df = dataframes_list[0].copy()
        
        df['number_sd'] = 0.0
        df['concentration_cm-3_sd'] = 0.0
        df['volume_nm^3_sd'] = 0.0
        df['area_nm^2_sd'] = 0.0
        
        df.rename(columns={
            'number': 'number_avg',
            'concentration_cm-3': 'concentration_cm-3_avg',
            'volume_nm^3': 'volume_nm^3_avg',
            'area_nm^2': 'area_nm^2_avg'
        }, inplace=True)
        
        df['num_replicates'] = 1
        df['source_files'] = filenames_list[0]
        
        return True, df
    
    averaged_dfs = []
    
    for scale in ['linear', 'logarithmic']:
        scale_dfs = []
        for df in dataframes_list:
            scale_data = df[df['scale'] == scale]
            if not scale_data.empty:
                scale_dfs.append(scale_data)
        
        if not scale_dfs:
            continue
        
        all_sizes = []
        for df in scale_dfs:
            all_sizes.extend(df['size_nm'].values)
        unique_sizes = sorted(set(all_sizes))
        
        averaged_data = []
        
        for size in unique_sizes:
            numbers = []
            concentrations = []
            volumes = []
            areas = []
            
            for df in scale_dfs:
                size_row = df[df['size_nm'] == size]
                if not size_row.empty:
                    numbers.append(size_row['number'].iloc[0])
                    concentrations.append(size_row['concentration_cm-3'].iloc[0])
                    volumes.append(size_row['volume_nm^3'].iloc[0])
                    areas.append(size_row['area_nm^2'].iloc[0])
            
            if numbers:
                row = {
                    'size_nm': size,
                    'number_avg': np.mean(numbers),
                    'number_sd': np.std(numbers, ddof=1) if len(numbers) > 1 else 0.0,
                    'concentration_cm-3_avg': np.mean(concentrations),
                    'concentration_cm-3_sd': np.std(concentrations, ddof=1) if len(concentrations) > 1 else 0.0,
                    'volume_nm^3_avg': np.mean(volumes),
                    'volume_nm^3_sd': np.std(volumes, ddof=1) if len(volumes) > 1 else 0.0,
                    'area_nm^2_avg': np.mean(areas),
                    'area_nm^2_sd': np.std(areas, ddof=1) if len(areas) > 1 else 0.0,
                    'scale': scale,
                    'num_replicates': len(numbers),
                    'source_files': '; '.join(filenames_list)
                }
                averaged_data.append(row)
        
        if averaged_data:
            scale_df = pd.DataFrame(averaged_data)
            averaged_dfs.append(scale_df)
    
    if not averaged_dfs:
        return False, "No data could be averaged"
    
    final_df = pd.concat(averaged_dfs, ignore_index=True)
    final_df = final_df.sort_values(['scale', 'size_nm']).reset_index(drop=True)
    
    return True, final_df


# ============================================================================
# CELL 04 - METADATA EXTRACTION & ANALYSIS
# ============================================================================

def extract_all_metadata_fields(content, filename):
    """
    Extract ALL possible metadata fields from file content using comprehensive regex patterns.
    (From Cell 04)
    """
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
        # NOTE: Removed median_*_d50 fields - calculated separately in Cell 05/07
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
    
    metadata = {}
    
    for key, pattern in metadata_patterns:
        match = re.search(pattern, content)
        if match:
            value = match.group(1).strip() if match.group(1) else ''
            if value and not value.lower() in ['none', 'null', '']:
                metadata[key] = value
    
    metadata['filename'] = filename
    
    # Extract persistentID (up to and including YYYYMMDD date)
    persistent_id, _ = extract_sample_info(filename)
    metadata['uniqueID'] = persistent_id
    
    original_file = metadata.get('original_file', '')
    if original_file:
        avi_filename = ntpath.basename(original_file)
        metadata['avi_filename'] = avi_filename
    
    return metadata


def extract_metadata_from_all_files(files_data):
    """
    Extract metadata from all files and organize by filename.
    (From Cell 04)
    """
    if not files_data:
        return False, "No files provided for metadata extraction"
    
    all_files_metadata = {}
    
    for filename, content, sections in files_data:
        metadata = extract_all_metadata_fields(content, filename)
        all_files_metadata[filename] = metadata
    
    return True, all_files_metadata


def analyze_field_differences(all_files_metadata):
    """
    Analyze which fields are identical vs. different across files.
    (From Cell 04)
    """
    if not all_files_metadata:
        return {}, {}, {}
    
    all_field_names = set()
    for metadata in all_files_metadata.values():
        all_field_names.update(metadata.keys())
    
    identical_fields = {}
    different_fields = {}
    field_analysis = {}
    
    for field_name in sorted(all_field_names):
        values = []
        files_with_field = []
        
        for filename, metadata in all_files_metadata.items():
            if field_name in metadata:
                values.append(metadata[field_name])
                files_with_field.append(filename)
        
        field_info = {
            'values': values,
            'files_with_field': files_with_field,
            'present_in_files': len(files_with_field),
            'total_files': len(all_files_metadata)
        }
        
        if len(set(values)) == 1:
            identical_fields[field_name] = values[0]
            field_info['status'] = 'identical'
        else:
            different_fields[field_name] = values
            field_info['status'] = 'different'
        
        field_analysis[field_name] = field_info
    
    return identical_fields, different_fields, field_analysis


def smart_format_field(field_name, values):
    """
    Apply smart formatting rules based on field type and content.
    Special rules:
    - Temperature: Alert if difference > 1°C
    - Dilution: Alert if any difference
    - Others: Alert if CV > 10%
    (From Cell 04)
    """
    file_specific_fields = [
        'filename', 'avi_filename', 'experiment', 'original_file', 'uniqueID', 
        'time', 'particle_drift_checked'
    ]
    
    text_first_fields = [
        'sample_info_1', 'sample_info_2', 'sample_info_3', 'remarks', 'general_remarks'
    ]
    
    sum_fields = [
        'number_of_traces'
    ]
    
    quality_control_fields = [
        'particle_drift_check_result', 'cell_check_result'
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
    
    elif field_name in quality_control_fields:
        unique_values = list(set(values))
        
        # Special check for particle_drift_check_result
        if field_name == 'particle_drift_check_result':
            acceptable_values = ['Good', 'Very Good']
            bad_values = [v for v in values if v not in acceptable_values]
            
            if bad_values:
                # Has "too high" or "Poor" or other non-good values
                return json.dumps(values), f"QC_ALERT_particle_drift_not_good"
            elif len(unique_values) == 1:
                return values[0], "qc_consistent"
            else:
                return json.dumps(values), f"QC_ALERT_inconsistent_values"
        
        # General QC logic for other fields
        elif len(unique_values) == 1:
            return values[0], "qc_consistent"
        else:
            return json.dumps(values), f"QC_ALERT_inconsistent_values"
    
    elif field_name == 'temperature':
        # Special handling: Alert if temperature differs by more than 1°C
        try:
            numeric_values = [float(v) for v in values]
            mean_val = np.mean(numeric_values)
            std_val = np.std(numeric_values, ddof=1) if len(numeric_values) > 1 else 0.0
            max_range = max(numeric_values) - min(numeric_values)
            
            if max_range > 1.0:
                notes = f"TEMPERATURE_ALERT_range_{max_range:.2f}°C"
            else:
                notes = f"mean_sd_of_{len(values)}"
            
            return f"{mean_val:.2f} ± {std_val:.2f}", notes
        except ValueError:
            return json.dumps(values), "non_numeric"
    
    elif field_name == 'dilution':
        # Special handling: Alert if any different
        unique_values = list(set(values))
        if len(unique_values) > 1:
            return json.dumps(values), "DILUTION_ALERT_different_values"
        else:
            try:
                return f"{float(values[0]):.2f}", "identical"
            except ValueError:
                return values[0], "identical"
    
    else:
        # Try to calculate mean ± SD for numeric fields
        try:
            numeric_values = [float(v) for v in values]
            mean_val = np.mean(numeric_values)
            std_val = np.std(numeric_values, ddof=1) if len(numeric_values) > 1 else 0.0
            
            # Check for concerning variations
            cv = (std_val / mean_val * 100) if mean_val != 0 else 0
            
            if cv > 10:  # More than 10% coefficient of variation
                notes = f"HIGH_VARIATION_CV_{cv:.1f}%"
            else:
                notes = f"mean_sd_of_{len(values)}"
            
            return f"{mean_val:.2f} ± {std_val:.2f}", notes
            
        except ValueError:
            # Non-numeric field - keep as array
            return json.dumps(values), "non_numeric_different"


def create_automated_metadata(all_files_metadata, identical_fields, different_fields, config=None):
    """
    Create essential standardized metadata for multi-file analysis.
    (From Cell 04)
    """
    filenames = list(all_files_metadata.keys())
    num_files = len(filenames)
    
    metadata = {}
    processing_notes = {}
    
    if 'uniqueID' in identical_fields:
        base_id = identical_fields['uniqueID']
    else:
        base_id = list(all_files_metadata.values())[0].get('uniqueID', 'unknown')
    
    if num_files > 1:
        unique_id = f"{base_id}_avg{num_files}"
    else:
        unique_id = base_id
    
    project_meta = {}
    if config and "project_metadata" in config:
        project_meta = config["project_metadata"]
    
    # Check for manual persistent ID override
    if config and "manual_persistent_id" in config and config["manual_persistent_id"]:
        unique_id = config["manual_persistent_id"]
    
    metadata['experimenter'] = project_meta.get('experimenter', 'Your_Initials')
    metadata['location'] = project_meta.get('location', 'Your_Lab_Location')
    metadata['project'] = project_meta.get('project', 'Your_Project_Name')
    metadata['meta_version'] = project_meta.get('meta_version', 'v03')
    metadata['pi'] = project_meta.get('pi', 'Your_PI_Initials')
    metadata['funding'] = project_meta.get('funding', 'none')
    metadata['persistentID'] = unique_id
    metadata['data_collection_method'] = project_meta.get('data_collection_method', 'NTA')
    metadata['nta_instrument'] = 'ZetaView'
    
    if 'analyze' in identical_fields:
        metadata['nta_software'] = f"ZetaView {identical_fields['analyze']}"
    elif 'analyze' in different_fields:
        metadata['nta_software'] = f"ZetaView {different_fields['analyze'][0]}"
    else:
        metadata['nta_software'] = 'ZetaView'
    
    metadata['nta_processed_file'] = f"Data_{unique_id}_PSD.txt"
    
    if 'sample' in identical_fields:
        metadata['sample'] = identical_fields['sample']
    elif 'sample' in different_fields:
        metadata['sample'] = different_fields['sample'][0]
    
    metadata['num_replicates'] = num_files
    metadata['source_files'] = json.dumps(filenames)
    
    # Fields to exclude (instrument D50s - we calculate better ones)
    excluded_fields = {
        'median_number_d50', 'median_concentration_d50', 'median_volume_d50',
        'detected_particles'  # Use average_number_of_particles instead
    }
    
    essential_fields = {
        'date', 'temperature', 'ph', 'dilution', 'laser_wavelength', 'electrolyte',
        'positions', 'cycles', 'fps', 
        'particle_drift_check_result', 'cell_check_result',
        'average_number_of_particles', 'number_of_traces',
        'conductivity', 'scattering_intensity', 'viscosity',
        'avi_filesize'
    }
    
    # Remove excluded fields from different_fields
    different_fields = {k: v for k, v in different_fields.items() if k not in excluded_fields}
    
    for field_name, value in identical_fields.items():
        if field_name in essential_fields:
            if field_name in ['temperature', 'ph', 'dilution', 'laser_wavelength', 'positions', 
                             'cycles', 'fps', 'average_number_of_particles', 'number_of_traces', 
                             'particle_drift_check_result', 'cell_check_result',
                             'conductivity', 'scattering_intensity', 'viscosity', 'avi_filesize']:
                if field_name in ['number_of_traces']:
                    metadata[f'nta_{field_name}_sum'] = value
                else:
                    metadata[f'nta_{field_name}'] = value
            else:
                metadata[field_name] = value
    
    quality_alerts = []
    high_variation_fields = []
    
    # Check for particle drift issues even if identical across files
    if 'particle_drift_check_result' in identical_fields:
        drift_value = identical_fields['particle_drift_check_result']
        if drift_value not in ['Good', 'Very Good']:
            quality_alerts.append(f"particle_drift_check_result: {drift_value} (not 'Good' or 'Very Good')")
    
    for field_name, values in different_fields.items():
        if field_name in essential_fields:
            formatted_value, notes = smart_format_field(field_name, values)
            
            if "QC_ALERT" in notes:
                quality_alerts.append(f"{field_name}: {formatted_value}")
            if "HIGH_VARIATION" in notes:
                high_variation_fields.append(f"{field_name}: {notes}")
            
            if field_name in ['temperature', 'ph', 'dilution', 'laser_wavelength', 'positions', 
                             'cycles', 'fps', 'average_number_of_particles', 'number_of_traces', 
                             'particle_drift_check_result', 'cell_check_result',
                             'conductivity', 'scattering_intensity', 'viscosity', 'avi_filesize']:
                if field_name in ['number_of_traces']:
                    metadata[f'nta_{field_name}_sum'] = formatted_value
                else:
                    metadata[f'nta_{field_name}'] = formatted_value
            else:
                metadata[field_name] = formatted_value
            
            processing_notes[field_name] = notes
    
    metadata['python_analysis'] = str(date.today())
    
    if quality_alerts:
        metadata['quality_control_alerts'] = json.dumps(quality_alerts)
    
    if high_variation_fields:
        metadata['high_variation_fields'] = json.dumps(high_variation_fields)
    
    return metadata, quality_alerts, high_variation_fields


def save_metadata_file(metadata, output_dir=None, config=None, include_headers=True):
    """
    Save metadata to a file in the specified output directory.
    Can be organized into logical sections with clear headers, or as plain tab-delimited.
    (From Cell 04)
    
    Parameters:
    metadata (dict): Metadata dictionary to save
    output_dir (str): Output directory path (optional)
    config (dict): Configuration dictionary (optional)
    include_headers (bool): If True, organize with section headers. If False, plain format for Excel import.
    """
    if output_dir is None:
        if config is not None and "directory" in config:
            base_dir = config["directory"]
            output_dir = os.path.join(base_dir, "metadata")
        else:
            output_dir = os.path.join(os.getcwd(), "metadata")
    
    try:
        os.makedirs(output_dir, exist_ok=True)
    except Exception as e:
        return False, f"Failed to create metadata directory: {str(e)}"
    
    unique_id = metadata.get('persistentID', 'unknown')
    
    # Use different filenames based on include_headers
    if include_headers:
        metadata_path = os.path.join(output_dir, f"Data_{unique_id}_metadata.txt")
    else:
        metadata_path = os.path.join(output_dir, f"Data_{unique_id}_metadata_excel.txt")
    
    try:
        if include_headers:
            # Define field order and sections (with headers)
            sections = {
                'SAMPLE & PROJECT INFORMATION': [
                    'persistentID', 'sample', 'electrolyte', 'date', 'num_replicates'
                ],
                'EXPERIMENTER & LOCATION': [
                    'experimenter', 'location', 'project', 'pi', 'funding'
                ],
                'INSTRUMENT & METHOD': [
                    'data_collection_method', 'nta_instrument', 'nta_software'
                ],
                'MEASUREMENT CONDITIONS': [
                    'nta_temperature', 'nta_ph', 'nta_conductivity', 
                    'nta_dilution', 'nta_viscosity'
                ],
                'MEASUREMENT PARAMETERS': [
                    'nta_laser_wavelength', 'nta_positions', 'nta_cycles', 'nta_fps'
                ],
                'DATA QUALITY & STATISTICS': [
                    'nta_average_number_of_particles', 'nta_number_of_traces_sum', 
                    'nta_scattering_intensity'
                ],
                'QUALITY CONTROL': [
                    'nta_particle_drift_check_result', 'nta_cell_check_result'
                ],
                'FILE REFERENCES': [
                    'nta_processed_file', 'source_files', 'meta_version'
                ],
                'CALCULATIONS': [
                    'nta_python_analysis', 'nta_metrics_scale',
                    'nta_specific_surface_area_m^2_per_cm^3',
                    'nta_total_particles_per_mL',
                    'nta_total_volume_uL_per_mL',
                    'nta_volume_percentage',
                    'nta_linear_number_d10',
                    'nta_linear_number_d50',
                    'nta_linear_number_d90',
                    'nta_linear_number_span',
                    'nta_linear_volume_d10',
                    'nta_linear_volume_d50',
                    'nta_linear_volume_d90',
                    'nta_linear_volume_span',
                    'nta_linear_surface_area_d10',
                    'nta_linear_surface_area_d50',
                    'nta_linear_surface_area_d90',
                    'nta_linear_surface_area_span'
                ],
                'QUALITY ALERTS': [
                    'quality_control_alerts', 'high_variation_fields'
                ]
            }
            
            # Write metadata with section headers
            with open(metadata_path, 'w') as f:
                for section, field_names in sections.items():
                    # Write section header
                    f.write(f"\n# ============ {section} ============\n")
                    
                    # Write fields in this section
                    for field_name in field_names:
                        if field_name in metadata:
                            value = metadata[field_name]
                            f.write(f"{field_name}\t{value}\t\n")
                
                # Write any remaining fields not in sections
                used_fields = set()
                for field_names in sections.values():
                    used_fields.update(field_names)
                
                remaining = set(metadata.keys()) - used_fields
                if remaining:
                    f.write(f"\n# ============ OTHER ============\n")
                    for field_name in sorted(remaining):
                        value = metadata[field_name]
                        f.write(f"{field_name}\t{value}\t\n")
        else:
            # Plain tab-delimited format (no headers - for Excel import)
            with open(metadata_path, 'w') as f:
                for key, value in sorted(metadata.items()):
                    f.write(f"{key}\t{value}\n")
        
        return True, metadata_path
    except Exception as e:
        return False, f"Failed to write metadata file: {str(e)}"


# ============================================================================
# CELL 05: DILUTION CORRECTION & NORMALIZATION WITH UNCERTAINTY PROPAGATION
# ============================================================================

def apply_dilution_correction_with_uncertainty(df, metadata=None, manual_dilution=None):
    """
    Apply dilution correction to all measured values with uncertainty propagation.
    For a dilution factor D: actual sample concentration = measured × D
    Uncertainty propagation: σ_actual = σ_measured × D
    """
    updated_df = df.copy()
    
    dilution_factor = 1.0
    dilution_source = "default (no dilution)"
    
    if manual_dilution is not None:
        try:
            dilution_factor = float(manual_dilution)
            dilution_source = "manually specified"
        except (ValueError, TypeError):
            return False, f"Invalid manual dilution factor: {manual_dilution}"
    elif metadata is not None and 'nta_dilution' in metadata:
        try:
            dilution_string = metadata['nta_dilution']
            dilution_factor = float(dilution_string.split('±')[0].strip()) if '±' in dilution_string else float(dilution_string)
            dilution_source = "metadata (nta_dilution)"
        except (ValueError, TypeError):
            print("Warning: Could not parse dilution factor from metadata, using default (1.0)")
    
    print(f"Applying dilution correction: factor = {dilution_factor} (source: {dilution_source})")
    
    # Apply dilution correction to concentration
    if 'concentration_cm-3_avg' in updated_df.columns:
        updated_df['particles_per_mL_avg'] = updated_df['concentration_cm-3_avg'] * dilution_factor
        if 'concentration_cm-3_sd' in updated_df.columns:
            updated_df['particles_per_mL_sd'] = updated_df['concentration_cm-3_sd'] * dilution_factor
        updated_df = updated_df.drop(['concentration_cm-3_avg', 'concentration_cm-3_sd'], axis=1, errors='ignore')
    else:
        return False, "Missing concentration_cm-3_avg column for dilution correction"
    
    # Apply dilution correction to volume
    if 'volume_nm^3_avg' in updated_df.columns:
        updated_df['volume_nm^3_per_mL_avg'] = updated_df['volume_nm^3_avg'] * dilution_factor
        if 'volume_nm^3_sd' in updated_df.columns:
            updated_df['volume_nm^3_per_mL_sd'] = updated_df['volume_nm^3_sd'] * dilution_factor
        updated_df = updated_df.drop(['volume_nm^3_avg', 'volume_nm^3_sd'], axis=1, errors='ignore')
    
    # Apply dilution correction to area
    if 'area_nm^2_avg' in updated_df.columns:
        updated_df['area_nm^2_per_mL_avg'] = updated_df['area_nm^2_avg'] * dilution_factor
        if 'area_nm^2_sd' in updated_df.columns:
            updated_df['area_nm^2_per_mL_sd'] = updated_df['area_nm^2_sd'] * dilution_factor
        updated_df = updated_df.drop(['area_nm^2_avg', 'area_nm^2_sd'], axis=1, errors='ignore')
    
    return True, updated_df


def normalize_distributions_with_uncertainty(df, size_column='size_nm'):
    """
    Normalize particle distributions by area under the curve with uncertainty propagation.
    Creates normalized number distributions from the averaged number data.
    """
    normalized_df = df.copy()
    
    for scale in normalized_df['scale'].unique():
        scale_mask = normalized_df['scale'] == scale
        scale_data = normalized_df[scale_mask].copy()
        
        if scale_data.empty or 'number_avg' not in scale_data.columns:
            continue
        
        scale_data = scale_data.sort_values(size_column)
        sizes = scale_data[size_column].values
        numbers_avg = scale_data['number_avg'].values
        
        if len(sizes) < 2:
            continue
        
        area_avg = np.trapz(numbers_avg, sizes)
        
        if area_avg > 0:
            normalized_df.loc[scale_mask, 'number_normalized_avg'] = numbers_avg / area_avg
            
            if 'number_sd' in scale_data.columns:
                numbers_sd = scale_data['number_sd'].values
                normalized_df.loc[scale_mask, 'number_normalized_sd'] = numbers_sd / area_avg
            else:
                normalized_df.loc[scale_mask, 'number_normalized_sd'] = 0.0
        else:
            normalized_df.loc[scale_mask, 'number_normalized_avg'] = 0.0
            normalized_df.loc[scale_mask, 'number_normalized_sd'] = 0.0
    
    return normalized_df


def calculate_cumulative_distributions_with_uncertainty(df, scale_column='scale'):
    """
    Calculate cumulative distributions with proper uncertainty propagation.
    For independent uncertainties: σ_cumsum[j] = √(Σ(i=0 to j) σ[i]²)
    """
    result_df = df.copy()
    
    for scale in result_df[scale_column].unique():
        scale_mask = result_df[scale_column] == scale
        scale_indices = result_df[scale_mask].sort_values('size_nm').index
        
        if len(scale_indices) == 0:
            continue
        
        # Normalized number distribution
        if 'number_normalized_avg' in result_df.columns:
            cumsum_avg = result_df.loc[scale_indices, 'number_normalized_avg'].cumsum()
            if cumsum_avg.iloc[-1] > 0:
                result_df.loc[scale_indices, 'number_normalized_cumsum_avg'] = cumsum_avg / cumsum_avg.iloc[-1]
            else:
                result_df.loc[scale_indices, 'number_normalized_cumsum_avg'] = 0
            
            if 'number_normalized_sd' in result_df.columns:
                normalized_var_cumsum = (result_df.loc[scale_indices, 'number_normalized_sd'] ** 2).cumsum()
                cumsum_sd = np.sqrt(normalized_var_cumsum)
                if cumsum_avg.iloc[-1] > 0:
                    result_df.loc[scale_indices, 'number_normalized_cumsum_sd'] = cumsum_sd / cumsum_avg.iloc[-1]
                else:
                    result_df.loc[scale_indices, 'number_normalized_cumsum_sd'] = 0
        
        # Absolute volume distribution
        if 'volume_nm^3_per_mL_avg' in result_df.columns:
            cumsum_avg = result_df.loc[scale_indices, 'volume_nm^3_per_mL_avg'].cumsum()
            result_df.loc[scale_indices, 'volume_nm^3_per_mL_cumsum_avg'] = cumsum_avg
            
            if 'volume_nm^3_per_mL_sd' in result_df.columns:
                volume_var_cumsum = (result_df.loc[scale_indices, 'volume_nm^3_per_mL_sd'] ** 2).cumsum()
                result_df.loc[scale_indices, 'volume_nm^3_per_mL_cumsum_sd'] = np.sqrt(volume_var_cumsum)
        
        # Absolute surface area distribution
        if 'area_nm^2_per_mL_avg' in result_df.columns:
            cumsum_avg = result_df.loc[scale_indices, 'area_nm^2_per_mL_avg'].cumsum()
            result_df.loc[scale_indices, 'area_nm^2_per_mL_cumsum_avg'] = cumsum_avg
            
            if 'area_nm^2_per_mL_sd' in result_df.columns:
                area_var_cumsum = (result_df.loc[scale_indices, 'area_nm^2_per_mL_sd'] ** 2).cumsum()
                result_df.loc[scale_indices, 'area_nm^2_per_mL_cumsum_sd'] = np.sqrt(area_var_cumsum)
    
    return result_df


def calculate_total_metrics_with_uncertainty(df, scale_column='scale'):
    """
    Calculate total metrics for each scale with proper uncertainty propagation.
    For totals: σ_total = √(Σ σ_i²)
    """
    results = {}
    
    for scale in df[scale_column].unique():
        scale_df = df[df[scale_column] == scale].sort_values('size_nm')
        scale_metrics = {}
        
        if not scale_df.empty:
            # Total particles per mL
            if 'particles_per_mL_avg' in scale_df.columns:
                total_particles_avg = scale_df['particles_per_mL_avg'].sum()
                scale_metrics['total_particles_per_mL_avg'] = total_particles_avg
                scale_metrics['total_particles_per_uL_avg'] = total_particles_avg / 1000
                
                if 'particles_per_mL_sd' in scale_df.columns:
                    total_particles_sd = np.sqrt((scale_df['particles_per_mL_sd'] ** 2).sum())
                    scale_metrics['total_particles_per_mL_sd'] = total_particles_sd
                    scale_metrics['total_particles_per_uL_sd'] = total_particles_sd / 1000
            
            # Total volume per mL
            if 'volume_nm^3_per_mL_avg' in scale_df.columns:
                total_volume_avg = scale_df['volume_nm^3_per_mL_avg'].sum()
                scale_metrics['total_volume_nm^3_per_mL_avg'] = total_volume_avg
                scale_metrics['total_volume_um^3_per_mL_avg'] = total_volume_avg / 1e9
                scale_metrics['total_volume_uL_per_mL_avg'] = total_volume_avg / 1e18
                scale_metrics['volume_percentage_avg'] = (total_volume_avg / 1e18) * 0.1
                
                if 'volume_nm^3_per_mL_sd' in scale_df.columns:
                    total_volume_sd = np.sqrt((scale_df['volume_nm^3_per_mL_sd'] ** 2).sum())
                    scale_metrics['total_volume_nm^3_per_mL_sd'] = total_volume_sd
                    scale_metrics['total_volume_um^3_per_mL_sd'] = total_volume_sd / 1e9
                    scale_metrics['total_volume_uL_per_mL_sd'] = total_volume_sd / 1e18
                    scale_metrics['volume_percentage_sd'] = (total_volume_sd / 1e18) * 0.1
            
            # Total surface area per mL
            if 'area_nm^2_per_mL_avg' in scale_df.columns:
                total_area_avg = scale_df['area_nm^2_per_mL_avg'].sum()
                scale_metrics['total_surface_area_nm^2_per_mL_avg'] = total_area_avg
                scale_metrics['total_surface_area_um^2_per_mL_avg'] = total_area_avg / 1e6
                scale_metrics['total_surface_area_cm^2_per_mL_avg'] = total_area_avg / 1e14
                
                if 'area_nm^2_per_mL_sd' in scale_df.columns:
                    total_area_sd = np.sqrt((scale_df['area_nm^2_per_mL_sd'] ** 2).sum())
                    scale_metrics['total_surface_area_nm^2_per_mL_sd'] = total_area_sd
                    scale_metrics['total_surface_area_um^2_per_mL_sd'] = total_area_sd / 1e6
                    scale_metrics['total_surface_area_cm^2_per_mL_sd'] = total_area_sd / 1e14
                
                # Specific surface area (derived, no uncertainty)
                if ('total_volume_nm^3_per_mL_avg' in scale_metrics and 
                    scale_metrics['total_volume_nm^3_per_mL_avg'] > 0):
                    ssa_1_per_nm = total_area_avg / scale_metrics['total_volume_nm^3_per_mL_avg']
                    scale_metrics['specific_surface_area_m^2_per_cm^3_avg'] = ssa_1_per_nm * 10
        
        results[scale] = scale_metrics
    
    return results


def add_metrics_to_metadata_with_uncertainty(metadata, metrics, scale='linear', num_replicates=None):
    """
    Add key metrics with uncertainties to the metadata dictionary.
    Only adds metrics for the specified scale (default: linear).
    
    Parameters:
    metadata (dict): Current metadata dictionary
    metrics (dict): Dictionary of calculated metrics with uncertainties
    scale (str): Which scale's metrics to use ('linear' or 'logarithmic')
    num_replicates (int): Number of replicates (optional)
    
    Returns:
    dict: Updated metadata dictionary with added metrics
    """
    updated_metadata = metadata.copy()
    
    # Use specified scale (default: linear)
    if scale not in metrics:
        if metrics:
            scale = list(metrics.keys())[0]
        else:
            return updated_metadata
    
    scale_metrics = metrics[scale]
    
    # Add key metrics to metadata with uncertainties
    if 'total_particles_per_mL_avg' in scale_metrics:
        avg_val = scale_metrics['total_particles_per_mL_avg']
        sd_val = scale_metrics.get('total_particles_per_mL_sd', 0)
        updated_metadata['nta_total_particles_per_mL'] = f"{avg_val:.2E} ± {sd_val:.2E}"
    
    if 'total_volume_uL_per_mL_avg' in scale_metrics:
        avg_val = scale_metrics['total_volume_uL_per_mL_avg']
        sd_val = scale_metrics.get('total_volume_uL_per_mL_sd', 0)
        updated_metadata['nta_total_volume_uL_per_mL'] = f"{avg_val:.4E} ± {sd_val:.4E}"
    
    if 'volume_percentage_avg' in scale_metrics:
        avg_val = scale_metrics['volume_percentage_avg']
        sd_val = scale_metrics.get('volume_percentage_sd', 0)
        # Use scientific notation for volume_percentage
        updated_metadata['nta_volume_percentage'] = f"{avg_val:.6E} ± {sd_val:.6E}"
    
    if 'specific_surface_area_m^2_per_cm^3_avg' in scale_metrics:
        avg_val = scale_metrics['specific_surface_area_m^2_per_cm^3_avg']
        updated_metadata['nta_specific_surface_area_m^2_per_cm^3'] = f"{avg_val:.2f}"
    
    # Add scale info (but NOT replicates - it's already in metadata)
    updated_metadata['nta_metrics_scale'] = scale
    
    return updated_metadata


# ============================================================================
# CELL 06: STATISTICS WITH BOUNDS-BASED UNCERTAINTY PROPAGATION
# ============================================================================

def interpolate_d_value_with_bounds(sizes, cumsum_avg, cumsum_sd, target_fraction):
    """
    Calculate D-value with asymmetric confidence bounds using bounds approach.
    
    Parameters:
    sizes (array): Size values (nm)
    cumsum_avg (array): Average cumulative distribution values
    cumsum_sd (array): Standard deviation of cumulative distribution values
    target_fraction (float): Target fraction (e.g., 0.1 for D10, 0.5 for D50)
    
    Returns:
    tuple: (d_value_avg, d_value_lower, d_value_upper)
    """
    # Ensure arrays are sorted by size
    sorted_indices = np.argsort(sizes)
    sizes_sorted = sizes[sorted_indices]
    cumsum_avg_sorted = cumsum_avg[sorted_indices]
    cumsum_sd_sorted = cumsum_sd[sorted_indices]
    
    # Calculate bounds
    cumsum_lower = cumsum_avg_sorted - cumsum_sd_sorted
    cumsum_upper = cumsum_avg_sorted + cumsum_sd_sorted
    
    # Ensure cumulative distributions are monotonic and in valid range [0,1]
    cumsum_avg_sorted = np.clip(cumsum_avg_sorted, 0, 1)
    cumsum_lower = np.clip(cumsum_lower, 0, 1)
    cumsum_upper = np.clip(cumsum_upper, 0, 1)
    
    # Ensure monotonicity by taking cumulative maximum
    cumsum_avg_sorted = np.maximum.accumulate(cumsum_avg_sorted)
    cumsum_lower = np.maximum.accumulate(cumsum_lower)
    cumsum_upper = np.maximum.accumulate(cumsum_upper)
    
    # Check if target fraction is achievable
    if target_fraction < cumsum_avg_sorted[0] or target_fraction > cumsum_avg_sorted[-1]:
        return np.nan, np.nan, np.nan
    
    # Interpolate D-values
    try:
        d_value_avg = np.interp(target_fraction, cumsum_avg_sorted, sizes_sorted)
        d_value_lower = np.interp(target_fraction, cumsum_upper, sizes_sorted)  # Note: swapped!
        d_value_upper = np.interp(target_fraction, cumsum_lower, sizes_sorted)  # Note: swapped!
        
        # The swapping is because:
        # - When cumsum is higher (cumsum + sd), we reach target fraction at smaller size → lower bound
        # - When cumsum is lower (cumsum - sd), we reach target fraction at larger size → upper bound
        
        return d_value_avg, d_value_lower, d_value_upper
        
    except Exception as e:
        print(f"Error in interpolation for target fraction {target_fraction}: {e}")
        return np.nan, np.nan, np.nan


def calculate_percentile_statistics_with_uncertainty(df, size_column='size_nm'):
    """
    Calculate percentile statistics (D10, D50, D90, span) with uncertainties for all
    cumulative distributions (number, volume, surface area) and scales.
    
    Parameters:
    df (DataFrame): DataFrame containing size and cumulative distribution data with uncertainties
    size_column (str): Column name for particle sizes
    
    Returns:
    dict: Nested dictionary of statistics by scale, distribution type, and metric
    """
    # Initialize results dictionary
    stats = {'linear': {}, 'logarithmic': {}}
    
    # Define cumulative distribution configurations
    cumsum_configs = [
        {
            'name': 'number',
            'avg_column': 'number_normalized_cumsum_avg',
            'sd_column': 'number_normalized_cumsum_sd'
        },
        {
            'name': 'volume',
            'avg_column': 'volume_nm^3_per_mL_cumsum_avg', 
            'sd_column': 'volume_nm^3_per_mL_cumsum_sd'
        },
        {
            'name': 'surface_area',
            'avg_column': 'area_nm^2_per_mL_cumsum_avg',
            'sd_column': 'area_nm^2_per_mL_cumsum_sd'
        }
    ]
    
    # Process each scale
    for scale in ['linear', 'logarithmic']:
        print(f"\nCalculating statistics for {scale.upper()} scale:")
        
        # Filter data for this scale
        scale_df = df[df['scale'] == scale].copy()
        if scale_df.empty:
            print(f"  No data available for {scale} scale")
            continue
            
        # Sort data by size
        scale_df = scale_df.sort_values(size_column)
        
        # Store skip reasons for debugging
        skip_reasons = []
        
        # Process each cumulative distribution type
        for config in cumsum_configs:
            name = config['name']
            avg_column = config['avg_column']
            sd_column = config['sd_column']
            
            print(f"  Processing {name}-weighted distribution...")
            
            # Skip if required columns don't exist
            if avg_column not in scale_df.columns or sd_column not in scale_df.columns:
                reason = f"Missing columns for {name}: {avg_column} or {sd_column}"
                print(f"    ✗ {reason}")
                skip_reasons.append(reason)
                if name == 'number':
                    print(f"      Available columns: {list(scale_df.columns)}")
                continue
            
            # Extract data arrays
            sizes = scale_df[size_column].values
            cumsum_avg = scale_df[avg_column].values
            cumsum_sd = scale_df[sd_column].values
            
            # Skip if all values are zero or NaN
            if np.all(cumsum_avg == 0) or np.all(np.isnan(cumsum_avg)):
                reason = f"All values are zero or NaN for {name}: min={np.nanmin(cumsum_avg)}, max={np.nanmax(cumsum_avg)}"
                print(f"    ✗ {reason}")
                skip_reasons.append(reason)
                continue
            
            # For absolute distributions (volume, surface area), normalize to 0-1 for D-value calculation
            if name in ['volume', 'surface_area']:
                max_cumsum = np.nanmax(cumsum_avg)
                if max_cumsum > 0:
                    cumsum_avg = cumsum_avg / max_cumsum
                    cumsum_sd = cumsum_sd / max_cumsum
                else:
                    reason = f"Maximum cumsum is zero for {name}"
                    print(f"    ✗ {reason}")
                    skip_reasons.append(reason)
                    continue
            
            # Calculate D10, D50, D90 with asymmetric bounds
            try:
                d10_avg, d10_lower, d10_upper = interpolate_d_value_with_bounds(sizes, cumsum_avg, cumsum_sd, 0.1)
                d50_avg, d50_lower, d50_upper = interpolate_d_value_with_bounds(sizes, cumsum_avg, cumsum_sd, 0.5)
                d90_avg, d90_lower, d90_upper = interpolate_d_value_with_bounds(sizes, cumsum_avg, cumsum_sd, 0.9)
                
                # Check if results are valid
                if np.isnan(d10_avg) or np.isnan(d50_avg) or np.isnan(d90_avg):
                    reason = f"D-values are NaN for {name}: d10={d10_avg}, d50={d50_avg}, d90={d90_avg}"
                    print(f"    ✗ {reason}")
                    skip_reasons.append(reason)
                    continue
                
                # Calculate span with bounds
                # span = (D90 - D10) / D50
                if not np.isnan(d10_avg) and not np.isnan(d50_avg) and not np.isnan(d90_avg) and d50_avg > 0:
                    span_avg = (d90_avg - d10_avg) / d50_avg
                    
                    # For span bounds, consider which combinations give min/max
                    possible_spans = [
                        (d90_lower - d10_upper) / d50_upper,
                        (d90_lower - d10_upper) / d50_lower,
                        (d90_upper - d10_lower) / d50_upper,
                        (d90_upper - d10_lower) / d50_lower
                    ]
                    
                    # Filter out invalid spans
                    valid_spans = [s for s in possible_spans if not np.isnan(s) and np.isfinite(s)]
                    
                    if valid_spans:
                        span_lower = min(valid_spans)
                        span_upper = max(valid_spans)
                    else:
                        span_lower, span_upper = np.nan, np.nan
                else:
                    span_avg, span_lower, span_upper = np.nan, np.nan, np.nan
                
                # Store results with bounds
                metrics = {
                    'D10_avg': d10_avg,
                    'D10_lower': d10_lower,
                    'D10_upper': d10_upper,
                    'D50_avg': d50_avg,
                    'D50_lower': d50_lower,
                    'D50_upper': d50_upper,
                    'D90_avg': d90_avg,
                    'D90_lower': d90_lower,
                    'D90_upper': d90_upper,
                    'span_avg': span_avg,
                    'span_lower': span_lower,
                    'span_upper': span_upper
                }
                
                stats[scale][name] = metrics
                
                print(f"    ✓ D10: {d10_avg:.2f} nm ({d10_lower:.2f} - {d10_upper:.2f})")
                print(f"    ✓ D50: {d50_avg:.2f} nm ({d50_lower:.2f} - {d50_upper:.2f})")
                print(f"    ✓ D90: {d90_avg:.2f} nm ({d90_lower:.2f} - {d90_upper:.2f})")
                print(f"    ✓ Span: {span_avg:.3f} ({span_lower:.3f} - {span_upper:.3f})")
                
            except Exception as e:
                reason = f"Exception for {name}: {str(e)}"
                print(f"    ✗ {reason}")
                skip_reasons.append(reason)
                import traceback
                traceback.print_exc()
        
        # Store skip reasons in stats for debugging
        if skip_reasons:
            stats[scale]['_skip_reasons'] = skip_reasons
    
    return stats


def add_key_statistics_to_metadata(metadata, stats):
    """
    Add key D-values and span for all distributions to metadata.
    
    Parameters:
    metadata (dict): Current metadata dictionary
    stats (dict): Dictionary with calculated statistics
    
    Returns:
    dict: Updated metadata dictionary
    """
    # Create a copy to avoid modifying the original
    updated_metadata = metadata.copy()
    
    # Add debug info temporarily
    debug_info = []
    debug_info.append(f"stats keys: {list(stats.keys())}")
    
    # Get linear statistics
    if 'linear' in stats:
        linear_stats = stats['linear']
        
        # Check for skip reasons
        if '_skip_reasons' in linear_stats:
            debug_info.append(f"SKIP REASONS: {linear_stats['_skip_reasons']}")
        
        debug_info.append(f"linear_stats keys: {list([k for k in linear_stats.keys() if k != '_skip_reasons'])}")
        
        # Define distribution types
        dist_types = [
            ('number', 'linear_number'),
            ('volume', 'linear_volume'),
            ('surface_area', 'linear_surface_area')
        ]
        
        # Add D-values and span for each distribution type
        for dist_key, field_prefix in dist_types:
            if dist_key in linear_stats:
                debug_info.append(f"✓ {dist_key} FOUND in linear_stats")
                dist_stats = linear_stats[dist_key]
                
                # Add D10, D50, D90
                for param in ['D10', 'D50', 'D90']:
                    avg_val = dist_stats.get(f'{param}_avg', np.nan)
                    lower_val = dist_stats.get(f'{param}_lower', np.nan)
                    upper_val = dist_stats.get(f'{param}_upper', np.nan)
                    
                    if not np.isnan(avg_val):
                        updated_metadata[f'nta_{field_prefix}_{param.lower()}'] = f"{avg_val:.2f} nm ({lower_val:.2f} - {upper_val:.2f})"
                        debug_info.append(f"  ✓ {field_prefix}_{param.lower()} = {avg_val:.2f}")
                    else:
                        updated_metadata[f'nta_{field_prefix}_{param.lower()}'] = "Not available"
                        debug_info.append(f"  ✗ {field_prefix}_{param.lower()} = NaN")
                
                # Add span with bounds
                span_avg = dist_stats.get('span_avg', np.nan)
                span_lower = dist_stats.get('span_lower', np.nan)
                span_upper = dist_stats.get('span_upper', np.nan)
                if not np.isnan(span_avg):
                    updated_metadata[f'nta_{field_prefix}_span'] = f"{span_avg:.3f} ({span_lower:.3f} - {span_upper:.3f})"
                    debug_info.append(f"  ✓ {field_prefix}_span = {span_avg:.3f}")
                else:
                    updated_metadata[f'nta_{field_prefix}_span'] = "Not available"
                    debug_info.append(f"  ✗ {field_prefix}_span = NaN")
            else:
                debug_info.append(f"✗ {dist_key} NOT FOUND in linear_stats - Check calculate_percentile_statistics_with_uncertainty function output above")
    else:
        debug_info.append("✗ 'linear' key not in stats!")
    
    # Add debug info to metadata temporarily (will show in Metadata tab)
    updated_metadata['_debug_d_value_calculation'] = " | ".join(debug_info)
    
    return updated_metadata


# ============================================================================
# MAIN ANALYZER CLASS
# ============================================================================

class NTAAnalyzer:
    """
    Main NTA analysis class integrating all cells 01-06.
    
    Complete workflow:
    1. Read multiple files (Cell 02)
    2. Extract and average distributions (Cell 03)
    3. Extract and analyze metadata from all files (Cell 04)
    4. Apply dilution correction, normalization, cumulative distributions (Cell 05)
    5. Calculate percentile statistics with uncertainty (D-values, span) (Cell 06)
    6. Generate standardized outputs
    """
    
    def __init__(self, config=None):
        self.config = config or CONFIG
        self.results = {}
    
    def process(self, filepaths):
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
            if dist[0]:  # success flag
                distributions.append(dist[1])
        
        if not distributions:
            raise Exception("Could not extract distribution data from any files")
        
        # Average distributions
        filenames = [os.path.basename(f[0]) for f in successful_files]
        success, avg_dist = average_replicate_data(distributions, filenames)
        
        if not success:
            raise Exception(f"Failed to average data: {avg_dist}")
        
        # ===== CELL 04: Metadata Analysis =====
        success, all_metadata = extract_metadata_from_all_files(successful_files)
        
        if not success:
            raise Exception(f"Failed to extract metadata: {all_metadata}")
        
        # Analyze differences
        identical, different, analysis = analyze_field_differences(all_metadata)
        
        # Create standardized metadata
        metadata, quality_alerts, high_variation = create_automated_metadata(
            all_metadata, identical, different, self.config
        )
        
        # ===== CELL 05: Dilution Correction, Normalization & Metrics =====
        print("\n" + "="*80)
        print("EXECUTING CELL 05: DILUTION CORRECTION & NORMALIZATION")
        print("="*80)
        
        # Step 1: Dilution correction
        print("\n1. Applying dilution correction with uncertainty propagation...")
        success, dilution_corrected_df = apply_dilution_correction_with_uncertainty(
            avg_dist, 
            metadata=metadata
        )
        if success:
            print("✓ Dilution correction completed")
            avg_dist = dilution_corrected_df
        else:
            print(f"✗ Dilution correction failed: {dilution_corrected_df}")
        
        # Step 2: Normalization
        print("\n2. Normalizing distributions with uncertainty propagation...")
        try:
            normalized_df = normalize_distributions_with_uncertainty(avg_dist)
            print("✓ Normalization completed")
            avg_dist = normalized_df
        except Exception as e:
            print(f"✗ Normalization failed: {e}")
        
        # Step 3: Cumulative distributions
        print("\n3. Calculating cumulative distributions with uncertainty...")
        try:
            final_df = calculate_cumulative_distributions_with_uncertainty(avg_dist)
            print("✓ Cumulative distributions calculated")
            avg_dist = final_df
        except Exception as e:
            print(f"✗ Cumulative distribution calculation failed: {e}")
        
        # Step 4: Total metrics
        print("\n4. Calculating total metrics with uncertainty...")
        try:
            total_metrics = calculate_total_metrics_with_uncertainty(avg_dist)
            print("✓ Total metrics calculated")
        except Exception as e:
            print(f"✗ Total metrics calculation failed: {e}")
            total_metrics = {}
        
        # Step 5: Add metrics to metadata for saving
        print("\n5. Adding metrics to metadata...")
        try:
            metadata = add_metrics_to_metadata_with_uncertainty(
                metadata, total_metrics, scale='linear', num_replicates=len(successful_files)
            )
            # Add analysis date
            metadata['nta_python_analysis'] = str(date.today())
            print("✓ Metrics added to metadata")
        except Exception as e:
            print(f"✗ Failed to add metrics to metadata: {e}")
        
        print("\n" + "="*80)
        print("CELL 05 PROCESSING COMPLETE")
        print("="*80 + "\n")
        
        # ===== CELL 06: STATISTICS WITH UNCERTAINTY PROPAGATION =====
        print("="*80)
        print("EXECUTING CELL 06: STATISTICS WITH UNCERTAINTY PROPAGATION")
        print("="*80)
        
        # Step 1: Calculate percentile statistics (D10, D50, D90, span)
        print("\n1. Calculating percentile statistics with uncertainty...")
        try:
            percentile_stats = calculate_percentile_statistics_with_uncertainty(avg_dist)
            print("✓ Percentile statistics calculated")
        except Exception as e:
            print(f"✗ Percentile statistics calculation failed: {e}")
            percentile_stats = {}
        
        # Step 2: Add key D-values to metadata
        print("\n2. Adding key D-values to metadata...")
        print(f"   percentile_stats available: {bool(percentile_stats)}")
        print(f"   percentile_stats keys: {list(percentile_stats.keys()) if percentile_stats else 'EMPTY'}")
        try:
            metadata_before = len(metadata)
            metadata = add_key_statistics_to_metadata(metadata, percentile_stats)
            metadata_after = len(metadata)
            print(f"   Metadata fields added: {metadata_after - metadata_before}")
            
            # Show D-values that were added
            d_values = [k for k in metadata.keys() if any(p in k.lower() for p in ['d10', 'd50', 'd90', 'span'])]
            print(f"   D-value fields in metadata: {len(d_values)}")
            for field in sorted(d_values):
                print(f"     - {field}: {metadata[field]}")
            
            print(f"✓ D-values added to metadata")
        except Exception as e:
            print(f"✗ Failed to add D-values to metadata: {e}")
            import traceback
            traceback.print_exc()
        
        print("\n" + "="*80)
        print("CELL 06 PROCESSING COMPLETE")
        print("="*80 + "\n")
        
        self.results = {
            'distribution': avg_dist,
            'metadata': metadata,
            'total_metrics': total_metrics,
            'percentile_stats': percentile_stats,
            'identical_fields': identical,
            'different_fields': different,
            'field_analysis': analysis,
            'all_file_metadata': all_metadata,
            'quality_alerts': quality_alerts,
            'high_variation_fields': high_variation,
            'failed_files': failed_files,
            'num_replicates': len(successful_files)
        }
        
        return self.results
    
    def save_outputs(self, output_dir):
        """Save all outputs with proper organization."""
        if 'distribution' not in self.results:
            raise Exception("No results. Run process() first.")
        
        os.makedirs(output_dir, exist_ok=True)
        
        unique_id = self.results['metadata'].get('persistentID', 'analysis')
        created_files = []
        
        # ===== Save Metadata (with headers - human readable) =====
        success, meta_result = save_metadata_file(self.results['metadata'], output_dir, self.config, include_headers=True)
        if success:
            created_files.append(meta_result)
        else:
            raise Exception(f"Failed to save metadata: {meta_result}")
        
        # ===== Save Metadata (without headers - Excel import friendly) =====
        success, meta_result_excel = save_metadata_file(self.results['metadata'], output_dir, self.config, include_headers=False)
        if success:
            created_files.append(meta_result_excel)
        
        
        # ===== Save Distribution Data (Split by Scale) =====
        dist = self.results['distribution']
        
        # Linear scale
        linear_data = dist[dist['scale'] == 'linear'].copy()
        if not linear_data.empty:
            linear_path = os.path.join(output_dir, f'Data_{unique_id}_PSD_LINEAR.txt')
            
            # Save data without headers
            linear_data_export = linear_data.drop('scale', axis=1)
            linear_data_export.to_csv(linear_path, sep='\t', index=False)
            
            created_files.append(linear_path)
        
        # Logarithmic scale
        log_data = dist[dist['scale'] == 'logarithmic'].copy()
        if not log_data.empty:
            log_path = os.path.join(output_dir, f'Data_{unique_id}_PSD_LOGARITHMIC.txt')
            
            # Save data without headers
            log_data_export = log_data.drop('scale', axis=1)
            log_data_export.to_csv(log_path, sep='\t', index=False)
            
            created_files.append(log_path)
        
        return created_files
