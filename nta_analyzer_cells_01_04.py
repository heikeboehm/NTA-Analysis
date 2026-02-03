"""
NTA Data Analysis - Complete Integrated Script (Cells 01-04)

This script combines all functions from your original notebook:
- Cell 01: Configuration & utilities
- Cell 02: File I/O (read, identify sections, validate)
- Cell 03: Data extraction & averaging
- Cell 04: Metadata extraction & analysis

All code is your EXACT notebook code, adapted for Streamlit.
"""

import os
import re
import json
import ntpath
from datetime import date

import pandas as pd
import numpy as np


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
            "funding": "Funding_Source",
            "data_collection_method": "NTA",
            "unit_of_analysis": '["nm", "nm^2", "nm^3"]',
            "keywords": "particle_size_distribution",
            "publications": "None",
        }
    }


def extract_sample_info(filename):
    """
    Extract sample information from filename for more readable display.
    (From Cell 01)
    """
    base_name = filename.replace("_rawdata.txt", "").replace("size_NTA", "")
    
    if base_name.startswith("Data_"):
        base_name = base_name[5:]
    
    date_match = re.search(r'_(\d{8})_', base_name)
    date_str = ""
    if date_match:
        date = date_match.group(1)
        date_str = f" (Date: {date[0:4]}-{date[4:6]}-{date[6:8]})"
    
    return f"{base_name}{date_str}"


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
    
    metadata = {}
    
    for key, pattern in metadata_patterns:
        match = re.search(pattern, content)
        if match:
            value = match.group(1).strip() if match.group(1) else ''
            if value and not value.lower() in ['none', 'null', '']:
                metadata[key] = value
    
    metadata['filename'] = filename
    
    base_name = os.path.splitext(filename)[0]
    if base_name.endswith("_rawdata"):
        base_name = base_name[:-8]
    if base_name.startswith("Data_"):
        base_name = base_name[5:]
    metadata['uniqueID'] = base_name
    
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
        'number_of_traces', 'detected_particles'
    ]
    
    instrument_determined_fields = [
        'median_number_d50', 'median_concentration_d50', 'median_volume_d50'
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
    
    elif field_name in instrument_determined_fields:
        return json.dumps(values), "instrument_determined_per_measurement"
    
    elif field_name in quality_control_fields:
        unique_values = list(set(values))
        if len(unique_values) == 1:
            return values[0], "qc_consistent"
        else:
            return json.dumps(values), f"QC_ALERT_inconsistent_values"
    
    else:
        try:
            if field_name == 'avi_filesize':
                numeric_values = []
                for v in values:
                    if isinstance(v, str) and 'MB' in v:
                        numeric_values.append(float(v.replace(' MB', '')))
                    else:
                        numeric_values.append(float(v))
                
                mean_val = np.mean(numeric_values)
                std_val = np.std(numeric_values, ddof=1) if len(numeric_values) > 1 else 0.0
                
                cv = (std_val / mean_val * 100) if mean_val != 0 else 0
                if cv > 10:
                    notes = f"HIGH_VARIATION_CV_{cv:.1f}%"
                else:
                    notes = f"mean_sd_of_{len(values)}"
                
                return f"{mean_val:.2f} ± {std_val:.2f} MB", notes
            
            numeric_values = [float(v) for v in values]
            mean_val = np.mean(numeric_values)
            std_val = np.std(numeric_values, ddof=1) if len(numeric_values) > 1 else 0.0
            
            cv = (std_val / mean_val * 100) if mean_val != 0 else 0
            
            if cv > 10:
                notes = f"HIGH_VARIATION_CV_{cv:.1f}%"
            else:
                notes = f"mean_sd_of_{len(values)}"
            
            return f"{mean_val:.2f} ± {std_val:.2f}", notes
            
        except ValueError:
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
    
    metadata['experimenter'] = project_meta.get('experimenter', 'Your_Initials')
    metadata['location'] = project_meta.get('location', 'Your_Lab_Location')
    metadata['project'] = project_meta.get('project', 'Your_Project_Name')
    metadata['meta_version'] = project_meta.get('meta_version', 'v03')
    metadata['pi'] = project_meta.get('pi', 'Your_PI_Initials')
    metadata['funding'] = project_meta.get('funding', 'Funding_Source')
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
    
    essential_fields = {
        'date', 'temperature', 'ph', 'dilution', 'laser_wavelength', 'electrolyte',
        'positions', 'cycles', 'fps', 
        'particle_drift_check_result', 'cell_check_result',
        'average_number_of_particles', 'number_of_traces', 'detected_particles',
        'conductivity', 'scattering_intensity', 'viscosity',
        'avi_filesize'
    }
    
    for field_name, value in identical_fields.items():
        if field_name in essential_fields:
            if field_name in ['temperature', 'ph', 'dilution', 'laser_wavelength', 'positions', 
                             'cycles', 'fps', 'average_number_of_particles', 'number_of_traces', 
                             'detected_particles', 'particle_drift_check_result', 'cell_check_result',
                             'conductivity', 'scattering_intensity', 'viscosity', 'avi_filesize']:
                if field_name in ['number_of_traces', 'detected_particles']:
                    metadata[f'nta_{field_name}_sum'] = value
                else:
                    metadata[f'nta_{field_name}'] = value
            else:
                metadata[field_name] = value
    
    quality_alerts = []
    high_variation_fields = []
    
    for field_name, values in different_fields.items():
        if field_name in essential_fields:
            formatted_value, notes = smart_format_field(field_name, values)
            
            if "QC_ALERT" in notes:
                quality_alerts.append(f"{field_name}: {formatted_value}")
            if "HIGH_VARIATION" in notes:
                high_variation_fields.append(f"{field_name}: {notes}")
            
            if field_name in ['temperature', 'ph', 'dilution', 'laser_wavelength', 'positions', 
                             'cycles', 'fps', 'average_number_of_particles', 'number_of_traces', 
                             'detected_particles', 'particle_drift_check_result', 'cell_check_result',
                             'conductivity', 'scattering_intensity', 'viscosity', 'avi_filesize']:
                if field_name in ['number_of_traces', 'detected_particles']:
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


def save_metadata_file(metadata, output_dir=None, config=None):
    """
    Save metadata to a file in the specified output directory.
    (From Cell 04)
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
    metadata_path = os.path.join(output_dir, f"Data_{unique_id}_metadata.txt")
    
    try:
        with open(metadata_path, 'w') as f:
            for key, value in metadata.items():
                f.write(f"{key}\t{value}\t\n")
        
        return True, metadata_path
    except Exception as e:
        return False, f"Failed to write metadata file: {str(e)}"


# ============================================================================
# MAIN ANALYZER CLASS
# ============================================================================

class NTAAnalyzer:
    """
    Main NTA analysis class integrating all cells 01-04.
    
    Complete workflow:
    1. Read multiple files (Cell 02)
    2. Extract and average distributions (Cell 03)
    3. Extract and analyze metadata from all files (Cell 04)
    4. Generate standardized outputs
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
        
        self.results = {
            'distribution': avg_dist,
            'metadata': metadata,
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
        """Save all outputs."""
        if 'distribution' not in self.results:
            raise Exception("No results. Run process() first.")
        
        os.makedirs(output_dir, exist_ok=True)
        
        unique_id = self.results['metadata'].get('persistentID', 'analysis')
        created_files = []
        
        # Save metadata
        meta_path = os.path.join(output_dir, f'Data_{unique_id}_metadata.txt')
        with open(meta_path, 'w') as f:
            for key, value in sorted(self.results['metadata'].items()):
                f.write(f"{key}\t{value}\t\n")
        created_files.append(meta_path)
        
        # Save distribution
        dist_path = os.path.join(output_dir, f'Data_{unique_id}_PSD.txt')
        self.results['distribution'].to_csv(dist_path, sep='\t', index=False)
        created_files.append(dist_path)
        
        return created_files
