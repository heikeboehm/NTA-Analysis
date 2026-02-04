import os
import re
import json
import ntpath
from datetime import date

import pandas as pd
import numpy as np
from scipy import integrate


"""
NTA Data Analysis - Multi-file Selection Module (Cell 01)

This module provides file selection for single or multiple NTA files
for replicate analysis with averaging and standard deviation calculations.
"""

import os
import re

if 'CONFIG' not in globals():
    CONFIG = {
        "directory": "Inbox",  # Update this path
        "file_identifier": ".txt",
        "output_subdirs": ["metadata", "processed"],
        "nta_concentration_calibration_factor": 4.61E+5,  # Default ZetaView calibration
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

def set_data_directory(directory_path):
    """
    Set or update the data directory in CONFIG.
    
    Parameters:
    directory_path (str): Path to directory containing NTA data
    
    Returns:
    bool: True if directory exists and was set, False otherwise
    """
    global CONFIG
    
    # Validate directory exists
    if not os.path.exists(directory_path):
        print(f"Error: Directory not found: {directory_path}")
        return False
    
    # Update CONFIG with new directory
    CONFIG["directory"] = directory_path
    
    # Ensure output subdirectories exist
    for subdir in CONFIG["output_subdirs"]:
        output_dir = os.path.join(directory_path, subdir)
        os.makedirs(output_dir, exist_ok=True)
    
    print(f"Data directory set to: {directory_path}")
    print(f"Output directories created/verified.")
    
    return True

def find_nta_files():
    """
    Find NTA data files in the current directory.
    
    Returns:
    list: List of NTA files, or empty list if none found
    """
    directory = CONFIG["directory"]
    file_identifier = CONFIG["file_identifier"]
    
    # Validate directory exists
    if not os.path.exists(directory):
        print(f"Error: Directory not found: {directory}")
        return []
    
    # Find all matching files
    all_files = os.listdir(directory)
    nta_files = [f for f in all_files if f.endswith(file_identifier)]
    
    return nta_files

def extract_sample_info(filename):
    """
    Extract sample information from filename for more readable display.
    
    Parameters:
    filename (str): Filename to analyze
    
    Returns:
    str: Formatted sample information
    """
    # Remove common prefixes and suffixes
    base_name = filename.replace("_rawdata.txt", "").replace("size_NTA", "")
    
    if base_name.startswith("Data_"):
        base_name = base_name[5:]
    
    # Look for date pattern (e.g., 20250311)
    date_match = re.search(r'_(\d{8})_', base_name)
    date_str = ""
    if date_match:
        date = date_match.group(1)
        date_str = f" (Date: {date[0:4]}-{date[4:6]}-{date[6:8]})"
    
    return f"{base_name}{date_str}"

def parse_indices(indices_input):
    """
    Parse file indices from string input.
    
    Parameters:
    indices_input (str): Comma-separated indices (e.g., "0,1,2" or "0")
    
    Returns:
    list: List of integer indices, or empty list if invalid
    """
    try:
        # Remove spaces and split by commas
        indices_str = indices_input.replace(" ", "").split(",")
        indices = [int(idx) for idx in indices_str if idx.strip()]
        return indices
    except ValueError:
        print(f"Error: Invalid indices format. Use comma-separated numbers (e.g., '0,1,2')")
        return []

def select_files(indices_input):
    """
    Select multiple files by indices for replicate analysis.
    
    Parameters:
    indices_input (str): Comma-separated string of indices (e.g., "0,1,2")
    
    Returns:
    tuple: (filenames_list, filepaths_list) or (None, None) if invalid
    """
    nta_files = find_nta_files()
    
    if not nta_files:
        print("No NTA files found in the current directory.")
        return None, None
    
    # Parse indices
    indices = parse_indices(indices_input)
    if not indices:
        return None, None
    
    # Validate indices
    invalid_indices = [idx for idx in indices if idx < 0 or idx >= len(nta_files)]
    if invalid_indices:
        print(f"Invalid file indices: {invalid_indices}. Must be between 0 and {len(nta_files)-1}")
        return None, None
    
    # Remove duplicates while preserving order
    unique_indices = []
    for idx in indices:
        if idx not in unique_indices:
            unique_indices.append(idx)
    
    # Select files
    selected_filenames = []
    selected_filepaths = []
    
    for idx in unique_indices:
        filename = nta_files[idx]
        filepath = os.path.join(CONFIG["directory"], filename)
        selected_filenames.append(filename)
        selected_filepaths.append(filepath)
    
    # Display selection
    print(f"\nSelected {len(selected_filenames)} file(s) for analysis:")
    for i, (idx, filename) in enumerate(zip(unique_indices, selected_filenames)):
        print(f"  File {i+1}: #{idx} - {filename}")
    
    # Set global variables for use in other cells (always as lists)
    global selected_filenames_list, selected_filepaths_list, num_replicates
    selected_filenames_list = selected_filenames
    selected_filepaths_list = selected_filepaths
    num_replicates = len(selected_filenames)
    
    print(f"\nNumber of replicates: {num_replicates}")
    if num_replicates == 1:
        print("Single file analysis - SD values will be 0")
    else:
        print("Multi-file analysis - averaging with SD calculations")
    
    return selected_filenames, selected_filepaths

# Set a custom directory (uncomment and modify if needed)
# set_data_directory("/path/to/your/data")

# List available files
print("=" * 80)
print("NTA DATA ANALYSIS - MULTI-FILE SELECTION")
print("=" * 80)
print(f"Current directory: {CONFIG['directory']}")
print()

nta_files = find_nta_files()
if nta_files:
    print(f"Found {len(nta_files)} NTA files:")
    for i, file in enumerate(nta_files):
        print(f"  {i}: {file}")
    
    print("\nTo select files belonging to the same sample for analysis:")
    print("  Single file: select_files('0')")
    print("  Multiple files: select_files('0,1,2')")
    print("  Example: select_files('0,2,5') for files 0, 2, and 5")
else:
    print("No NTA files found in the current directory.")
    print("Please set a different directory using set_data_directory()")

print("=" * 80)

# Initialize global variables for file selection (always as lists)
selected_filenames_list = None
selected_filepaths_list = None
num_replicates = 0

# Default selection - uncomment to automatically select the first file
# if nta_files:
#     select_files('0')

"""
NTA Data Analysis - Multi-file I/O Module (Cell 02)

This module handles:
1. Reading multiple NTA data files
2. Identifying key sections in file content
3. Validating file structure before processing
4. Handling errors gracefully by skipping problematic files

"""

import os
import re


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
    filepath (str): Path to the NTA file
    
    Returns:
    tuple: (success_flag, result)
        - If successful: (True, (content, sections))
        - If failed: (False, error_message)
    """
    # Read the file
    success, result = read_nta_file(filepath)
    if not success:
        return False, result
    content = result
    
    # Identify sections
    success, result = identify_sections(content)
    if not success:
        return False, result
    sections = result
    
    # Validate structure
    success, message = validate_file_structure(content, sections)
    if not success:
        return False, message
    
    return True, (content, sections)


def process_multiple_file_contents(filepaths):
    """
    Process multiple NTA files, skipping any that fail with error messages.
    
    Parameters:
    filepaths (list): List of file paths to process
    
    Returns:
    tuple: (file_data_list, error_summary)
        - file_data_list: List of (filename, content, sections) for successful files
        - error_summary: Dict with failed files and their error messages
    """
    successful_files = []
    failed_files = {}
    
    print(f"Processing {len(filepaths)} file(s)...")
    
    for filepath in filepaths:
        filename = os.path.basename(filepath)
        print(f"  Processing: {filename}")
        
        success, result = process_single_file_content(filepath)
        
        if success:
            content, sections = result
            successful_files.append((filename, content, sections))
            print(f"    ✓ Success")
        else:
            failed_files[filename] = result
            print(f"    ✗ Failed: {result}")
    
    # Summary
    print(f"\nProcessing complete:")
    print(f"  Successful: {len(successful_files)} files")
    print(f"  Failed: {len(failed_files)} files")
    
    if failed_files:
        print("\nFailed files:")
        for filename, error in failed_files.items():
            print(f"  {filename}: {error}")
    
    return successful_files, failed_files


def preview_file_content(content, max_chars=500):
    """
    Display a preview of the file content to verify it was read correctly.
    
    Parameters:
    content (str): File content to preview
    max_chars (int): Maximum number of characters to display
    
    Returns:
    str: Preview text for display
    """
    preview = content[:max_chars]
    return "-" * 50 + "\n" + preview + "\n" + "-" * 50


# Execute file processing if files were selected in Cell 00
if 'selected_filepaths_list' in globals() and selected_filepaths_list:
    print("=" * 80)
    print("PROCESSING SELECTED NTA FILES")
    print("=" * 80)
    
    # Process all selected files
    successful_files, failed_files = process_multiple_file_contents(selected_filepaths_list)
    
    if successful_files:
        print(f"\nSuccessfully processed {len(successful_files)} file(s):")
        for filename, content, sections in successful_files:
            print(f"  {filename}")
            print(f"    Linear section at position: {sections['linear_start']}")
            print(f"    Logarithmic section at position: {sections['logarithmic_start']}")
        
        # Show preview of first file
        if successful_files:
            first_filename, first_content, first_sections = successful_files[0]
            print(f"\nFile Preview ({first_filename}):")
            print(preview_file_content(first_content))
        
        # Store results for use in subsequent cells
        current_files_data = successful_files
        current_failed_files = failed_files
        
        print("\nFile processing completed successfully!")
        print("Ready for data extraction and averaging (Cell 03)")
        
    else:
        print("\nERROR: No files could be processed successfully.")
        if failed_files:
            print("All files failed with errors (see details above).")
        
else:
    print("No files selected. Please run Cell 00 (file selection) first.")
    print("Use: select_files('0,1,2') to select files for analysis.")

"""
NTA Data Analysis - Data Extraction and Averaging Module (Cell 03)

This module handles:
1. Extracting distribution data from multiple NTA files
2. Averaging raw particle counts bin-by-bin across replicates
3. Calculating standard deviations for each size bin
4. Creating a single averaged dataset for downstream analysis

"""

import pandas as pd
import numpy as np
import re


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
    
    # For logarithmic section, skip past the separator line
    if is_log_section:
        sep_pos = section.find("-1.000E+0")
        if sep_pos != -1:
            sep_line_end = section.find('\n', sep_pos)
            if sep_line_end != -1:
                section = section[sep_line_end + 1:]
            else:
                section = section[sep_pos + 20:]
        
        # Sometimes there's a second "Size Distribution" header
        second_header = section.find("Size Distribution")
        if second_header != -1:
            section = section[second_header:]
    
    # Find the header line
    header_match = re.search(r'Size / nm\s+Number\s+Concentration.+', section)
    
    if not header_match and is_log_section:
        # For log section, try more relaxed pattern or use default
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
    
    # Extract data lines
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
    
    Parameters:
    header_line (str): Header line with column names
    data_lines (list): List of data lines to parse
    scale_type (str): Type of scale ('linear' or 'logarithmic')
    
    Returns:
    tuple: (success_flag, result)
        - If successful: (True, DataFrame)
        - If failed: (False, error_message)
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
            except ValueError as e:
                print(f"Error parsing line: {e}")
                continue
    
    if not parsed_data:
        return False, f"Failed to parse any data lines for {scale_type} scale"
    
    return True, pd.DataFrame(parsed_data)


def extract_single_file_distribution(content, sections, filename):
    """
    Extract distribution data from a single file.
    
    Parameters:
    content (str): File content
    sections (dict): Dictionary with section positions
    filename (str): Filename for reference
    
    Returns:
    tuple: (success_flag, result)
        - If successful: (True, combined_df)
        - If failed: (False, error_message)
    """
    # Extract linear data
    lin_start = sections['linear_start']
    log_start = sections['logarithmic_start']
    
    success, lin_result = extract_data_section(content, lin_start, log_start, is_log_section=False)
    if not success:
        return False, f"Failed to extract linear data from {filename}: {lin_result}"
    
    lin_header, lin_data_lines = lin_result
    success, lin_df = parse_data_lines(lin_header, lin_data_lines, 'linear')
    if not success:
        return False, f"Failed to parse linear data from {filename}: {lin_df}"
    
    # Extract logarithmic data
    success, log_result = extract_data_section(content, log_start, is_log_section=True)
    if not success:
        print(f"Warning: Could not extract logarithmic data from {filename}: {log_result}")
        combined_df = lin_df.copy()
    else:
        log_header, log_data_lines = log_result
        success, log_df = parse_data_lines(log_header, log_data_lines, 'logarithmic')
        if not success:
            print(f"Warning: Could not parse logarithmic data from {filename}: {log_df}")
            combined_df = lin_df.copy()
        else:
            combined_df = pd.concat([lin_df, log_df], ignore_index=True)
    
    # Add filename for tracking
    combined_df['source_file'] = filename
    
    return True, combined_df


def average_replicate_data(dataframes_list, filenames_list):
    """
    Average distribution data across multiple replicates bin-by-bin.
    
    Parameters:
    dataframes_list (list): List of DataFrames from individual files
    filenames_list (list): List of filenames for reference
    
    Returns:
    tuple: (success_flag, result)
        - If successful: (True, averaged_df)
        - If failed: (False, error_message)
    """
    if not dataframes_list:
        return False, "No dataframes to average"
    
    if len(dataframes_list) == 1:
        # Single file case - add SD columns with zeros
        df = dataframes_list[0].copy()
        
        # Add SD columns (all zeros for single file)
        df['number_sd'] = 0.0
        df['concentration_cm-3_sd'] = 0.0
        df['volume_nm^3_sd'] = 0.0
        df['area_nm^2_sd'] = 0.0
        
        # Rename value columns to indicate they're averages
        df.rename(columns={
            'number': 'number_avg',
            'concentration_cm-3': 'concentration_cm-3_avg',
            'volume_nm^3': 'volume_nm^3_avg',
            'area_nm^2': 'area_nm^2_avg'
        }, inplace=True)
        
        # Add replicate info
        df['num_replicates'] = 1
        df['source_files'] = filenames_list[0]
        
        print(f"Single file analysis: {filenames_list[0]}")
        print("Standard deviation values set to 0")
        
        return True, df
    
    # Multiple files case - perform averaging
    print(f"Averaging data from {len(dataframes_list)} files:")
    for filename in filenames_list:
        print(f"  {filename}")
    
    # Process each scale separately
    averaged_dfs = []
    
    for scale in ['linear', 'logarithmic']:
        scale_dfs = []
        for df in dataframes_list:
            scale_data = df[df['scale'] == scale]
            if not scale_data.empty:
                scale_dfs.append(scale_data)
        
        if not scale_dfs:
            continue
        
        # Find common size bins across all files
        all_sizes = []
        for df in scale_dfs:
            all_sizes.extend(df['size_nm'].values)
        unique_sizes = sorted(set(all_sizes))
        
        # Create averaged data for this scale
        averaged_data = []
        
        for size in unique_sizes:
            # Collect values for this size across all files
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
            
            # Calculate averages and standard deviations
            if numbers:  # Only if we have data for this size
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
    
    # Combine linear and logarithmic scales
    final_df = pd.concat(averaged_dfs, ignore_index=True)
    
    # Sort by scale and size
    final_df = final_df.sort_values(['scale', 'size_nm']).reset_index(drop=True)
    
    print(f"Averaging completed:")
    print(f"  Total size bins: {len(final_df)}")
    print(f"  Linear scale bins: {len(final_df[final_df['scale'] == 'linear'])}")
    print(f"  Log scale bins: {len(final_df[final_df['scale'] == 'logarithmic'])}")
    
    return True, final_df


# Execute data extraction and averaging if files were processed in Cell 02
if 'current_files_data' in globals() and current_files_data:
    print("=" * 80)
    print("EXTRACTING AND AVERAGING DISTRIBUTION DATA")
    print("=" * 80)
    
    # Extract data from each file
    individual_dataframes = []
    filenames = []
    extraction_errors = []
    
    for filename, content, sections in current_files_data:
        print(f"Extracting data from: {filename}")
        
        success, result = extract_single_file_distribution(content, sections, filename)
        
        if success:
            individual_dataframes.append(result)
            filenames.append(filename)
            print(f"  ✓ Extracted {len(result)} data points")
        else:
            extraction_errors.append((filename, result))
            print(f"  ✗ Failed: {result}")
    
    if individual_dataframes:
        # Average the data across replicates
        print(f"\nAveraging data from {len(individual_dataframes)} successful file(s)...")
        
        success, averaged_df = average_replicate_data(individual_dataframes, filenames)
        
        if success:
            # Store results for downstream cells (mimicking old format)
            current_distribution_df = averaged_df
            current_file_content = None  # Not applicable for averaged data
            current_file_sections = None  # Not applicable for averaged data
            
            # Add uniqueID based on first filename
            first_filename = filenames[0]
            # Extract base name for uniqueID
            base_name = first_filename.replace('_rawdata.txt', '').replace('.txt', '')
            if base_name.startswith('Data_'):
                base_name = base_name[5:]
            
            # For multiple files, add indication
            if len(filenames) > 1:
                uniqueID = f"{base_name}_avg{len(filenames)}"
            else:
                uniqueID = base_name
            
            current_distribution_df['uniqueID'] = uniqueID
            
            print(f"\nData extraction and averaging completed successfully!")
            print(f"Dataset ID: {uniqueID}")
            print(f"Ready for metadata extraction (Cell 04)")
            
            # Display preview
            print(f"\nAveraged Data Preview:")
            display(current_distribution_df.head())
            
        else:
            print(f"\nERROR: Failed to average data: {averaged_df}")
    
    else:
        print(f"\nERROR: No files could be processed for data extraction.")
        if extraction_errors:
            print("Extraction errors:")
            for filename, error in extraction_errors:
                print(f"  {filename}: {error}")

else:
    print("No processed files found. Please run Cell 02 (file processing) first.")
    print("Make sure to run Cell 00 (file selection) before Cell 02.")

"""
NTA Data Analysis - Automated Modular Metadata System (Cell 04)

This module handles:
1. Automated extraction of ALL metadata fields from multiple files
2. Automatic detection of differences vs. identical values
3. Smart array creation for differing values
4. Detailed conflict reporting for analysis

"""

import os
import re
import json
import ntpath  
from datetime import date


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
        # Additional patterns for more comprehensive extraction
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
        
        # Try to find AVI file size if file exists
        if 'CONFIG' in globals():
            directory = CONFIG.get('directory', '')
            full_avi_path = os.path.join(directory, avi_filename)
            if os.path.exists(full_avi_path) and os.path.isfile(full_avi_path):
                size_bytes = os.path.getsize(full_avi_path)
                metadata['avi_filesize'] = f"{size_bytes / (1024 * 1024):.2f} MB"
    
    return metadata


def extract_metadata_from_all_files(files_data):
    """
    Extract metadata from all files and organize by filename.
    
    Parameters:
    files_data (list): List of (filename, content, sections) tuples
    
    Returns:
    tuple: (success_flag, result)
        - If successful: (True, all_files_metadata_dict)
        - If failed: (False, error_message)
    """
    if not files_data:
        return False, "No files provided for metadata extraction"
    
    all_files_metadata = {}
    
    print(f"Extracting ALL metadata fields from {len(files_data)} file(s):")
    
    for filename, content, sections in files_data:
        print(f"  Processing: {filename}")
        
        metadata = extract_all_metadata_fields(content, filename)
        all_files_metadata[filename] = metadata
        
        print(f"    ✓ Extracted {len(metadata)} metadata fields")
    
    print(f"\nExtracted metadata from {len(all_files_metadata)} files")
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
    
    print(f"\nAnalyzing {len(all_field_names)} unique metadata fields:")
    
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
    
    # Report concise summary
    print(f"  Identical fields: {len(identical_fields)}")
    print(f"  Different fields: {len(different_fields)}")
    
    # Only show alerts, not all differences
    quality_issues = []
    for field_name, values in different_fields.items():
        if field_name == 'particle_drift_check_result':
            unique_values = list(set(values))
            if len(unique_values) > 1:
                quality_issues.append(f"QC Alert - {field_name}: {values}")
    
    if quality_issues:
        print(f"\nQuality Control Issues Detected:")
        for issue in quality_issues:
            print(f"  {issue}")
    
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
    import numpy as np
    
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
    
    # Instrument-determined fields - keep as arrays with note
    instrument_determined_fields = [
        'median_number_d50', 'median_concentration_d50', 'median_volume_d50'
    ]
    
    # Quality control fields - flag differences
    quality_control_fields = [
        'particle_drift_check_result', 'cell_check_result'
    ]
    
    notes = ""
    
    if field_name in file_specific_fields:
        # Keep as JSON array
        return json.dumps(values), "file_specific"
    
    elif field_name in text_first_fields:
        # Use first value only
        return values[0], f"using_first_of_{len(values)}"
    
    elif field_name in sum_fields:
        # Sum all values
        try:
            numeric_values = [float(v) for v in values]
            total = sum(numeric_values)
            return f"{total:.0f}", f"sum_of_{len(values)}_measurements"
        except ValueError:
            return json.dumps(values), "non_numeric_sum_field"
    
    elif field_name in instrument_determined_fields:
        # Keep as array with explanatory note
        return json.dumps(values), "instrument_determined_per_measurement"
    
    elif field_name in quality_control_fields:
        # Check if all values are the same
        unique_values = list(set(values))
        if len(unique_values) == 1:
            return values[0], "qc_consistent"
        else:
            return json.dumps(values), f"QC_ALERT_inconsistent_values"
    
    else:
        # Try to calculate mean ± SD for numeric fields
        try:
            # Special handling for file sizes with units
            if field_name == 'avi_filesize':
                # Extract numeric values from strings like "219.24 MB"
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
            
            # Regular numeric processing for other fields
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
    
    Parameters:
    all_files_metadata (dict): All files' metadata
    identical_fields (dict): Fields with identical values
    different_fields (dict): Fields with different values
    config (dict): Configuration dictionary
    
    Returns:
    dict: Essential standardized metadata dictionary with ordered fields
    """
    filenames = list(all_files_metadata.keys())
    num_files = len(filenames)
    
    # Start with ordered metadata structure
    metadata = {}
    processing_notes = {}
    
    # Create uniqueID with averaging suffix
    if 'uniqueID' in identical_fields:
        base_id = identical_fields['uniqueID']
    else:
        # Take from first file if not identical
        base_id = list(all_files_metadata.values())[0].get('uniqueID', 'unknown')
    
    if num_files > 1:
        unique_id = f"{base_id}_avg{num_files}"
    else:
        unique_id = base_id
    
    # Check for custom persistent_ID from config
    if config and "project_metadata" in config:
        project_meta = config["project_metadata"]
        if "persistent_ID" in project_meta:
            custom_id = project_meta["persistent_ID"]
            if num_files > 1:
                unique_id = f"{custom_id}_avg{num_files}"
            else:
                unique_id = custom_id
            print(f"Using custom persistent ID: {unique_id}")
    else:
        project_meta = {}
    
    # SECTION 1: CORE IDENTIFICATION (first in file)
    metadata['experimenter'] = project_meta.get('experimenter', 'SH/HB')
    metadata['location'] = project_meta.get('location', 'HD_MPImF_CBP_R0.106')
    metadata['project'] = project_meta.get('project', 'LEAF')
    metadata['meta_version'] = project_meta.get('meta_version', 'v02')
    metadata['pi'] = project_meta.get('pi', 'HB')
    metadata['funding'] = project_meta.get('funding', 'MPG')
    metadata['persistentID'] = unique_id
    metadata['data_collection_method'] = project_meta.get('data_collection_method', 'NTA')
    metadata['nta_instrument'] = 'ZetaView'
    
    # NTA software version
    if 'analyze' in identical_fields:
        metadata['nta_software'] = f"ZetaView {identical_fields['analyze']}"
    elif 'analyze' in different_fields:
        metadata['nta_software'] = f"ZetaView {different_fields['analyze'][0]}"
    else:
        metadata['nta_software'] = 'ZetaView'
    
    metadata['nta_processed_file'] = f"Data_{unique_id}_PSD.txt"
    
    # Sample info
    if 'sample' in identical_fields:
        metadata['sample'] = identical_fields['sample']
    elif 'sample' in different_fields:
        metadata['sample'] = different_fields['sample'][0]  # Use first sample name
    
    # SECTION 2: MULTI-FILE INFO
    metadata['num_replicates'] = num_files
    metadata['source_files'] = json.dumps(filenames)
    
    # SECTION 3: MEASUREMENT PARAMETERS
    # Define essential fields to save (others will be kept in memory only)
    essential_fields = {
        # Key measurement parameters
        'date', 'temperature', 'ph', 'dilution', 'laser_wavelength', 'electrolyte',
        'positions', 'cycles', 'fps', 
        # Quality control
        'particle_drift_check_result', 'cell_check_result',
        # Key results (not instrument D50s - we calculate better ones later)
        'average_number_of_particles', 'number_of_traces', 'detected_particles',
        'conductivity', 'scattering_intensity', 'viscosity',
        # File info (simplified)
        'avi_filesize'
    }
    
    # Process identical fields (only essential ones)
    for field_name, value in identical_fields.items():
        if field_name in essential_fields:
            # Add nta_ prefix for measurement-related fields
            if field_name in ['temperature', 'ph', 'dilution', 'laser_wavelength', 'positions', 
                             'cycles', 'fps', 'average_number_of_particles', 'number_of_traces', 
                             'detected_particles', 'particle_drift_check_result', 'cell_check_result',
                             'conductivity', 'scattering_intensity', 'viscosity', 'avi_filesize']:
                # Special naming for summed fields
                if field_name in ['number_of_traces', 'detected_particles']:
                    metadata[f'nta_{field_name}_sum'] = value
                else:
                    metadata[f'nta_{field_name}'] = value
            else:
                metadata[field_name] = value
    
    # Process different fields with smart formatting (only essential ones)
    quality_alerts = []
    high_variation_fields = []
    
    for field_name, values in different_fields.items():
        if field_name in essential_fields:
            # Apply smart formatting
            formatted_value, notes = smart_format_field(field_name, values)
            
            # Track quality issues
            if "QC_ALERT" in notes:
                quality_alerts.append(f"{field_name}: {formatted_value}")
            if "HIGH_VARIATION" in notes:
                high_variation_fields.append(f"{field_name}: {notes}")
            
            # Add nta_ prefix for measurement-related fields
            if field_name in ['temperature', 'ph', 'dilution', 'laser_wavelength', 'positions', 
                             'cycles', 'fps', 'average_number_of_particles', 'number_of_traces', 
                             'detected_particles', 'particle_drift_check_result', 'cell_check_result',
                             'conductivity', 'scattering_intensity', 'viscosity', 'avi_filesize']:
                # Special naming for summed fields
                if field_name in ['number_of_traces', 'detected_particles']:
                    metadata[f'nta_{field_name}_sum'] = formatted_value
                else:
                    metadata[f'nta_{field_name}'] = formatted_value
            else:
                metadata[field_name] = formatted_value
            
            # Store processing note (in memory only)
            processing_notes[field_name] = notes
    
    # SECTION 4: ADDITIONAL REFERENCES
    # Note: nta_plot_file will be added later when plots are actually generated
    metadata['python_analysis'] = str(date.today())
    
    # SECTION 5: QUALITY ALERTS (only if present)
    if quality_alerts:
        metadata['quality_control_alerts'] = json.dumps(quality_alerts)
        print(f"\n⚠ QUALITY CONTROL ALERTS:")
        for alert in quality_alerts:
            print(f"  {alert}")
    
    if high_variation_fields:
        metadata['high_variation_fields'] = json.dumps(high_variation_fields)
        print(f"\n⚠ HIGH VARIATION DETECTED:")
        for field in high_variation_fields:
            print(f"  {field}")
    
    # Store detailed analysis in memory only (not in saved file)
    global current_processing_notes, current_all_fields_metadata, current_original_differences
    current_processing_notes = processing_notes
    current_all_fields_metadata = all_files_metadata  # Full detailed metadata
    current_original_differences = different_fields    # All original differences
    
    return metadata


def save_metadata_file(metadata, output_dir=None, config=None):
    """
    Save metadata to a file in the specified output directory.
    """
    # Determine output directory
    if output_dir is None:
        if config is not None and "directory" in config:
            base_dir = config["directory"]
            output_dir = os.path.join(base_dir, "metadata")
        else:
            output_dir = os.path.join(os.getcwd(), "metadata")
    
    # Ensure directory exists
    try:
        os.makedirs(output_dir, exist_ok=True)
    except Exception as e:
        return False, f"Failed to create metadata directory: {str(e)}"
    
    # Create filepath
    unique_id = metadata.get('persistentID', 'unknown')
    metadata_path = os.path.join(output_dir, f"Data_{unique_id}_metadata.txt")
    
    try:
        # Write metadata
        with open(metadata_path, 'w') as f:
            for key, value in metadata.items():
                f.write(f"{key}\t{value}\t\n")
        
        return True, metadata_path
    except Exception as e:
        return False, f"Failed to write metadata file: {str(e)}"


# Execute automated metadata extraction if files were processed
if 'current_files_data' in globals() and current_files_data:
    print("=" * 80)
    print("AUTOMATED METADATA EXTRACTION AND ANALYSIS")
    print("=" * 80)
    
    # Step 1: Extract all metadata from all files
    success, all_files_metadata = extract_metadata_from_all_files(current_files_data)
    
    if success:
        # Step 2: Analyze which fields are identical vs. different
        identical_fields, different_fields, field_analysis = analyze_field_differences(all_files_metadata)
        
        # Step 3: Create standardized metadata with automated rules
        standard_metadata = create_automated_metadata(
            all_files_metadata, identical_fields, different_fields, CONFIG
        )
        
        # Step 4: Save metadata file
        success, metadata_path = save_metadata_file(standard_metadata, config=CONFIG)
        
        if success:
            print(f"\nAutomated metadata processing completed successfully!")
            print(f"Saved comprehensive metadata to: {metadata_path}")
            
            # Store for downstream cells
            current_metadata = standard_metadata
            
            print(f"\nMetadata Summary:")
            print(f"  persistentID: {standard_metadata['persistentID']}")
            print(f"  num_replicates: {standard_metadata['num_replicates']}")
            print(f"  total_fields_extracted: {len(standard_metadata)}")
            
            # Store detailed analysis for inspection (but don't clutter output)
            current_field_analysis = field_analysis
            current_identical_fields = identical_fields
            current_different_fields = different_fields
            # current_processing_notes, current_all_fields_metadata, and current_original_differences
            # are stored globally by create_automated_metadata()
            
            print(f"\nDetailed analysis available in memory:")
            print(f"  current_all_fields_metadata (full extraction from all files)")
            print(f"  current_original_differences (all field differences)")
            print(f"  current_processing_notes (how each field was processed)")
            
        else:
            print(f"ERROR: Failed to save metadata: {metadata_path}")
    
    else:
        print(f"ERROR: Failed to extract metadata: {all_files_metadata}")

else:
    print("No processed files found. Please run Cell 02 and Cell 03 first.")

"""
NTA Data Analysis - Core Calculations with Uncertainty Propagation (Cell 05)

This module provides core calculation functions for NTA data analysis with rigorous
uncertainty propagation through averaging of replicates:

1. Dilution correction for all measurements
2. Normalization of particle distributions with uncertainty propagation
3. Calculation of cumulative distributions with proper error propagation
4. Total metrics calculation with uncertainty combination

These functions provide the mathematical foundation for analyzing averaged
number, volume, and surface area distributions from replicate NTA measurements.
"""

import numpy as np
import pandas as pd
from scipy import integrate


def apply_dilution_correction_with_uncertainty(df, metadata=None, manual_dilution=None):
    """
    Apply dilution correction to all measured values with uncertainty propagation.
    
    For a dilution factor D, the actual sample concentration = measured × D
    Uncertainty propagation: σ_actual = σ_measured × D
    
    Parameters:
    df (DataFrame): DataFrame containing averaged distribution data
    metadata (dict): Metadata dictionary that may contain dilution info
    manual_dilution (float): Manual dilution factor override (optional)
    
    Returns:
    tuple: (success_flag, result)
        - If successful: (True, updated_df)
        - If failed: (False, error_message)
    """
    # Create a copy to avoid modifying the original dataframe
    updated_df = df.copy()
    
    # Determine dilution factor
    dilution_factor = 1.0  # Default: no dilution
    dilution_source = "default (no dilution)"
    
    if manual_dilution is not None:
        try:
            dilution_factor = float(manual_dilution)
            dilution_source = "manually specified"
        except (ValueError, TypeError):
            return False, f"Invalid manual dilution factor: {manual_dilution}"
    
    elif metadata is not None:
        # Try nta_dilution field
        if 'nta_dilution' in metadata:
            try:
                dilution_string = metadata['nta_dilution']
                # Handle format like "10.0" or other formats
                dilution_factor = float(dilution_string.split('±')[0].strip()) if '±' in dilution_string else float(dilution_string)
                dilution_source = "metadata (nta_dilution)"
            except (ValueError, TypeError):
                print("Warning: Could not parse dilution factor from metadata, using default (1.0)")
    
    print(f"Applying dilution correction: factor = {dilution_factor} (source: {dilution_source})")
    
    # Apply dilution correction to concentration and rename
    if 'concentration_cm-3_avg' in updated_df.columns:
        updated_df['particles_per_mL_avg'] = updated_df['concentration_cm-3_avg'] * dilution_factor
        if 'concentration_cm-3_sd' in updated_df.columns:
            updated_df['particles_per_mL_sd'] = updated_df['concentration_cm-3_sd'] * dilution_factor
        
        # Remove old concentration columns
        updated_df = updated_df.drop(['concentration_cm-3_avg', 'concentration_cm-3_sd'], axis=1, errors='ignore')
    else:
        return False, "Missing concentration_cm-3_avg column for dilution correction"
    
    # Apply dilution correction to volume and rename to per_mL
    if 'volume_nm^3_avg' in updated_df.columns:
        updated_df['volume_nm^3_per_mL_avg'] = updated_df['volume_nm^3_avg'] * dilution_factor
        if 'volume_nm^3_sd' in updated_df.columns:
            updated_df['volume_nm^3_per_mL_sd'] = updated_df['volume_nm^3_sd'] * dilution_factor
        
        # Remove old volume columns
        updated_df = updated_df.drop(['volume_nm^3_avg', 'volume_nm^3_sd'], axis=1, errors='ignore')
    
    # Apply dilution correction to area and rename to per_mL
    if 'area_nm^2_avg' in updated_df.columns:
        updated_df['area_nm^2_per_mL_avg'] = updated_df['area_nm^2_avg'] * dilution_factor
        if 'area_nm^2_sd' in updated_df.columns:
            updated_df['area_nm^2_per_mL_sd'] = updated_df['area_nm^2_sd'] * dilution_factor
        
        # Remove old area columns
        updated_df = updated_df.drop(['area_nm^2_avg', 'area_nm^2_sd'], axis=1, errors='ignore')
    
    return True, updated_df


def normalize_distributions_with_uncertainty(df, size_column='size_nm'):
    """
    Normalize particle distributions by area under the curve with uncertainty propagation.
    
    This creates normalized number distributions from the averaged number data.
    
    Parameters:
    df (DataFrame): DataFrame containing averaged particle distribution data
    size_column (str): Name of the column containing size values
    
    Returns:
    DataFrame: Updated dataframe with normalized columns and uncertainties
    """
    # Create a copy to avoid modifying the original dataframe
    normalized_df = df.copy()
    
    # Process each scale separately
    for scale in normalized_df['scale'].unique():
        scale_mask = normalized_df['scale'] == scale
        scale_data = normalized_df[scale_mask].copy()
        
        if scale_data.empty or 'number_avg' not in scale_data.columns:
            continue
        
        # Sort by size for correct integration
        scale_data = scale_data.sort_values(size_column)
        
        # Calculate area under the curve using trapezoidal rule
        sizes = scale_data[size_column].values
        numbers_avg = scale_data['number_avg'].values
        
        if len(sizes) < 2:
            continue
        
        # Calculate area for normalization
        area_avg = np.trapz(numbers_avg, sizes)
        
        if area_avg > 0:
            # Normalize the averages
            normalized_df.loc[scale_mask, 'number_normalized_avg'] = numbers_avg / area_avg
            
            # Normalize the standard deviations (uncertainty propagation)
            if 'number_sd' in scale_data.columns:
                numbers_sd = scale_data['number_sd'].values
                normalized_df.loc[scale_mask, 'number_normalized_sd'] = numbers_sd / area_avg
            else:
                normalized_df.loc[scale_mask, 'number_normalized_sd'] = 0.0
        else:
            # If area is zero, set normalized values to zero
            normalized_df.loc[scale_mask, 'number_normalized_avg'] = 0.0
            normalized_df.loc[scale_mask, 'number_normalized_sd'] = 0.0
    
    return normalized_df


def calculate_cumulative_distributions_with_uncertainty(df, scale_column='scale'):
    """
    Calculate cumulative distributions with proper uncertainty propagation.
    
    For independent uncertainties, cumulative uncertainties are calculated as:
    σ_cumsum[j] = √(Σ(i=0 to j) σ[i]²)
    
    Parameters:
    df (DataFrame): DataFrame containing particle distribution data with uncertainties
    scale_column (str): Name of the column distinguishing scale types
    
    Returns:
    DataFrame: Updated dataframe with cumulative distribution columns and uncertainties
    """
    # Create a copy to avoid modifying the original
    result_df = df.copy()
    
    # Process each scale type separately
    for scale in result_df[scale_column].unique():
        # Filter for current scale
        scale_mask = result_df[scale_column] == scale
        scale_indices = result_df[scale_mask].sort_values('size_nm').index
        
        if len(scale_indices) == 0:
            continue
        
        # 1. Normalized number distribution (shape analysis)
        if 'number_normalized_avg' in result_df.columns:
            # Calculate cumulative sum
            cumsum_avg = result_df.loc[scale_indices, 'number_normalized_avg'].cumsum()
            
            # Normalize to 0-1 range (should be close to 1 already due to normalization)
            if cumsum_avg.iloc[-1] > 0:
                result_df.loc[scale_indices, 'number_normalized_cumsum_avg'] = cumsum_avg / cumsum_avg.iloc[-1]
            else:
                result_df.loc[scale_indices, 'number_normalized_cumsum_avg'] = 0
            
            # Calculate uncertainty in cumulative sum using error propagation
            if 'number_normalized_sd' in result_df.columns:
                # For cumulative sum: σ_cumsum[j] = √(Σ(i=0 to j) σ[i]²)
                normalized_var_cumsum = (result_df.loc[scale_indices, 'number_normalized_sd'] ** 2).cumsum()
                cumsum_sd = np.sqrt(normalized_var_cumsum)
                
                # Normalize the uncertainty as well
                if cumsum_avg.iloc[-1] > 0:
                    result_df.loc[scale_indices, 'number_normalized_cumsum_sd'] = cumsum_sd / cumsum_avg.iloc[-1]
                else:
                    result_df.loc[scale_indices, 'number_normalized_cumsum_sd'] = 0
        
        # 2. Absolute volume distribution
        if 'volume_nm^3_per_mL_avg' in result_df.columns:
            # Calculate cumulative sum for averages
            cumsum_avg = result_df.loc[scale_indices, 'volume_nm^3_per_mL_avg'].cumsum()
            result_df.loc[scale_indices, 'volume_nm^3_per_mL_cumsum_avg'] = cumsum_avg
            
            # Calculate uncertainty in cumulative sum
            if 'volume_nm^3_per_mL_sd' in result_df.columns:
                volume_var_cumsum = (result_df.loc[scale_indices, 'volume_nm^3_per_mL_sd'] ** 2).cumsum()
                result_df.loc[scale_indices, 'volume_nm^3_per_mL_cumsum_sd'] = np.sqrt(volume_var_cumsum)
        
        # 3. Absolute surface area distribution
        if 'area_nm^2_per_mL_avg' in result_df.columns:
            # Calculate cumulative sum for averages
            cumsum_avg = result_df.loc[scale_indices, 'area_nm^2_per_mL_avg'].cumsum()
            result_df.loc[scale_indices, 'area_nm^2_per_mL_cumsum_avg'] = cumsum_avg
            
            # Calculate uncertainty in cumulative sum
            if 'area_nm^2_per_mL_sd' in result_df.columns:
                area_var_cumsum = (result_df.loc[scale_indices, 'area_nm^2_per_mL_sd'] ** 2).cumsum()
                result_df.loc[scale_indices, 'area_nm^2_per_mL_cumsum_sd'] = np.sqrt(area_var_cumsum)
                
    return result_df


def calculate_total_metrics_with_uncertainty(df, scale_column='scale'):
    """
    Calculate total metrics for each scale with proper uncertainty propagation.
    
    For totals across bins, uncertainties are combined as: σ_total = √(Σ σ_i²)
    For derived metrics, simple calculations are used without uncertainty propagation.
    
    Parameters:
    df (DataFrame): DataFrame containing particle distribution data with uncertainties
    scale_column (str): Name of the column distinguishing scale types
    
    Returns:
    dict: Dictionary with total metrics and uncertainties for each scale
    """
    # Initialize results structure
    results = {}
    
    # Process each scale type separately
    for scale in df[scale_column].unique():
        # Filter for current scale and sort by size
        scale_df = df[df[scale_column] == scale].sort_values('size_nm')
        
        # Initialize metrics for this scale
        scale_metrics = {}
        
        # Only calculate if we have data for this scale
        if not scale_df.empty:
            
            # 1. Total particles per mL
            if 'particles_per_mL_avg' in scale_df.columns:
                total_particles_avg = scale_df['particles_per_mL_avg'].sum()
                scale_metrics['total_particles_per_mL_avg'] = total_particles_avg
                scale_metrics['total_particles_per_uL_avg'] = total_particles_avg / 1000
                
                # Calculate uncertainty: σ_total = √(Σ σ_i²)
                if 'particles_per_mL_sd' in scale_df.columns:
                    total_particles_sd = np.sqrt((scale_df['particles_per_mL_sd'] ** 2).sum())
                    scale_metrics['total_particles_per_mL_sd'] = total_particles_sd
                    scale_metrics['total_particles_per_uL_sd'] = total_particles_sd / 1000
            
            # 2. Total volume per mL
            if 'volume_nm^3_per_mL_avg' in scale_df.columns:
                total_volume_avg = scale_df['volume_nm^3_per_mL_avg'].sum()
                scale_metrics['total_volume_nm^3_per_mL_avg'] = total_volume_avg
                scale_metrics['total_volume_um^3_per_mL_avg'] = total_volume_avg / 1e9  # nm³ to μm³
                scale_metrics['total_volume_uL_per_mL_avg'] = total_volume_avg / 1e18  # nm³ to μL
                scale_metrics['volume_percentage_avg'] = (total_volume_avg / 1e18) * 0.1  # percentage
                
                # Calculate uncertainty
                if 'volume_nm^3_per_mL_sd' in scale_df.columns:
                    total_volume_sd = np.sqrt((scale_df['volume_nm^3_per_mL_sd'] ** 2).sum())
                    scale_metrics['total_volume_nm^3_per_mL_sd'] = total_volume_sd
                    scale_metrics['total_volume_um^3_per_mL_sd'] = total_volume_sd / 1e9
                    scale_metrics['total_volume_uL_per_mL_sd'] = total_volume_sd / 1e18
                    scale_metrics['volume_percentage_sd'] = (total_volume_sd / 1e18) * 0.1
            
            # 3. Total surface area per mL
            if 'area_nm^2_per_mL_avg' in scale_df.columns:
                total_area_avg = scale_df['area_nm^2_per_mL_avg'].sum()
                scale_metrics['total_surface_area_nm^2_per_mL_avg'] = total_area_avg
                scale_metrics['total_surface_area_um^2_per_mL_avg'] = total_area_avg / 1e6  # nm² to μm²
                scale_metrics['total_surface_area_cm^2_per_mL_avg'] = total_area_avg / 1e14  # nm² to cm²
                
                # Calculate uncertainty
                if 'area_nm^2_per_mL_sd' in scale_df.columns:
                    total_area_sd = np.sqrt((scale_df['area_nm^2_per_mL_sd'] ** 2).sum())
                    scale_metrics['total_surface_area_nm^2_per_mL_sd'] = total_area_sd
                    scale_metrics['total_surface_area_um^2_per_mL_sd'] = total_area_sd / 1e6
                    scale_metrics['total_surface_area_cm^2_per_mL_sd'] = total_area_sd / 1e14
                
                # 4. Specific surface area (derived metric, no uncertainty propagation for now)
                if ('total_volume_nm^3_per_mL_avg' in scale_metrics and 
                    scale_metrics['total_volume_nm^3_per_mL_avg'] > 0):
                    # Surface area (nm²) / volume (nm³) = 1/nm
                    ssa_1_per_nm = total_area_avg / scale_metrics['total_volume_nm^3_per_mL_avg']
                    # Convert to m²/cm³ (standard unit)
                    scale_metrics['specific_surface_area_m^2_per_cm^3_avg'] = ssa_1_per_nm * 10
        
        # Store metrics for this scale
        results[scale] = scale_metrics
    
    return results


def add_metrics_to_metadata_with_uncertainty(metadata, metrics, scale='linear'):
    """
    Add key metrics with uncertainties to the metadata dictionary.
    
    Parameters:
    metadata (dict): Current metadata dictionary
    metrics (dict): Dictionary of calculated metrics with uncertainties
    scale (str): Which scale's metrics to use ('linear' or 'logarithmic')
    
    Returns:
    dict: Updated metadata dictionary
    """
    # Create a copy to avoid modifying the original
    updated_metadata = metadata.copy()
    
    # Use linear scale metrics by default, but check if available
    if scale not in metrics:
        # Fall back to any available scale
        if metrics:
            scale = list(metrics.keys())[0]
        else:
            return updated_metadata  # Return unchanged if no metrics available
    
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
        updated_metadata['nta_volume_percentage'] = f"{avg_val:.6f} ± {sd_val:.6f}"
    
    if 'total_surface_area_cm^2_per_mL_avg' in scale_metrics:
        avg_val = scale_metrics['total_surface_area_cm^2_per_mL_avg']
        sd_val = scale_metrics.get('total_surface_area_cm^2_per_mL_sd', 0)
        updated_metadata['nta_total_surface_area_cm^2_per_mL'] = f"{avg_val:.4E} ± {sd_val:.4E}"
    
    if 'specific_surface_area_m^2_per_cm^3_avg' in scale_metrics:
        avg_val = scale_metrics['specific_surface_area_m^2_per_cm^3_avg']
        # No uncertainty for derived metrics yet
        updated_metadata['nta_specific_surface_area_m^2_per_cm^3'] = f"{avg_val:.2f}"
    
    # Add a scale indicator and number of replicates
    updated_metadata['nta_metrics_scale'] = scale
    if 'num_replicates' in updated_metadata:
        updated_metadata['nta_metrics_replicates'] = updated_metadata['num_replicates']
    
    return updated_metadata



# ============================================================================
# STREAMLIT ORCHESTRATOR - Calls cells in sequence
# ============================================================================

def run_analysis_pipeline(filepaths, config=None):
    """
    Run the complete analysis pipeline (Cells 02-05) in sequence.
    
    This replaces running cells manually in Jupyter.
    
    Parameters:
    filepaths (list): List of file paths to analyze
    config (dict): Configuration dictionary (optional)
    
    Returns:
    dict: Results dictionary with all outputs
    """
    
    if config is None:
        config = CONFIG
    
    # CELL 02: Process files
    print("\n" + "="*80)
    print("CELL 02: FILE I/O")
    print("="*80)
    successful_files, failed_files = process_multiple_file_contents(filepaths)
    
    if not successful_files:
        raise Exception("No files could be processed")
    
    # CELL 03: Extract and average distributions
    print("\n" + "="*80)
    print("CELL 03: DATA EXTRACTION & AVERAGING")
    print("="*80)
    distributions = []
    for filename, content, sections in successful_files:
        success, dist = extract_single_file_distribution(content, sections, filename)
        if success:
            distributions.append(dist)
    
    if not distributions:
        raise Exception("Could not extract distribution data")
    
    filenames = [f[0] for f in successful_files]
    success, avg_dist = average_replicate_data(distributions, filenames)
    if not success:
        raise Exception(f"Failed to average: {avg_dist}")
    
    # CELL 04: Metadata
    print("\n" + "="*80)
    print("CELL 04: METADATA EXTRACTION")
    print("="*80)
    success, all_metadata = extract_metadata_from_all_files(successful_files)
    if not success:
        raise Exception(f"Failed to extract metadata: {all_metadata}")
    
    identical, different, analysis = analyze_field_differences(all_metadata)
    metadata = create_automated_metadata(all_metadata, identical, different, config)
    
    # CELL 05: Core calculations
    print("\n" + "="*80)
    print("CELL 05: CORE CALCULATIONS WITH UNCERTAINTY PROPAGATION")
    print("="*80)
    
    # Step 1: Dilution correction
    print("\n1. Dilution correction...")
    success, dilution_df = apply_dilution_correction_with_uncertainty(
        avg_dist, metadata=metadata
    )
    if not success:
        print(f"Warning: {dilution_df}")
        processed_dist = avg_dist
    else:
        processed_dist = dilution_df
    
    # Step 2: Normalization
    print("2. Normalization...")
    normalized_df = normalize_distributions_with_uncertainty(processed_dist)
    processed_dist = normalized_df
    
    # Step 3: Cumulative distributions
    print("3. Cumulative distributions...")
    cumsum_df = calculate_cumulative_distributions_with_uncertainty(processed_dist)
    processed_dist = cumsum_df
    
    # Step 4: Total metrics
    print("4. Total metrics...")
    total_metrics = calculate_total_metrics_with_uncertainty(processed_dist)
    
    # Step 5: Add to metadata
    metadata = add_metrics_to_metadata_with_uncertainty(metadata, total_metrics, scale='linear')
    
    print("\n" + "="*80)
    print("ANALYSIS COMPLETE")
    print("="*80)
    
    return {
        'distribution': processed_dist,
        'metadata': metadata,
        'total_metrics': total_metrics,
        'identical_fields': identical,
        'different_fields': different,
        'field_analysis': analysis,
        'all_file_metadata': all_metadata,
        'failed_files': failed_files,
        'num_replicates': len(successful_files)
    }
