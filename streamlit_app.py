"""
NTA Data Analysis - Streamlit Web Application
Updated to use combined analyzer with Cells 01-06
Includes metadata discrepancy detection, warnings, metrics, and statistics
"""

import streamlit as st
import pandas as pd
import numpy as np
import io
import tempfile
import os
from nta_analyzer_cells_01_04 import NTAAnalyzer, CONFIG

# Set page config
st.set_page_config(
    page_title="NTA Analysis",
    page_icon="beaker",
    layout="wide",
    initial_sidebar_state="expanded"
)

st.title("NTA Particle Size Analysis")
st.markdown("---")

# Initialize session state
if 'analyzer' not in st.session_state:
    st.session_state.analyzer = None
if 'results' not in st.session_state:
    st.session_state.results = None

# Sidebar configuration
with st.sidebar:
    st.header("Configuration")
    
    st.subheader("Project Information")
    experimenter = st.text_input(
        "Experimenter",
        value=CONFIG["project_metadata"]["experimenter"],
        help="Your initials or name"
    )
    location = st.text_input(
        "Lab Location",
        value=CONFIG["project_metadata"]["location"],
        help="Where the analysis is performed"
    )
    project = st.text_input(
        "Project Name",
        value=CONFIG["project_metadata"]["project"],
        help="Your project name"
    )
    pi = st.text_input(
        "Principal Investigator",
        value=CONFIG["project_metadata"].get("pi", "Your_PI_Initials"),
        help="PI initials or name"
    )
    funding = st.text_input(
        "Funding Source",
        value=CONFIG["project_metadata"].get("funding", "none"),
        help="Funding source (optional, default: none)",
        placeholder="none"
    )
    
    st.subheader("Data Identification")
    manual_persistent_id = st.text_input(
        "Manual Persistent ID (optional)",
        value="",
        help="Override auto-generated persistentID. Leave empty to auto-generate from filename",
        placeholder="Leave empty for auto-generation"
    )
    if manual_persistent_id:
        st.info(f"Will use: {manual_persistent_id}")
    
    # Update CONFIG with user values
    CONFIG["project_metadata"]["experimenter"] = experimenter
    CONFIG["project_metadata"]["location"] = location
    CONFIG["project_metadata"]["project"] = project
    CONFIG["project_metadata"]["pi"] = pi
    CONFIG["project_metadata"]["funding"] = funding if funding else "none"
    
    # Store manual persistent ID if provided
    if manual_persistent_id:
        CONFIG["manual_persistent_id"] = manual_persistent_id
    
    st.info("Configuration changes apply when you click 'Analyze Files'")
    
    # Reset button
    if st.button("Reset All"):
        st.session_state.clear()
        st.rerun()

# Main content
st.header("Upload NTA Files")
st.markdown("Upload one or more NTA data files (.txt format)")

uploaded_files = st.file_uploader(
    "Choose NTA files",
    type="txt",
    accept_multiple_files=True,
    help="Select one or more NTA data files for analysis"
)

if uploaded_files:
    if st.button("Analyze Files", key="analyze_btn", type="primary"):
        with st.spinner("Processing files..."):
            try:
                # Create temp directory
                with tempfile.TemporaryDirectory() as temp_dir:
                    # Save uploaded files to temp directory
                    temp_files = []
                    for uploaded_file in uploaded_files:
                        temp_path = os.path.join(temp_dir, uploaded_file.name)
                        with open(temp_path, 'wb') as f:
                            f.write(uploaded_file.getbuffer())
                        temp_files.append(temp_path)
                    
                    # Run analysis
                    analyzer = NTAAnalyzer(config=CONFIG)
                    results = analyzer.process(temp_files)
                    
                    # Store in session state
                    st.session_state.analyzer = analyzer
                    st.session_state.results = results
                
                st.success("Analysis completed!")
                
            except Exception as e:
                st.error(f"Error during analysis: {str(e)}")
                st.stop()

# Display results if analysis was successful
if st.session_state.results:
    results = st.session_state.results
    
    st.markdown("---")
    st.header("Analysis Results")
    
    # Tabs for different views
    tab0, tab1, tab2, tab3, tab4, tab5 = st.tabs([
        "Overview",
        "Distribution Data",
        "Metadata",
        "Warnings",
        "Metrics",
        "Download"
    ])
    
    # TAB 0: Overview (Key Metrics)
    with tab0:
        st.subheader("Analysis Overview")
        
        # Key metrics in columns
        col1, col2, col3 = st.columns(3)
        
        with col1:
            particles = results['metadata'].get('nta_total_particles_per_mL', 'N/A')
            st.metric("Total Particles/mL", particles)
        
        with col2:
            volume = results['metadata'].get('nta_total_volume_uL_per_mL', 'N/A')
            st.metric("Total Volume (uL/mL)", volume)
        
        with col3:
            ssa = results['metadata'].get('nta_specific_surface_area_m^2_per_cm^3', 'N/A')
            st.metric("Specific Surface Area", ssa)
        
        # D50 and Span - most important size metrics
        st.subheader("Size Distribution")
        
        col1, col2 = st.columns(2)
        
        with col1:
            d50 = results['metadata'].get('nta_linear_number_d50', 'N/A')
            st.info(f"**D50 (Number-weighted):** {d50}")
        
        with col2:
            span = results['metadata'].get('nta_linear_number_span', 'N/A')
            st.info(f"**Span (Distribution width):** {span}")
        
        # Additional info
        st.subheader("Analysis Information")
        
        col1, col2, col3 = st.columns(3)
        
        with col1:
            st.metric("Files Analyzed", results['num_replicates'])
        
        with col2:
            date = results['metadata'].get('nta_python_analysis', 'N/A')
            st.metric("Analysis Date", date)
        
        with col3:
            scale = results['metadata'].get('nta_metrics_scale', 'N/A')
            st.metric("Scale", scale)
    
    # TAB 1: Distribution Data
    with tab1:
        st.subheader("Distribution Data (PSD)")
        
        # Show metrics
        col1, col2, col3 = st.columns(3)
        with col1:
            st.metric("Files", results['num_replicates'])
        with col2:
            st.metric("Traces", results['metadata'].get('nta_number_of_traces_sum', 'N/A'))
        with col3:
            if results['quality_alerts']:
                st.metric("Alert", "High drift")
            elif results['high_variation_fields']:
                st.metric("Status", "High variation")
            else:
                st.metric("Status", "Good")
        
        # Get distribution data
        dist = results['distribution']
        linear_data = dist[dist['scale'] == 'linear'].copy()
        log_data = dist[dist['scale'] == 'logarithmic'].copy()
        
        # Create tabs for linear and logarithmic
        tab_linear, tab_log = st.tabs(["Linear Scale", "Logarithmic Scale"])
        
        with tab_linear:
            st.write(f"Linear scale: {len(linear_data)} size bins")
            
            # Column information
            with st.expander("Column Information (Linear Scale)"):
                st.write(f"Total columns: {len(linear_data.columns)}")
                col_groups = {
                    'Size': ['size_nm'],
                    'Raw Counts': ['number_avg', 'number_sd'],
                    'Metadata': ['num_replicates', 'source_files'],
                    'Concentration (dilution-corrected)': ['particles_per_mL_avg', 'particles_per_mL_sd'],
                    'Volume (dilution-corrected)': ['volume_nm^3_per_mL_avg', 'volume_nm^3_per_mL_sd'],
                    'Surface Area (dilution-corrected)': ['area_nm^2_per_mL_avg', 'area_nm^2_per_mL_sd'],
                    'Normalized Number [CRITICAL]': ['number_normalized_avg', 'number_normalized_sd'],
                    'Cumulative Normalized [CRITICAL]': ['number_normalized_cumsum_avg', 'number_normalized_cumsum_sd'],
                    'Cumulative Volume': ['volume_nm^3_per_mL_cumsum_avg', 'volume_nm^3_per_mL_cumsum_sd'],
                    'Cumulative Surface Area': ['area_nm^2_per_mL_cumsum_avg', 'area_nm^2_per_mL_cumsum_sd']
                }
                
                for group_name, cols in col_groups.items():
                    available = [c for c in cols if c in linear_data.columns]
                    missing = [c for c in cols if c not in linear_data.columns]
                    
                    status = "OK" if not missing else "MISSING"
                    st.write(f"**{group_name}** [{status}]")
                    for col in available:
                        st.write(f"  - {col}")
                    if missing:
                        st.warning(f"  Missing: {missing}")
            
            # Column selector
            st.write("Select columns to display:")
            all_cols = list(linear_data.columns)
            
            # Pre-select important columns
            important_cols = ['size_nm', 'number_normalized_avg', 'number_normalized_cumsum_avg', 
                            'volume_nm^3_per_mL_avg', 'volume_nm^3_per_mL_cumsum_avg',
                            'area_nm^2_per_mL_avg', 'area_nm^2_per_mL_cumsum_avg']
            
            selected_cols = st.multiselect(
                "Columns",
                all_cols,
                default=[c for c in important_cols if c in all_cols],
                key="linear_cols"
            )
            
            if selected_cols:
                st.dataframe(linear_data[selected_cols].head(20), use_container_width=True)
                
                # Download button for this view
                csv = linear_data[selected_cols].to_csv(index=False, sep='\t')
                st.download_button(
                    label="Download Linear PSD (TSV)",
                    data=csv,
                    file_name=f"{results['metadata'].get('persistentID', 'data')}_PSD_LINEAR_preview.txt",
                    mime="text/plain"
                )
        
        with tab_log:
            st.write(f"Logarithmic scale: {len(log_data)} size bins")
            
            # Column information
            with st.expander("Column Information (Logarithmic Scale)"):
                st.write(f"Total columns: {len(log_data.columns)}")
                col_groups = {
                    'Size': ['size_nm'],
                    'Raw Counts': ['number_avg', 'number_sd'],
                    'Metadata': ['num_replicates', 'source_files'],
                    'Concentration (dilution-corrected)': ['particles_per_mL_avg', 'particles_per_mL_sd'],
                    'Volume (dilution-corrected)': ['volume_nm^3_per_mL_avg', 'volume_nm^3_per_mL_sd'],
                    'Surface Area (dilution-corrected)': ['area_nm^2_per_mL_avg', 'area_nm^2_per_mL_sd'],
                    'Normalized Number [CRITICAL]': ['number_normalized_avg', 'number_normalized_sd'],
                    'Cumulative Normalized [CRITICAL]': ['number_normalized_cumsum_avg', 'number_normalized_cumsum_sd'],
                    'Cumulative Volume': ['volume_nm^3_per_mL_cumsum_avg', 'volume_nm^3_per_mL_cumsum_sd'],
                    'Cumulative Surface Area': ['area_nm^2_per_mL_cumsum_avg', 'area_nm^2_per_mL_cumsum_sd']
                }
                
                for group_name, cols in col_groups.items():
                    available = [c for c in cols if c in log_data.columns]
                    missing = [c for c in cols if c not in log_data.columns]
                    
                    status = "OK" if not missing else "MISSING"
                    st.write(f"**{group_name}** [{status}]")
                    for col in available:
                        st.write(f"  - {col}")
                    if missing:
                        st.warning(f"  Missing: {missing}")
            
            # Column selector
            st.write("Select columns to display:")
            all_cols = list(log_data.columns)
            
            # Pre-select important columns
            important_cols = ['size_nm', 'number_normalized_avg', 'number_normalized_cumsum_avg', 
                            'volume_nm^3_per_mL_avg', 'volume_nm^3_per_mL_cumsum_avg',
                            'area_nm^2_per_mL_avg', 'area_nm^2_per_mL_cumsum_avg']
            
            selected_cols = st.multiselect(
                "Columns",
                all_cols,
                default=[c for c in important_cols if c in all_cols],
                key="log_cols"
            )
            
            if selected_cols:
                st.dataframe(log_data[selected_cols].head(20), use_container_width=True)
                
                # Download button for this view
                csv = log_data[selected_cols].to_csv(index=False, sep='\t')
                st.download_button(
                    label="Download Logarithmic PSD (TSV)",
                    data=csv,
                    file_name=f"{results['metadata'].get('persistentID', 'data')}_PSD_LOGARITHMIC_preview.txt",
                    mime="text/plain"
                )
    
    # TAB 2: Metadata
    with tab2:
        st.subheader("Metadata")
        
        metadata = results['metadata']
        
        # Define sections in order (same as in save_metadata_file)
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
            ]
        }
        
        # Build ordered list matching sections
        metadata_rows = []
        for section_name, field_names in sections.items():
            for field_name in field_names:
                if field_name in metadata:
                    metadata_rows.append({'Field': field_name, 'Value': metadata[field_name]})
        
        metadata_df = pd.DataFrame(metadata_rows)
        st.dataframe(metadata_df, use_container_width=True, hide_index=True)
    
    # TAB 3: Warnings & Discrepancies
    with tab3:
        # Check if there are any concerning issues
        has_alerts = bool(results['quality_alerts'])
        has_variation = bool(results['high_variation_fields'])
        
        if not has_alerts and not has_variation:
            # All good!
            st.success("No quality issues detected!")
            st.write("""
            Your measurement looks good:
            - No quality control alerts
            - No high variation in measurement parameters
            - All data consistent across replicates
            """)
        else:
            st.subheader("Concerning Items")
            
            # Quality control alerts
            if results['quality_alerts']:
                st.error("Quality Control Alerts")
                for alert in results['quality_alerts']:
                    st.write(f"- {alert}")
                st.write("**Recommendation:** Review measurement conditions and consider if data is suitable for publication.")
            
            # High variation fields
            if results['high_variation_fields']:
                st.warning("High Variation Between Replicates")
                for field in results['high_variation_fields']:
                    st.write(f"- {field}")
                st.write("**Recommendation:** Check sample consistency, mixing, and instrument stability.")
    
    # TAB 4: Metrics (Cell 05 & 06 Results)
    with tab4:
        st.subheader("Metrics")
        
        # Add description of span
        st.info(
            "**Span** measures particle size distribution width: "
            "Span = (D90 - D10) / D50. "
            "Lower span = narrower distribution (more uniform), "
            "Higher span = broader distribution (more polydisperse)."
        )
        
        if 'total_metrics' in results and results['total_metrics']:
            metrics_dict = results['total_metrics']
            
            # Build formatted metrics list - only LINEAR scale
            metrics_rows = []
            
            # Use LINEAR scale only (not logarithmic)
            if 'linear' in metrics_dict:
                scale_metrics = metrics_dict['linear']
                
                # Total particles per mL
                if 'total_particles_per_mL_avg' in scale_metrics:
                    avg = scale_metrics['total_particles_per_mL_avg']
                    sd = scale_metrics.get('total_particles_per_mL_sd', 0)
                    metrics_rows.append({
                        'Field': f'nta_total_particles_per_mL',
                        'Value': f'{avg:.2E} +/- {sd:.2E}'
                    })
                
                # Total volume per mL (in uL)
                if 'total_volume_uL_per_mL_avg' in scale_metrics:
                    avg = scale_metrics['total_volume_uL_per_mL_avg']
                    sd = scale_metrics.get('total_volume_uL_per_mL_sd', 0)
                    metrics_rows.append({
                        'Field': f'nta_total_volume_uL_per_mL',
                        'Value': f'{avg:.4E} +/- {sd:.4E}'
                    })
                
                # Volume percentage
                if 'volume_percentage_avg' in scale_metrics:
                    avg = scale_metrics['volume_percentage_avg']
                    sd = scale_metrics.get('volume_percentage_sd', 0)
                    metrics_rows.append({
                        'Field': f'nta_volume_percentage',
                        'Value': f'{avg:.6E} +/- {sd:.6E}'
                    })
                
                # Specific surface area
                if 'specific_surface_area_m^2_per_cm^3_avg' in scale_metrics:
                    avg = scale_metrics['specific_surface_area_m^2_per_cm^3_avg']
                    metrics_rows.append({
                        'Field': f'nta_specific_surface_area_m^2_per_cm^3',
                        'Value': f'{avg:.2f}'
                    })
            
            # Add metrics_scale
            metrics_rows.append({
                'Field': 'nta_metrics_scale',
                'Value': 'linear'
            })
            
            # Add all D-values from percentile_stats (number, volume, surface area)
            if 'percentile_stats' in results and 'linear' in results['percentile_stats']:
                linear_stats = results['percentile_stats']['linear']
                
                # Define distribution types in order
                dist_types = [
                    ('number', 'Number-weighted (normalized)'),
                    ('volume', 'Volume-weighted'),
                    ('surface_area', 'Surface area-weighted')
                ]
                
                missing_types = []
                for dist_key, dist_label in dist_types:
                    if dist_key in linear_stats:
                        stats = linear_stats[dist_key]
                        
                        # Add distribution header
                        metrics_rows.append({
                            'Field': f'--- {dist_label} ---',
                            'Value': ''
                        })
                        
                        # Add D10, D50, D90, Span
                        for param in ['D10', 'D50', 'D90']:
                            avg_val = stats.get(f'{param}_avg')
                            lower_val = stats.get(f'{param}_lower')
                            upper_val = stats.get(f'{param}_upper')
                            
                            if avg_val is not None and not np.isnan(avg_val):
                                metrics_rows.append({
                                    'Field': f'{dist_key}_{param.lower()}',
                                    'Value': f'{avg_val:.2f} nm ({lower_val:.2f} - {upper_val:.2f})'
                                })
                        
                        # Add span
                        span_avg = stats.get('span_avg')
                        span_lower = stats.get('span_lower')
                        span_upper = stats.get('span_upper')
                        
                        if span_avg is not None and not np.isnan(span_avg):
                            metrics_rows.append({
                                'Field': f'{dist_key}_span',
                                'Value': f'{span_avg:.3f} ({span_lower:.3f} - {span_upper:.3f})'
                            })
                    else:
                        missing_types.append(dist_label)
                
                # Show warning if any distribution types are missing
                if missing_types:
                    st.warning(f"Missing D-values for: {', '.join(missing_types)}")
                    st.info("This may occur if required cumulative distribution columns are missing. Check the downloaded PSD file.")
            
            if metrics_rows:
                metrics_df = pd.DataFrame(metrics_rows)
                st.dataframe(metrics_df, use_container_width=True, hide_index=True)
        else:
            st.info("No metrics calculated yet")
    
    # TAB 5: Download
    with tab5:
        st.subheader("Download")
        
        output_dir = tempfile.mkdtemp()
        
        try:
            created_files = st.session_state.analyzer.save_outputs(output_dir)
            
            # Individual download buttons
            for filepath in sorted(created_files):
                filename = os.path.basename(filepath)
                with open(filepath, 'r') as f:
                    file_content = f.read()
                
                st.download_button(
                    label=f"Download {filename}",
                    data=file_content,
                    file_name=filename,
                    mime="text/plain"
                )
            
            # Download all as zip
            st.markdown("---")
            
            import zipfile
            import io
            
            zip_buffer = io.BytesIO()
            with zipfile.ZipFile(zip_buffer, 'w', zipfile.ZIP_DEFLATED) as zip_file:
                for filepath in created_files:
                    filename = os.path.basename(filepath)
                    with open(filepath, 'r') as f:
                        zip_file.writestr(filename, f.read())
            
            zip_buffer.seek(0)
            sample_id = results['metadata'].get('persistentID', 'analysis')
            
            st.download_button(
                label="Download All (ZIP)",
                data=zip_buffer.getvalue(),
                file_name=f"{sample_id}_all.zip",
                mime="application/zip"
            )
            
        except Exception as e:
            st.error(f"Error: {str(e)}")
else:
    st.info("Upload NTA files and click 'Analyze Files'")

st.markdown("---")
st.markdown("**NTA Analysis** | Cells 01-06 Implementation")
