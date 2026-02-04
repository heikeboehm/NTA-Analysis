"""
NTA Data Analysis - Streamlit Web Application
Updated to use combined analyzer with Cells 01-04
Includes metadata discrepancy detection and warnings
"""

import streamlit as st
import pandas as pd
import io
import tempfile
import os
from nta_analyzer_cells_01_04 import NTAAnalyzer, CONFIG

# Set page config
st.set_page_config(
    page_title="NTA Analysis",
    page_icon="🧪",
    layout="wide",
    initial_sidebar_state="expanded"
)

st.title("🧪 NTA Particle Size Analysis")
st.markdown("---")

# Initialize session state
if 'analyzer' not in st.session_state:
    st.session_state.analyzer = None
if 'results' not in st.session_state:
    st.session_state.results = None

# Sidebar configuration
with st.sidebar:
    st.header("⚙️ Configuration")
    
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
        st.info(f"ℹ️ Will use: {manual_persistent_id}")
    
    # Update CONFIG with user values
    CONFIG["project_metadata"]["experimenter"] = experimenter
    CONFIG["project_metadata"]["location"] = location
    CONFIG["project_metadata"]["project"] = project
    CONFIG["project_metadata"]["pi"] = pi
    CONFIG["project_metadata"]["funding"] = funding if funding else "none"
    
    # Store manual persistent ID if provided
    if manual_persistent_id:
        CONFIG["manual_persistent_id"] = manual_persistent_id
    
    st.info("ℹ️ Configuration changes apply when you click 'Analyze Files'")
    
    # Reset button
    if st.button("🔄 Reset All"):
        st.session_state.clear()
        st.rerun()

# Main content
st.header("📤 Upload NTA Files")
st.markdown("Upload one or more NTA data files (.txt format)")

uploaded_files = st.file_uploader(
    "Choose NTA files",
    type="txt",
    accept_multiple_files=True,
    help="Select one or more NTA data files for analysis"
)

if uploaded_files:
    if st.button("🔍 Analyze Files", key="analyze_btn", type="primary"):
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
                
                st.success("✅ Analysis completed!")
                
            except Exception as e:
                st.error(f"❌ Error during analysis: {str(e)}")
                st.stop()

# Display results if analysis was successful
if st.session_state.results:
    results = st.session_state.results
    
    st.markdown("---")
    st.header("📊 Analysis Results")
    
    # Tabs for different views
    tab1, tab2, tab3, tab4, tab5 = st.tabs([
        "📈 Distribution Data",
        "🔍 Metadata",
        "⚠️ Warnings",
        "📊 Metrics",
        "💾 Download"
    ])
    
    # TAB 1: Distribution Data
    with tab1:
        st.subheader("Distribution Data")
        
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
                st.metric("Status", "✓ Good")
        
        # Show first 9 rows of linear data
        st.subheader("Data Preview (Linear)")
        linear_data = results['distribution'][results['distribution']['scale'] == 'linear']
        st.dataframe(linear_data.head(9), use_container_width=True)
    
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
                'meta_version', 'python_analysis'
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
            st.success("✅ No quality issues detected!")
            st.write("""
            Your measurement looks good:
            - No quality control alerts
            - No high variation in measurement parameters
            - All data consistent across replicates
            """)
        else:
            st.subheader("⚠️ Concerning Items")
            
            # Quality control alerts
            if results['quality_alerts']:
                st.error("🚨 **Quality Control Alerts**")
                for alert in results['quality_alerts']:
                    st.write(f"- {alert}")
                st.write("**Recommendation:** Review measurement conditions and consider if data is suitable for publication.")
            
            # High variation fields
            if results['high_variation_fields']:
                st.warning("📊 **High Variation Between Replicates**")
                for field in results['high_variation_fields']:
                    st.write(f"- {field}")
                st.write("**Recommendation:** Check sample consistency, mixing, and instrument stability.")
    
    # TAB 4: Metrics (Cell 05 Results)
    with tab4:
        st.subheader("Metrics")
        
        if 'total_metrics' in results and results['total_metrics']:
            metrics_dict = results['total_metrics']
            
            # Build formatted metrics list matching original metadata style
            metrics_rows = []
            
            for scale, scale_metrics in metrics_dict.items():
                # Total particles per mL
                if 'total_particles_per_mL_avg' in scale_metrics:
                    avg = scale_metrics['total_particles_per_mL_avg']
                    sd = scale_metrics.get('total_particles_per_mL_sd', 0)
                    metrics_rows.append({
                        'Field': f'nta_total_particles_per_mL',
                        'Value': f'{avg:.2E} ± {sd:.2E}'
                    })
                
                # Total volume per mL (in uL)
                if 'total_volume_uL_per_mL_avg' in scale_metrics:
                    avg = scale_metrics['total_volume_uL_per_mL_avg']
                    sd = scale_metrics.get('total_volume_uL_per_mL_sd', 0)
                    metrics_rows.append({
                        'Field': f'nta_total_volume_uL_per_mL',
                        'Value': f'{avg:.4E} ± {sd:.4E}'
                    })
                
                # Volume percentage
                if 'volume_percentage_avg' in scale_metrics:
                    avg = scale_metrics['volume_percentage_avg']
                    sd = scale_metrics.get('volume_percentage_sd', 0)
                    metrics_rows.append({
                        'Field': f'nta_volume_percentage',
                        'Value': f'{avg:.6E} ± {sd:.6E}'
                    })
                
                # Total surface area (in cm²)
                if 'total_surface_area_cm^2_per_mL_avg' in scale_metrics:
                    avg = scale_metrics['total_surface_area_cm^2_per_mL_avg']
                    sd = scale_metrics.get('total_surface_area_cm^2_per_mL_sd', 0)
                    metrics_rows.append({
                        'Field': f'nta_total_surface_area_cm^2_per_mL',
                        'Value': f'{avg:.4E} ± {sd:.4E}'
                    })
                
                # Specific surface area
                if 'specific_surface_area_m^2_per_cm^3_avg' in scale_metrics:
                    avg = scale_metrics['specific_surface_area_m^2_per_cm^3_avg']
                    metrics_rows.append({
                        'Field': f'nta_specific_surface_area_m^2_per_cm^3',
                        'Value': f'{avg:.2f}'
                    })
            
            # Add metadata about metrics
            metrics_rows.append({
                'Field': 'nta_metrics_scale',
                'Value': 'linear'
            })
            metrics_rows.append({
                'Field': 'nta_metrics_replicates',
                'Value': str(results['num_replicates'])
            })
            
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
                    label=f"📥 {filename}",
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
                label="📦 Download All (ZIP)",
                data=zip_buffer.getvalue(),
                file_name=f"{sample_id}_all.zip",
                mime="application/zip"
            )
            
        except Exception as e:
            st.error(f"Error: {str(e)}")
else:
    st.info("Upload NTA files and click 'Analyze Files'")

st.markdown("---")
st.markdown("**NTA Analysis** | Cells 01-04")
