"""
NTA Data Analysis - Streamlit Web Application
Uses analyzer with Cells 01-05
"""

import streamlit as st
import pandas as pd
import tempfile
import os
import sys

# Add current directory to path to find nta_analyzer
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from nta_analyzer_cells_01_05 import run_analysis_pipeline, CONFIG

# Set page config
st.set_page_config(
    page_title="NTA Analysis",
    page_icon="microscope",
    layout="wide",
    initial_sidebar_state="expanded"
)

st.title("NTA Particle Size Analysis")
st.markdown("---")

# Initialize session state
if 'results' not in st.session_state:
    st.session_state.results = None

# Sidebar configuration
with st.sidebar:
    st.header("Configuration")
    
    st.subheader("Project Information")
    experimenter = st.text_input(
        "Experimenter",
        value=CONFIG["project_metadata"]["experimenter"]
    )
    location = st.text_input(
        "Lab Location",
        value=CONFIG["project_metadata"]["location"]
    )
    project = st.text_input(
        "Project Name",
        value=CONFIG["project_metadata"]["project"]
    )
    pi = st.text_input(
        "Principal Investigator",
        value=CONFIG["project_metadata"]["pi"]
    )
    funding = st.text_input(
        "Funding Source",
        value=CONFIG["project_metadata"]["funding"]
    )
    
    st.info("Configuration changes apply when you click 'Analyze Files'")
    
    if st.button("Reset All"):
        st.session_state.results = None
        st.rerun()

# File upload section
st.subheader("Upload NTA Files")

uploaded_files = st.file_uploader(
    "Select NTA data files",
    type="txt",
    accept_multiple_files=True,
    help="Upload one or more NTA raw data files for analysis"
)

if st.button("Analyze Files", type="primary"):
    if not uploaded_files:
        st.error("Please upload at least one file")
    else:
        with st.spinner("Analyzing files..."):
            try:
                # Create temporary files
                temp_dir = tempfile.mkdtemp()
                temp_files = []
                
                for uploaded_file in uploaded_files:
                    temp_path = os.path.join(temp_dir, uploaded_file.name)
                    with open(temp_path, 'wb') as f:
                        f.write(uploaded_file.getbuffer())
                    temp_files.append(temp_path)
                
                # Update config with sidebar values
                config = CONFIG.copy()
                config["project_metadata"]["experimenter"] = experimenter
                config["project_metadata"]["location"] = location
                config["project_metadata"]["project"] = project
                config["project_metadata"]["pi"] = pi
                config["project_metadata"]["funding"] = funding
                
                # Run analysis
                results = run_analysis_pipeline(temp_files, config=config)
                st.session_state.results = results
                
                st.success("Analysis completed successfully!")
                
            except Exception as e:
                st.error(f"Error during analysis: {str(e)}")
                import traceback
                st.error(traceback.format_exc())

# Display results if analysis was successful
if st.session_state.results is not None:
    results = st.session_state.results
    
    st.markdown("---")
    st.subheader("Results")
    
    # Create tabs for results
    tab1, tab2, tab3, tab4 = st.tabs([
        "Distribution Data",
        "Metadata",
        "Warnings",
        "Download"
    ])
    
    # TAB 1: Distribution Data
    with tab1:
        st.subheader("Particle Size Distribution")
        
        # Summary metrics
        col1, col2, col3 = st.columns(3)
        with col1:
            st.metric("Files Analyzed", results['num_replicates'])
        with col2:
            total_rows = len(results['distribution'])
            st.metric("Data Points", total_rows)
        with col3:
            st.metric("Scales", "Linear + Log")
        
        # Preview data
        st.write("**Linear Scale Data (first 9 rows):**")
        linear_data = results['distribution'][results['distribution']['scale'] == 'linear'].head(9)
        st.dataframe(linear_data, use_container_width=True)
    
    # TAB 2: Metadata
    with tab2:
        st.subheader("Sample Metadata")
        
        metadata = results['metadata']
        metadata_df = pd.DataFrame(list(metadata.items()), columns=['Field', 'Value'])
        st.dataframe(metadata_df, use_container_width=True)
    
    # TAB 3: Warnings
    with tab3:
        st.subheader("Quality Checks")
        
        field_analysis = results['field_analysis']
        if field_analysis:
            has_alerts = False
            for field, analysis in field_analysis.items():
                if analysis.get('alert'):
                    st.warning(f"**{field}**: {analysis.get('alert_message', 'Alert detected')}")
                    has_alerts = True
            if not has_alerts:
                st.success("No quality issues detected!")
        else:
            st.success("No quality issues detected!")
    
    # TAB 4: Download
    with tab4:
        st.subheader("Download Results")
        
        metadata = results['metadata']
        dist_df = results['distribution'].copy()
        sample_id = metadata.get('persistentID', 'analysis')
        
        col1, col2 = st.columns(2)
        
        # Download Metadata
        with col1:
            metadata_content = ""
            for key, value in sorted(metadata.items()):
                metadata_content += f"{key}\t{value}\n"
            
            st.download_button(
                label="Download Metadata",
                data=metadata_content,
                file_name=f"Data_{sample_id}_metadata.txt",
                mime="text/plain"
            )
        
        # Download PSD Data
        with col2:
            psd_content = dist_df.to_csv(sep='\t', index=False)
            
            st.download_button(
                label="Download PSD Data",
                data=psd_content,
                file_name=f"Data_{sample_id}_PSD_COMPREHENSIVE.txt",
                mime="text/plain"
            )
else:
    st.info("Upload NTA files and click 'Analyze Files' to begin")

st.markdown("---")
st.markdown("**NTA Analysis** | Cells 01-05 | Dilution correction, normalization, cumulative distributions")
