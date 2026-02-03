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
    CONFIG["project_metadata"]["funding"] = funding if funding else "none"
    
    # Store manual persistent ID if provided
    if manual_persistent_id:
        CONFIG["manual_persistent_id"] = manual_persistent_id

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
    tab1, tab2, tab3, tab4 = st.tabs([
        "📈 Distribution Data",
        "🔍 Metadata",
        "⚠️ Warnings",
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
        
        # Simple table of all metadata
        metadata_df = pd.DataFrame([
            {'Field': k, 'Value': v} 
            for k, v in sorted(metadata.items())
        ])
        
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
    
    # TAB 4: Download
    with tab4:
        st.subheader("Download")
        
        output_dir = tempfile.mkdtemp()
        
        try:
            created_files = st.session_state.analyzer.save_outputs(output_dir)
            
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
            
        except Exception as e:
            st.error(f"Error: {str(e)}")
else:
    st.info("Upload NTA files and click 'Analyze Files'")

st.markdown("---")
st.markdown("**NTA Analysis** | Cells 01-04")
