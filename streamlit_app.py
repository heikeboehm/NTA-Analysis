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
    
    # Update CONFIG with user values
    CONFIG["project_metadata"]["experimenter"] = experimenter
    CONFIG["project_metadata"]["location"] = location
    CONFIG["project_metadata"]["project"] = project

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
    st.write(f"✓ {len(uploaded_files)} file(s) selected")
    
    # Show file list
    with st.expander("📋 File List", expanded=True):
        for i, file in enumerate(uploaded_files, 1):
            st.write(f"{i}. {file.name}")
    
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
                
                st.success("✅ Analysis completed successfully!")
                
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
        "⚠️ Discrepancies",
        "📋 Details",
        "💾 Download"
    ])
    
    # TAB 1: Distribution Data
    with tab1:
        st.subheader("Particle Size Distribution")
        
        # Key statistics
        col1, col2, col3, col4 = st.columns(4)
        with col1:
            st.metric("Files Processed", results['num_replicates'])
        with col2:
            st.metric("Size Bins", len(results['distribution']))
        with col3:
            st.metric("Size Range", f"{results['distribution']['size_nm'].min():.1f} - {results['distribution']['size_nm'].max():.0f} nm")
        with col4:
            st.metric("Scales", f"{results['distribution']['scale'].nunique()}")
        
        # Distribution preview
        st.subheader("Distribution Data Preview")
        st.dataframe(results['distribution'].head(10), use_container_width=True)
        
        # Show statistics by scale
        st.subheader("Statistics by Scale")
        for scale in ['linear', 'logarithmic']:
            scale_data = results['distribution'][results['distribution']['scale'] == scale]
            if not scale_data.empty:
                st.write(f"**{scale.capitalize()} Scale:** {len(scale_data)} bins")
    
    # TAB 2: Metadata
    with tab2:
        st.subheader("Standardized Metadata")
        
        # Show key metadata
        col1, col2 = st.columns(2)
        with col1:
            st.write("**Project Information**")
            project_fields = ['persistentID', 'sample', 'experimenter', 'location', 'project']
            for field in project_fields:
                if field in results['metadata']:
                    st.write(f"- **{field}:** {results['metadata'][field]}")
        
        with col2:
            st.write("**NTA Parameters**")
            nta_fields = [k for k in results['metadata'].keys() if k.startswith('nta_')]
            for field in sorted(nta_fields)[:5]:
                st.write(f"- **{field}:** {results['metadata'][field]}")
        
        # Full metadata table
        st.subheader("All Metadata Fields")
        metadata_df = pd.DataFrame(
            list(results['metadata'].items()),
            columns=['Field', 'Value']
        )
        st.dataframe(metadata_df, use_container_width=True)
    
    # TAB 3: Metadata Discrepancies & Warnings
    with tab3:
        st.subheader("⚠️ Metadata Discrepancies & Quality Alerts")
        
        # Quality control alerts
        if results['quality_alerts']:
            st.warning("🚨 **Quality Control Alerts Detected**")
            for alert in results['quality_alerts']:
                st.write(f"- {alert}")
        else:
            st.info("✅ No quality control alerts")
        
        # High variation fields
        if results['high_variation_fields']:
            st.warning("📊 **High Variation Detected**")
            for field in results['high_variation_fields']:
                st.write(f"- {field}")
        else:
            st.info("✅ No high variation detected")
        
        # Different fields summary
        if results['different_fields']:
            st.subheader("Fields with Differences Across Files")
            st.write(f"**{len(results['different_fields'])} field(s) differ across files**")
            
            # Show which fields differ
            col1, col2 = st.columns(2)
            with col1:
                st.write("**Different Fields:**")
                for field in sorted(results['different_fields'].keys())[:10]:
                    values = results['different_fields'][field]
                    st.write(f"- {field}: {len(set(values))} unique value(s)")
            
            with col2:
                st.write("**Identical Fields:**")
                st.write(f"Total: {len(results['identical_fields'])} field(s)")
                for field in sorted(results['identical_fields'].keys())[:10]:
                    st.write(f"- {field}")
        else:
            st.info("✅ All extracted fields are identical across files")
        
        # Detailed comparison
        st.subheader("Detailed Field Comparison")
        
        if results['different_fields']:
            selected_field = st.selectbox(
                "Select a field to see values from each file:",
                sorted(results['different_fields'].keys())
            )
            
            if selected_field:
                values = results['different_fields'][selected_field]
                filenames = list(results['all_file_metadata'].keys())
                
                comparison_data = []
                for filename, value in zip(filenames, values):
                    comparison_data.append({
                        'File': filename,
                        'Value': value
                    })
                
                comparison_df = pd.DataFrame(comparison_data)
                st.dataframe(comparison_df, use_container_width=True)
                
                # Show explanation
                unique_values = set(values)
                if len(unique_values) > 1:
                    st.info(
                        f"⚠️ This field has **{len(unique_values)} different value(s)** across files:\n"
                        f"- {', '.join(str(v) for v in unique_values)}"
                    )
    
    # TAB 4: Details
    with tab4:
        st.subheader("Detailed Analysis Information")
        
        # File processing summary
        st.write("**File Processing Summary**")
        st.write(f"- Files processed successfully: {results['num_replicates']}")
        if results['failed_files']:
            st.write(f"- Files with errors: {len(results['failed_files'])}")
            for filename, error in results['failed_files'].items():
                st.write(f"  - {filename}: {error}")
        
        # Field analysis
        st.write("**Field Analysis**")
        field_analysis = results['field_analysis']
        
        analysis_summary = {
            'Identical': len(results['identical_fields']),
            'Different': len(results['different_fields']),
            'Total': len(field_analysis)
        }
        
        analysis_df = pd.DataFrame(
            list(analysis_summary.items()),
            columns=['Status', 'Count']
        )
        st.dataframe(analysis_df, use_container_width=True)
        
        # All files' metadata comparison
        st.subheader("Per-File Metadata Extraction")
        for filename, metadata in results['all_file_metadata'].items():
            with st.expander(f"📄 {filename}"):
                meta_df = pd.DataFrame(
                    list(metadata.items()),
                    columns=['Field', 'Value']
                )
                st.dataframe(meta_df, use_container_width=True)
    
    # TAB 5: Download
    with tab5:
        st.subheader("💾 Download Results")
        
        # Create output directory
        output_dir = tempfile.mkdtemp()
        
        # Save outputs
        try:
            created_files = st.session_state.analyzer.save_outputs(output_dir)
            
            # Create download buttons
            st.write("**Available files for download:**")
            
            for filepath in created_files:
                filename = os.path.basename(filepath)
                with open(filepath, 'rb') as f:
                    file_content = f.read()
                
                st.download_button(
                    label=f"📥 Download {filename}",
                    data=file_content,
                    file_name=filename,
                    mime="text/plain"
                )
            
            # Distribution as CSV
            dist_csv = results['distribution'].to_csv(index=False)
            st.download_button(
                label="📥 Download Distribution Data (CSV)",
                data=dist_csv,
                file_name=f"{results['metadata'].get('persistentID', 'analysis')}_distribution.csv",
                mime="text/csv"
            )
            
            # Metadata as JSON
            metadata_json = pd.DataFrame(
                list(results['metadata'].items()),
                columns=['Field', 'Value']
            ).to_json(orient='records')
            st.download_button(
                label="📥 Download Metadata (JSON)",
                data=metadata_json,
                file_name=f"{results['metadata'].get('persistentID', 'analysis')}_metadata.json",
                mime="application/json"
            )
            
            st.success("✅ All files ready for download")
            
        except Exception as e:
            st.error(f"Error preparing downloads: {str(e)}")
else:
    st.info("👆 Upload NTA files and click 'Analyze Files' to begin")

st.markdown("---")
st.markdown("**NTA Analysis App** | Cells 01-04 Integrated | © 2025")
