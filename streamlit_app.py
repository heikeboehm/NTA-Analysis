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
        "⚠️ Warnings",
        "📋 Details",
        "💾 Download"
    ])
    
    # TAB 1: Distribution Data
    with tab1:
        st.subheader("Particle Size Distribution")
        
        # Calculate quality status
        quality_status = "✓ GOOD"
        if results['quality_alerts']:
            quality_status = "❌ POOR"
        elif results['high_variation_fields']:
            quality_status = "⚠️ FAIR"
        
        # Get detected traces from metadata
        detected_traces = results['metadata'].get('nta_number_of_traces_sum', 'N/A')
        
        # Key statistics
        col1, col2, col3 = st.columns(3)
        with col1:
            st.metric("Files Processed", results['num_replicates'])
        with col2:
            st.metric("Detected Traces", detected_traces)
        with col3:
            st.metric("Quality Status", quality_status)
        
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
        
        # Scale selection guide
        st.info("""
        **📊 How to choose your scale:**
        
        - **LINEAR Scale**: Use when plotting with equal-width bins on X-axis
          - Best for: Statistical analysis, focus on small particles
          - When: Creating your own plots in spreadsheet/Python
        
        - **LOGARITHMIC Scale**: Use when plotting with log-spaced bins on X-axis  
          - Best for: Publication plots, wide size ranges, log-normal distributions
          - When: Standard choice for most applications
        
        👉 **Not sure?** Start with **LOGARITHMIC** - it's the publication standard.
        """)
        
        st.markdown("---")
        
        # Create output directory
        output_dir = tempfile.mkdtemp()
        
        # Save outputs
        try:
            created_files = st.session_state.analyzer.save_outputs(output_dir)
            
            # Organize files by type
            metadata_files = [f for f in created_files if 'metadata' in f.lower()]
            distribution_files = [f for f in created_files if 'psd' in f.lower()]
            
            if distribution_files:
                st.subheader("📈 Distribution Data Files (by Scale)")
                st.write("Choose the scale appropriate for your analysis:")
                
                for filepath in sorted(distribution_files):
                    filename = os.path.basename(filepath)
                    with open(filepath, 'r') as f:
                        file_content = f.read()
                    
                    # Show which scale this is
                    if 'LINEAR' in filename:
                        scale_label = "📊 LINEAR Scale"
                        scale_desc = "Equal bin widths on X-axis"
                    else:
                        scale_label = "📈 LOGARITHMIC Scale ⭐ Recommended"
                        scale_desc = "Equal spacing in log space - standard for publication"
                    
                    col1, col2 = st.columns([3, 1])
                    with col1:
                        st.write(f"**{scale_label}**")
                        st.caption(scale_desc)
                    with col2:
                        st.download_button(
                            label="📥 Download",
                            data=file_content,
                            file_name=filename,
                            mime="text/plain",
                            key=f"download_{filename}"
                        )
                    
                    # Show file preview
                    with st.expander(f"👁️ Preview: {filename}"):
                        # Show first few lines with header info
                        preview_lines = file_content.split('\n')[:12]
                        st.code('\n'.join(preview_lines), language='text')
            
            if metadata_files:
                st.subheader("📋 Metadata File")
                
                for filepath in metadata_files:
                    filename = os.path.basename(filepath)
                    with open(filepath, 'r') as f:
                        file_content = f.read()
                    
                    col1, col2 = st.columns([3, 1])
                    with col1:
                        st.write(f"**{filename}**")
                        st.caption("Standardized metadata with extraction notes")
                    with col2:
                        st.download_button(
                            label="📥 Download",
                            data=file_content,
                            file_name=filename,
                            mime="text/plain",
                            key=f"download_meta_{filename}"
                        )
            
            st.markdown("---")
            st.subheader("📦 Alternative Formats (for reference)")
            
            col1, col2 = st.columns(2)
            
            with col1:
                # Combined distribution as CSV (for reference)
                dist_csv = results['distribution'].to_csv(index=False)
                st.download_button(
                    label="📥 All Data Combined (CSV)",
                    data=dist_csv,
                    file_name=f"{results['metadata'].get('persistentID', 'analysis')}_distribution_combined.csv",
                    mime="text/csv",
                    help="Contains both LINEAR and LOGARITHMIC data in one file"
                )
            
            with col2:
                # Metadata as JSON
                metadata_json = pd.DataFrame(
                    list(results['metadata'].items()),
                    columns=['Field', 'Value']
                ).to_json(orient='records', indent=2)
                st.download_button(
                    label="📥 Metadata (JSON)",
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
