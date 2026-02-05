"""
NTA DATA ANALYSIS - STREAMLIT WEB APPLICATION

Interactive web application for analyzing Nanoparticle Tracking Analysis (NTA) data
using Streamlit. Provides file upload, processing, visualization, and download capabilities.
"""

import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import io
import os
from datetime import datetime
from pathlib import Path

# Import the NTA analysis module
from nta_analysis import (
    CONFIG,
    read_nta_file,
    identify_sections,
    extract_single_file_distribution,
    average_replicate_data,
    extract_metadata_from_all_files,
    analyze_field_differences,
    create_automated_metadata,
    apply_dilution_correction_with_uncertainty,
    calculate_normalized_distributions,
    calculate_cumulative_distributions,
    calculate_comprehensive_statistics,
    create_number_weighted_plot,
    create_volume_weighted_plot,
    create_surface_area_weighted_plot,
    create_raw_count_plot,
    create_count_vs_surface_area_plot,
    create_count_vs_volume_plot,
)

# ============================================================================
# PAGE CONFIG & STYLING
# ============================================================================

st.set_page_config(
    page_title="NTA Data Analysis",
    page_icon="🧪",
    layout="wide",
    initial_sidebar_state="expanded"
)

# Custom CSS for better styling
st.markdown("""
    <style>
    .metric-card {
        background-color: #f0f2f6;
        padding: 20px;
        border-radius: 10px;
        text-align: center;
        margin: 10px 0;
    }
    .metric-title {
        font-size: 14px;
        color: #666;
        margin-bottom: 10px;
    }
    .metric-value {
        font-size: 28px;
        font-weight: bold;
        color: #1f77b4;
    }
    .warning-box {
        background-color: #fff3cd;
        border-left: 4px solid #ffc107;
        padding: 10px;
        margin: 10px 0;
        border-radius: 4px;
    }
    .success-box {
        background-color: #d4edda;
        border-left: 4px solid #28a745;
        padding: 10px;
        margin: 10px 0;
        border-radius: 4px;
    }
    </style>
""", unsafe_allow_html=True)

# ============================================================================
# SESSION STATE INITIALIZATION
# ============================================================================

def initialize_session_state():
    """Initialize session state variables"""
    
    # Metadata fields - pre-filled from CONFIG
    if 'experimenter' not in st.session_state:
        st.session_state.experimenter = CONFIG['project_metadata'].get('experimenter', '')
    
    if 'project' not in st.session_state:
        st.session_state.project = CONFIG['project_metadata'].get('project', '')
    
    if 'location' not in st.session_state:
        st.session_state.location = CONFIG['project_metadata'].get('location', '')
    
    if 'pi' not in st.session_state:
        st.session_state.pi = CONFIG['project_metadata'].get('pi', '')
    
    # Analysis-specific fields
    if 'sample_id' not in st.session_state:
        st.session_state.sample_id = ''
    
    if 'dilution_factor' not in st.session_state:
        st.session_state.dilution_factor = 1.0
    
    # Analysis results storage
    if 'analysis_complete' not in st.session_state:
        st.session_state.analysis_complete = False
    
    if 'distribution_data' not in st.session_state:
        st.session_state.distribution_data = None
    
    if 'metadata' not in st.session_state:
        st.session_state.metadata = None
    
    if 'statistics' not in st.session_state:
        st.session_state.statistics = None
    
    if 'plots' not in st.session_state:
        st.session_state.plots = {}
    
    if 'uploaded_files' not in st.session_state:
        st.session_state.uploaded_files = None

initialize_session_state()

# ============================================================================
# SIDEBAR - FILE UPLOAD & METADATA
# ============================================================================

st.sidebar.title("🧪 NTA Analysis Control")
st.sidebar.markdown("---")

# File upload section
st.sidebar.subheader("📁 File Upload")
uploaded_files = st.sidebar.file_uploader(
    "Upload NTA data files (.txt)",
    type=["txt"],
    accept_multiple_files=True,
    help="Upload one or more ZetaView NTA .txt files for analysis"
)

if uploaded_files:
    st.session_state.uploaded_files = uploaded_files
    st.sidebar.success(f"✅ {len(uploaded_files)} file(s) uploaded")

st.sidebar.markdown("---")

# Metadata section
st.sidebar.subheader("👤 User Metadata")

col1, col2 = st.sidebar.columns([3, 1])
with col1:
    st.session_state.experimenter = st.text_input(
        "Experimenter Initials",
        value=st.session_state.experimenter,
        key="exp_input"
    )

with col2:
    if st.button("🔄 Reset", help="Reset all metadata to defaults", key="reset_button"):
        st.session_state.experimenter = CONFIG['project_metadata'].get('experimenter', '')
        st.session_state.project = CONFIG['project_metadata'].get('project', '')
        st.session_state.location = CONFIG['project_metadata'].get('location', '')
        st.session_state.pi = CONFIG['project_metadata'].get('pi', '')
        st.session_state.sample_id = ''
        st.session_state.dilution_factor = 1.0
        st.rerun()

st.session_state.project = st.sidebar.text_input(
    "Project Name",
    value=st.session_state.project,
    key="proj_input"
)

st.session_state.location = st.sidebar.text_input(
    "Lab Location",
    value=st.session_state.location,
    key="loc_input"
)

st.session_state.pi = st.sidebar.text_input(
    "Principal Investigator",
    value=st.session_state.pi,
    key="pi_input"
)

st.sidebar.markdown("---")

# Analysis-specific parameters
st.sidebar.subheader("🔬 Analysis Parameters")

st.session_state.sample_id = st.sidebar.text_input(
    "Sample ID",
    value=st.session_state.sample_id,
    key="sample_input"
)

st.session_state.dilution_factor = st.sidebar.number_input(
    "Dilution Factor",
    value=st.session_state.dilution_factor,
    min_value=0.1,
    step=0.5,
    key="dilution_input"
)

st.sidebar.markdown("---")

# Run analysis button
if st.sidebar.button("▶️ RUN ANALYSIS", key="run_button", type="primary", use_container_width=True):
    if not uploaded_files:
        st.sidebar.error("❌ Please upload at least one file")
    else:
        st.session_state.analysis_complete = False
        st.session_state.run_analysis = True

# ============================================================================
# MAIN CONTENT AREA
# ============================================================================

st.title("🧪 NTA Data Analysis Tool")
st.markdown("Comprehensive analysis of Nanoparticle Tracking Analysis data")

# Processing section
if 'run_analysis' in st.session_state and st.session_state.run_analysis:
    st.markdown("---")
    
    with st.spinner("⏳ Processing files..."):
        try:
            # Step 1: Read files
            st.write("📖 Reading files...")
            
            files_data = []
            for uploaded_file in uploaded_files:
                # Convert uploaded file to temporary file-like object
                file_content = uploaded_file.read().decode('latin1')
                
                # Identify sections
                success, sections = identify_sections(file_content)
                if success:
                    files_data.append((uploaded_file.name, file_content, sections))
                else:
                    st.warning(f"⚠️ Could not identify sections in {uploaded_file.name}")
            
            if not files_data:
                st.error("❌ Could not process any files")
            else:
                # Step 2: Extract distribution data
                st.write("📊 Extracting distribution data...")
                
                distribution_dfs = []
                for filename, content, sections in files_data:
                    success, result = extract_single_file_distribution(content, sections, filename)
                    if success:
                        distribution_dfs.append(result)
                    else:
                        st.warning(f"⚠️ Could not extract data from {filename}")
                
                if not distribution_dfs:
                    st.error("❌ Could not extract data from any files")
                else:
                    # Step 3: Average replicates
                    st.write("📈 Averaging replicates...")
                    
                    filenames = [f[0] for f in files_data]
                    success, distribution_df = average_replicate_data(distribution_dfs, filenames)
                    
                    if not success:
                        st.error(f"❌ Error averaging data: {distribution_df}")
                    else:
                        # Step 4: Extract metadata
                        st.write("📋 Extracting metadata...")
                        
                        success, all_files_metadata = extract_metadata_from_all_files(files_data)
                        if success:
                            identical_fields, different_fields, field_analysis = analyze_field_differences(all_files_metadata)
                            metadata = create_automated_metadata(all_files_metadata, identical_fields, different_fields, CONFIG)
                            
                            # Update with user metadata
                            metadata['experimenter'] = st.session_state.experimenter
                            metadata['project'] = st.session_state.project
                            metadata['location'] = st.session_state.location
                            metadata['pi'] = st.session_state.pi
                            metadata['sample_id'] = st.session_state.sample_id
                        else:
                            st.warning("⚠️ Could not extract metadata, using defaults")
                            metadata = {'experimenter': st.session_state.experimenter}
                        
                        # Step 5: Apply dilution correction
                        st.write("🔢 Applying dilution correction...")
                        
                        success, distribution_df = apply_dilution_correction_with_uncertainty(
                            distribution_df,
                            metadata=metadata,
                            manual_dilution=st.session_state.dilution_factor
                        )
                        
                        if not success:
                            st.warning(f"⚠️ {distribution_df}")
                        
                        # Step 6: Calculate distributions
                        st.write("📊 Calculating normalized distributions...")
                        
                        success, distribution_df = calculate_normalized_distributions(distribution_df, metadata)
                        
                        if success:
                            success, distribution_df = calculate_cumulative_distributions(distribution_df)
                        
                        if not success:
                            st.warning(f"⚠️ Error calculating distributions: {distribution_df}")
                        
                        # Step 7: Calculate statistics
                        st.write("📉 Calculating statistics...")
                        
                        success, statistics = calculate_comprehensive_statistics(distribution_df)
                        
                        if not success:
                            st.warning(f"⚠️ Error calculating statistics: {statistics}")
                        else:
                            # Store results in session state
                            st.session_state.distribution_data = distribution_df
                            st.session_state.metadata = metadata
                            st.session_state.statistics = statistics
                            st.session_state.analysis_complete = True
                            st.session_state.run_analysis = False
                            
                            st.success("✅ Analysis complete!")
                            st.rerun()
        
        except Exception as e:
            st.error(f"❌ Error during analysis: {str(e)}")
            st.session_state.run_analysis = False

# ============================================================================
# RESULTS TABS
# ============================================================================

if st.session_state.analysis_complete:
    st.markdown("---")
    
    tab_summary, tab_distributions, tab_statistics, tab_plots, tab_metadata, tab_download = st.tabs(
        ["📊 Summary", "📈 Distributions", "📉 Statistics", "🎨 Plots", "📋 Metadata", "⬇️ Download"]
    )
    
    # ========================================================================
    # SUMMARY TAB
    # ========================================================================
    with tab_summary:
        st.subheader("Analysis Summary")
        
        # Key metrics
        col1, col2, col3, col4 = st.columns(4)
        
        with col1:
            num_files = st.session_state.metadata.get('num_replicates', 'N/A')
            st.metric("Number of Files", num_files)
        
        with col2:
            total_tracks = st.session_state.metadata.get('nta_number_of_traces_sum', 'N/A')
            st.metric("Total Tracks", total_tracks)
        
        with col3:
            # Find concentration value from distribution data
            dist_df = st.session_state.distribution_data
            if 'concentration_cm-3_per_mL_avg' in dist_df.columns:
                concentration = dist_df['concentration_cm-3_per_mL_avg'].sum()
                st.metric("Concentration", f"{concentration:.2e}")
            else:
                st.metric("Concentration", "N/A")
        
        with col4:
            # Find D50 from statistics
            if st.session_state.statistics:
                stats = st.session_state.statistics
                if 'linear' in stats and 'number' in stats['linear']:
                    d50 = stats['linear']['number'].get('D50_avg', 'N/A')
                    span = stats['linear']['number'].get('span_avg', 'N/A')
                    if isinstance(d50, (int, float)):
                        st.metric("D50 (nm)", f"{d50:.2f}")
                    else:
                        st.metric("D50 (nm)", d50)
                else:
                    st.metric("D50 (nm)", "N/A")
            else:
                st.metric("D50 (nm)", "N/A")
        
        st.markdown("---")
        
        # Quality Control Summary
        st.subheader("🔍 Quality Control")
        
        qc_alerts = st.session_state.metadata.get('quality_control_alerts', '[]')
        if qc_alerts and qc_alerts != '[]':
            st.markdown('<div class="warning-box"><strong>⚠️ QC Alerts:</strong><br>' + qc_alerts.replace('[', '').replace(']', '').replace("'", '') + '</div>', unsafe_allow_html=True)
        else:
            st.markdown('<div class="success-box"><strong>✅ No QC alerts</strong> - Data looks good!</div>', unsafe_allow_html=True)
        
        st.markdown("---")
        
        # Sample Information
        st.subheader("📋 Sample Information")
        
        col1, col2 = st.columns(2)
        
        with col1:
            st.write(f"**Experimenter:** {st.session_state.experimenter}")
            st.write(f"**Project:** {st.session_state.project}")
            st.write(f"**Location:** {st.session_state.location}")
        
        with col2:
            st.write(f"**PI:** {st.session_state.pi}")
            st.write(f"**Sample ID:** {st.session_state.sample_id}")
            st.write(f"**Dilution Factor:** {st.session_state.dilution_factor}")
        
        st.markdown("---")
        
        # Simple visualization
        st.subheader("📊 Quick Visualization")
        
        dist_df = st.session_state.distribution_data
        
        # Find linear scale number-weighted data
        linear_data = dist_df[(dist_df['scale'] == 'linear')]
        
        if not linear_data.empty and 'size_nm' in linear_data.columns and 'number_normalized_avg' in linear_data.columns:
            fig, ax = plt.subplots(figsize=(10, 5))
            
            ax.bar(linear_data['size_nm'], linear_data['number_normalized_avg'], 
                   color='#1f77b4', alpha=0.7, edgecolor='black', linewidth=0.5)
            
            ax.set_xlabel('Size (nm)', fontsize=11)
            ax.set_ylabel('Normalized Count', fontsize=11)
            ax.set_title('Number-Weighted Distribution (Linear Scale)', fontsize=12, fontweight='bold')
            ax.grid(axis='y', alpha=0.3)
            
            plt.tight_layout()
            st.pyplot(fig)
            plt.close(fig)
    
    # ========================================================================
    # DISTRIBUTIONS TAB
    # ========================================================================
    with tab_distributions:
        st.subheader("Distribution Data Tables")
        
        dist_df = st.session_state.distribution_data
        
        # Filter options
        col1, col2 = st.columns(2)
        
        with col1:
            scale_filter = st.radio("Scale", ["linear", "logarithmic"], horizontal=True)
        
        with col2:
            # Get available columns
            value_cols = [col for col in dist_df.columns if '_avg' in col]
            
            if value_cols:
                st.info(f"✓ Found {len(dist_df)} data points in {scale_filter} scale")
        
        # Filter data
        filtered_df = dist_df[dist_df['scale'] == scale_filter].copy()
        
        if not filtered_df.empty:
            # Select columns to display
            display_cols = ['size_nm']
            display_cols += [col for col in filtered_df.columns if '_avg' in col or '_sd' in col]
            display_cols = [col for col in display_cols if col in filtered_df.columns]
            
            st.dataframe(
                filtered_df[display_cols].sort_values('size_nm'),
                use_container_width=True,
                height=400
            )
            
            # Download as CSV
            csv_data = filtered_df[display_cols].to_csv(index=False)
            st.download_button(
                label="📥 Download as CSV",
                data=csv_data,
                file_name=f"distribution_{scale_filter}.csv",
                mime="text/csv"
            )
        else:
            st.warning(f"No data available for {scale_filter} scale")
    
    # ========================================================================
    # STATISTICS TAB
    # ========================================================================
    with tab_statistics:
        st.subheader("Statistical Summary")
        
        if st.session_state.statistics:
            stats = st.session_state.statistics
            
            # Create a nice summary table
            stats_list = []
            
            for scale in ['linear', 'logarithmic']:
                if scale in stats:
                    for dist_type in ['number', 'volume', 'surface_area']:
                        if dist_type in stats[scale]:
                            stat_dict = stats[scale][dist_type]
                            
                            stats_list.append({
                                'Scale': scale.capitalize(),
                                'Distribution': dist_type.replace('_', ' ').title(),
                                'D10': f"{stat_dict.get('D10_avg', 'N/A'):.2f}" if isinstance(stat_dict.get('D10_avg'), (int, float)) else 'N/A',
                                'D50': f"{stat_dict.get('D50_avg', 'N/A'):.2f}" if isinstance(stat_dict.get('D50_avg'), (int, float)) else 'N/A',
                                'D90': f"{stat_dict.get('D90_avg', 'N/A'):.2f}" if isinstance(stat_dict.get('D90_avg'), (int, float)) else 'N/A',
                                'Span': f"{stat_dict.get('span_avg', 'N/A'):.3f}" if isinstance(stat_dict.get('span_avg'), (int, float)) else 'N/A',
                            })
            
            if stats_list:
                stats_df = pd.DataFrame(stats_list)
                st.dataframe(stats_df, use_container_width=True)
                
                # Download detailed statistics
                stats_csv = stats_df.to_csv(index=False)
                st.download_button(
                    label="📥 Download Statistics",
                    data=stats_csv,
                    file_name="statistics.csv",
                    mime="text/csv"
                )
        else:
            st.warning("No statistics available")
    
    # ========================================================================
    # PLOTS TAB
    # ========================================================================
    with tab_plots:
        st.subheader("Distribution Visualizations")
        
        st.info("Generating plots... This may take a moment for the first time.")
        
        # Create plots on demand
        col1, col2 = st.columns(2)
        
        with col1:
            if st.button("Generate Number-Weighted Plot"):
                dist_df = st.session_state.distribution_data
                for is_log in [False, True]:
                    try:
                        fig, _ = create_number_weighted_plot(
                            dist_df,
                            is_log_scale=is_log,
                            stats_dict=st.session_state.statistics,
                            uniqueID=st.session_state.sample_id or "sample",
                            metadata=st.session_state.metadata
                        )
                        
                        if fig:
                            scale_name = "Log" if is_log else "Linear"
                            st.pyplot(fig)
                            
                            # Download button
                            buf = io.BytesIO()
                            fig.savefig(buf, format='png', dpi=300, bbox_inches='tight')
                            buf.seek(0)
                            
                            st.download_button(
                                label=f"📥 Download {scale_name} (PNG)",
                                data=buf,
                                file_name=f"number_weighted_{scale_name.lower()}.png",
                                mime="image/png"
                            )
                    except Exception as e:
                        st.warning(f"Could not generate {scale_name} plot: {str(e)}")
        
        with col2:
            if st.button("Generate Volume-Weighted Plot"):
                dist_df = st.session_state.distribution_data
                for is_log in [False, True]:
                    try:
                        fig, _ = create_volume_weighted_plot(
                            dist_df,
                            is_log_scale=is_log,
                            stats_dict=st.session_state.statistics,
                            uniqueID=st.session_state.sample_id or "sample",
                            metadata=st.session_state.metadata
                        )
                        
                        if fig:
                            scale_name = "Log" if is_log else "Linear"
                            st.pyplot(fig)
                            
                            buf = io.BytesIO()
                            fig.savefig(buf, format='png', dpi=300, bbox_inches='tight')
                            buf.seek(0)
                            
                            st.download_button(
                                label=f"📥 Download {scale_name} (PNG)",
                                data=buf,
                                file_name=f"volume_weighted_{scale_name.lower()}.png",
                                mime="image/png"
                            )
                    except Exception as e:
                        st.warning(f"Could not generate {scale_name} plot: {str(e)}")
    
    # ========================================================================
    # METADATA TAB
    # ========================================================================
    with tab_metadata:
        st.subheader("Complete Metadata")
        
        if st.session_state.metadata:
            metadata_df = pd.DataFrame(
                [(k, v) for k, v in st.session_state.metadata.items()],
                columns=['Field', 'Value']
            )
            
            st.dataframe(metadata_df, use_container_width=True, height=500)
            
            # Download as CSV
            csv_data = metadata_df.to_csv(index=False)
            st.download_button(
                label="📥 Download Metadata",
                data=csv_data,
                file_name="metadata.csv",
                mime="text/csv"
            )
        else:
            st.warning("No metadata available")
    
    # ========================================================================
    # DOWNLOAD TAB
    # ========================================================================
    with tab_download:
        st.subheader("Download Analysis Results")
        
        st.markdown("""
        Export your complete analysis results:
        - **CSV files** with all data and metadata
        - **PNG images** of all visualizations
        - **Summary report** with key findings
        """)
        
        col1, col2, col3 = st.columns(3)
        
        with col1:
            st.subheader("📊 Data")
            
            # Distribution data
            dist_df = st.session_state.distribution_data
            csv_data = dist_df.to_csv(index=False)
            st.download_button(
                label="📥 Distribution Data (CSV)",
                data=csv_data,
                file_name=f"{st.session_state.sample_id or 'sample'}_distribution.csv",
                mime="text/csv",
                key="dist_csv"
            )
            
            # Metadata
            if st.session_state.metadata:
                metadata_df = pd.DataFrame(
                    [(k, v) for k, v in st.session_state.metadata.items()],
                    columns=['Field', 'Value']
                )
                csv_data = metadata_df.to_csv(index=False)
                st.download_button(
                    label="📥 Metadata (CSV)",
                    data=csv_data,
                    file_name=f"{st.session_state.sample_id or 'sample'}_metadata.csv",
                    mime="text/csv",
                    key="meta_csv"
                )
        
        with col2:
            st.subheader("📉 Statistics")
            
            if st.session_state.statistics:
                # Convert statistics to DataFrame for download
                stats_list = []
                stats = st.session_state.statistics
                
                for scale in stats:
                    for dist_type in stats[scale]:
                        stat_dict = stats[scale][dist_type]
                        row = {
                            'Scale': scale,
                            'Distribution': dist_type,
                            **stat_dict
                        }
                        stats_list.append(row)
                
                stats_df = pd.DataFrame(stats_list)
                csv_data = stats_df.to_csv(index=False)
                st.download_button(
                    label="📥 Statistics (CSV)",
                    data=csv_data,
                    file_name=f"{st.session_state.sample_id or 'sample'}_statistics.csv",
                    mime="text/csv",
                    key="stats_csv"
                )
        
        with col3:
            st.subheader("🎨 Plots")
            
            st.info("Click 'Generate' in the Plots tab to create visualizations, then download them from there.")
        
        st.markdown("---")
        st.info("💡 All files include a timestamp and your metadata for complete traceability.")

else:
    # No analysis yet
    st.info("👆 Upload NTA files and click 'RUN ANALYSIS' to begin")
    st.markdown("""
    ## How to use:
    
    1. **Upload files** - Use the sidebar to upload one or more NTA .txt files
    2. **Enter metadata** - Customize experimenter, project, and sample information
    3. **Run analysis** - Click the 'RUN ANALYSIS' button to process your data
    4. **View results** - Check the Summary, Distributions, Statistics, and Plots tabs
    5. **Download** - Export your results in CSV, PNG, and PDF formats
    
    The app will automatically:
    - Read and validate your NTA data files
    - Average multiple replicates with uncertainty quantification
    - Extract comprehensive metadata
    - Calculate distribution statistics (D10, D50, D90, span)
    - Generate publication-quality visualizations
    """)
