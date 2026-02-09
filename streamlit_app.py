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

# Import core analysis functions
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
    normalize_distributions_with_uncertainty,
    calculate_cumulative_distributions_with_uncertainty,
    calculate_percentile_statistics_with_uncertainty,
    calculate_concentration_totals_with_uncertainty,
    format_statistics_with_bounds,
)

# Import plotting functions
try:
    import nta_plots
    PLOTTING_AVAILABLE = True
except ImportError:
    PLOTTING_AVAILABLE = False

# Import download manager
try:
    import nta_download_manager
    DOWNLOAD_MANAGER_AVAILABLE = True
except ImportError:
    DOWNLOAD_MANAGER_AVAILABLE = False

# ============================================================================
# PAGE CONFIG & STYLING
# ============================================================================

st.set_page_config(
    page_title="NTA Data Analysis",
    page_icon="🧪",
    layout="wide",
    initial_sidebar_state="expanded"
)

st.markdown("""
    <style>
    .metric-card {
        background-color: #f0f2f6;
        padding: 20px;
        border-radius: 10px;
        text-align: center;
        margin: 10px 0;
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

# ============================================================================
# UTILITY FUNCTIONS
# ============================================================================

def metadata_to_tsv(metadata_dict):
    """Convert metadata dictionary to tab-separated text format"""
    lines = []
    for key, value in metadata_dict.items():
        lines.append(f"{key}\t{value}\t")
    return "\n".join(lines)



def distribution_to_tsv(df):
    """Convert distribution dataframe to tab-separated text format"""
    return df.to_csv(sep='\t', index=False)


def statistics_to_tsv(statistics_dict):
    """Convert statistics dictionary to tab-separated text format"""
    import pandas as pd
    stats_list = []
    
    for scale in statistics_dict:
        if isinstance(statistics_dict[scale], dict):
            for dist_type in statistics_dict[scale]:
                stat_dict = statistics_dict[scale][dist_type]
                row = {
                    'Scale': scale,
                    'Distribution': dist_type,
                    **stat_dict
                }
                stats_list.append(row)
    
    if stats_list:
        stats_df = pd.DataFrame(stats_list)
        return stats_df.to_csv(sep='\t', index=False)
    return ""
# SESSION STATE INITIALIZATION
# ============================================================================

def initialize_session_state():
    """Initialize session state variables"""
    if 'experimenter' not in st.session_state:
        st.session_state.experimenter = CONFIG['project_metadata'].get('experimenter', '')
    if 'project' not in st.session_state:
        st.session_state.project = CONFIG['project_metadata'].get('project', '')
    if 'location' not in st.session_state:
        st.session_state.location = CONFIG['project_metadata'].get('location', '')
    if 'pi' not in st.session_state:
        st.session_state.pi = CONFIG['project_metadata'].get('pi', '')
    if 'analysis_complete' not in st.session_state:
        st.session_state.analysis_complete = False
    if 'distribution_data' not in st.session_state:
        st.session_state.distribution_data = None
    if 'metadata' not in st.session_state:
        st.session_state.metadata = None
    if 'statistics' not in st.session_state:
        st.session_state.statistics = None
    if 'uploaded_files' not in st.session_state:
        st.session_state.uploaded_files = None
    if 'plot_output_dir' not in st.session_state:
        st.session_state.plot_output_dir = None
    if 'generated_plots' not in st.session_state:
        st.session_state.generated_plots = []

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

st.session_state.experimenter = st.sidebar.text_input(
    "Experimenter Initials",
    value=st.session_state.experimenter,
    key="exp_input"
)

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

# Run analysis button
if st.sidebar.button("▶️  RUN ANALYSIS", key="run_button", type="primary", use_container_width=True):
    if not uploaded_files:
        st.sidebar.error("❌ Please upload at least one file")
    else:
        st.session_state.analysis_complete = False
        st.session_state.run_analysis = True

# ============================================================================
# MAIN CONTENT AREA
# ============================================================================

st.title("🧪 NTA Data Analysis Tool")
st.markdown("Nanoparticle Tracking Analysis Data Processing")

# Processing section
if 'run_analysis' in st.session_state and st.session_state.run_analysis:
    st.markdown("---")
    
    with st.spinner("³ Processing files..."):
        try:
            # Step 1: Read files
            st.write("📖 Reading files...")
            
            files_data = []
            for uploaded_file in uploaded_files:
                file_content = uploaded_file.read().decode('latin1')
                success, sections = identify_sections(file_content)
                if success:
                    files_data.append((uploaded_file.name, file_content, sections))
                else:
                    st.warning(f"⚠️  Could not identify sections in {uploaded_file.name}")
            
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
                        st.warning(f"⚠️  Could not extract data from {filename}")
                
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
                        st.write("📊 Extracting metadata...")
                        
                        success, all_files_metadata = extract_metadata_from_all_files(files_data)
                        if success:
                            identical_fields, different_fields, field_analysis = analyze_field_differences(all_files_metadata)
                            metadata = create_automated_metadata(all_files_metadata, identical_fields, different_fields, CONFIG)
                            
                            # Update with user metadata
                            metadata['experimenter'] = st.session_state.experimenter
                            metadata['project'] = st.session_state.project
                            metadata['location'] = st.session_state.location
                            metadata['pi'] = st.session_state.pi
                            
                            sample_id_from_meta = metadata.get('sample', metadata.get('persistentID', 'unknown'))
                            
                            dilution_str = metadata.get('nta_dilution', '1.0')
                            try:
                                dilution_factor = float(dilution_str.split('Â±')[0].strip()) if 'Â±' in str(dilution_str) else float(dilution_str)
                            except (ValueError, TypeError):
                                dilution_factor = 1.0
                            
                            st.write(f"📌 Using dilution factor: {dilution_factor}")
                            st.write(f"📌 Sample ID: {sample_id_from_meta}")
                        else:
                            st.warning("âš  Could not extract metadata, using defaults")
                            metadata = {'experimenter': st.session_state.experimenter}
                        
                        # Step 5: Apply dilution correction
                        st.write("📢 Applying dilution correction...")
                        
                        success, distribution_df = apply_dilution_correction_with_uncertainty(
                            distribution_df,
                            metadata=metadata,
                            manual_dilution=dilution_factor
                        )
                        
                        if not success:
                            st.warning(f"⚠️  {distribution_df}")
                        
                        # Step 6: Calculate distributions
                        st.write("📊 Calculating normalized distributions...")
                        
                        success, distribution_df = normalize_distributions_with_uncertainty(distribution_df, metadata)
                        
                        if success:
                            success, distribution_df = calculate_cumulative_distributions_with_uncertainty(distribution_df)
                        
                        if not success:
                            st.warning(f"⚠️  Error calculating distributions: {distribution_df}")
                        
                        # Step 7: Calculate statistics
                        st.write("📉 Calculating statistics...")
                        
                        success, statistics = calculate_percentile_statistics_with_uncertainty(distribution_df)
                        
                        if not success:
                            st.warning(f"⚠️  Error calculating statistics: {statistics}")
                        else:
                            # Step 7b: Calculate concentration totals
                            st.write("📊 Calculating concentration totals...")
                            
                            success, totals = calculate_concentration_totals_with_uncertainty(distribution_df)
                            
                            if success:
                                # Add totals to metadata
                                metadata['nta_total_particles_per_mL'] = f"{totals['nta_total_particles_per_mL_avg']:.2e} ± {totals['nta_total_particles_per_mL_sd']:.2e}"
                                metadata['nta_total_volume_uL_per_mL'] = f"{totals['nta_total_volume_uL_per_mL_avg']:.4e} ± {totals['nta_total_volume_uL_per_mL_sd']:.4e}"
                                metadata['nta_total_surface_area_cm^2_per_mL'] = f"{totals['nta_total_surface_area_cm2_per_mL_avg']:.4e} ± {totals['nta_total_surface_area_cm2_per_mL_sd']:.4e}"
                                metadata['nta_volume_percentage'] = f"{totals['nta_volume_percentage_avg']:.6f} ± {totals['nta_volume_percentage_sd']:.6f}"
                                metadata['nta_specific_surface_area_m^2_per_cm^3'] = f"{totals['nta_specific_surface_area_m2_per_cm3']:.2f}"
                                
                                st.write(f"Total particles: {totals['nta_total_particles_per_mL_avg']:.2e}")
                                st.write(f"Total volume: {totals['nta_total_volume_uL_per_mL_avg']:.4e} µL/mL")
                                st.write(f"Total surface area: {totals['nta_total_surface_area_cm2_per_mL_avg']:.4e} cm²/mL")
                            else:
                                st.warning(f"Error: {totals}")
                            
                            # Step 7c: Format statistics with bounds
                            st.write("📉 Formatting statistics...")
                            formatted_stats = format_statistics_with_bounds(statistics)
                            metadata.update(formatted_stats)
                            
                            # Add scale and replicate info
                            metadata['nta_metrics_scale'] = 'linear'
                            metadata['nta_metrics_replicates'] = metadata.get('num_replicates', 'N/A')
                            
                            # Store results
                            st.session_state.distribution_data = distribution_df
                            st.session_state.metadata = metadata
                            st.session_state.statistics = statistics
                            st.session_state.analysis_complete = True
                            st.session_state.run_analysis = False
                            
                            st.success("Analysis complete!")
                            st.rerun()

                            st.rerun()
        
        except Exception as e:
            st.error(f"❌ Error during analysis: {str(e)}")
            st.session_state.run_analysis = False

# ============================================================================
# RESULTS TABS
# ============================================================================

if st.session_state.analysis_complete:
    st.markdown("---")
    
    tab_summary, tab_distributions, tab_statistics, tab_plots, tab_download = st.tabs(
        ["📊 Summary", "📈 Distributions", "📉 Statistics", "📊 Plots", "📥 Download"]
    )
    
    # SUMMARY TAB
    with tab_summary:
        st.subheader("📊 Analysis Summary")
        
        col1, col2, col3, col4 = st.columns(4)
        
        with col1:
            num_files = st.session_state.metadata.get('num_replicates', 'N/A')
            st.metric("Number of Files", num_files)
        
        with col2:
            total_tracks = st.session_state.metadata.get('nta_number_of_traces_sum', 'N/A')
            st.metric("Total Tracks", total_tracks)
        
        with col3:
            dist_df = st.session_state.distribution_data
            if 'concentration_cm-3_per_mL_avg' in dist_df.columns:
                linear_df = dist_df[dist_df['scale'] == 'linear']
                concentration = linear_df['concentration_cm-3_per_mL_avg'].sum() if not linear_df.empty else 0
                st.metric("Concentration (particles/mL)", f"{concentration:.2e}")
            else:
                st.metric("Concentration (particles/mL)", "N/A")
        
        with col4:
            if st.session_state.statistics:
                stats = st.session_state.statistics
                if 'linear' in stats and 'number' in stats['linear']:
                    d50 = stats['linear']['number'].get('D50_avg', 'N/A')
                    if isinstance(d50, (int, float)):
                        st.metric("D50 (nm)", f"{d50:.2f}")
                    else:
                        st.metric("D50 (nm)", d50)
                else:
                    st.metric("D50 (nm)", "N/A")
            else:
                st.metric("D50 (nm)", "N/A")
        
        st.markdown("---")
        st.subheader("📋 Sample Information")
        
        col1, col2 = st.columns(2)
        
        with col1:
            st.write(f"**Experimenter:** {st.session_state.experimenter}")
            st.write(f"**Project:** {st.session_state.project}")
            st.write(f"**Location:** {st.session_state.location}")
        
        with col2:
            st.write(f"**PI:** {st.session_state.pi}")
            dilution_val = st.session_state.metadata.get('nta_dilution', 'N/A')
            st.write(f"**Dilution Factor:** {dilution_val}")
        
        st.markdown("---")
        st.subheader("Quick Visualization")
        
        dist_df = st.session_state.distribution_data
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
    
    # DISTRIBUTIONS TAB
    with tab_distributions:
        st.subheader("📊 Distribution Data Tables")
        
        dist_df = st.session_state.distribution_data
        scale_filter = st.radio("Scale", ["linear", "logarithmic"], horizontal=True)
        
        filtered_df = dist_df[dist_df['scale'] == scale_filter].copy()
        
        if not filtered_df.empty:
            display_cols = ['size_nm']
            display_cols += [col for col in filtered_df.columns if '_avg' in col or '_sd' in col]
            display_cols = [col for col in display_cols if col in filtered_df.columns]
            
            st.dataframe(
                filtered_df[display_cols].sort_values('size_nm'),
                use_container_width=True
            )
        else:
            st.warning(f"No data available for {scale_filter} scale")
    
    # STATISTICS TAB
    with tab_statistics:
        st.subheader("📉 Statistical Summary")
        
        if st.session_state.statistics and isinstance(st.session_state.statistics, dict):
            stats = st.session_state.statistics
            
            stats_list = []
            
            try:
                for scale in ['linear', 'logarithmic']:
                    if scale in stats and isinstance(stats[scale], dict):
                        for dist_type in ['number', 'volume', 'surface_area']:
                            if dist_type in stats[scale]:
                                stat_dict = stats[scale][dist_type]
                                
                                # Format D values with confidence intervals
                                def format_d_value(avg_val, lower_val, upper_val):
                                    if isinstance(avg_val, (int, float)):
                                        return f"{avg_val:.2f} ({lower_val:.2f} - {upper_val:.2f})"
                                    return "N/A"
                                
                                d10_str = format_d_value(
                                    stat_dict.get('D10_avg', 'N/A'),
                                    stat_dict.get('D10_lower', 'N/A'),
                                    stat_dict.get('D10_upper', 'N/A')
                                ) if all(isinstance(stat_dict.get(k), (int, float)) for k in ['D10_avg', 'D10_lower', 'D10_upper']) else 'N/A'
                                
                                d50_str = format_d_value(
                                    stat_dict.get('D50_avg', 'N/A'),
                                    stat_dict.get('D50_lower', 'N/A'),
                                    stat_dict.get('D50_upper', 'N/A')
                                ) if all(isinstance(stat_dict.get(k), (int, float)) for k in ['D50_avg', 'D50_lower', 'D50_upper']) else 'N/A'
                                
                                d90_str = format_d_value(
                                    stat_dict.get('D90_avg', 'N/A'),
                                    stat_dict.get('D90_lower', 'N/A'),
                                    stat_dict.get('D90_upper', 'N/A')
                                ) if all(isinstance(stat_dict.get(k), (int, float)) for k in ['D90_avg', 'D90_lower', 'D90_upper']) else 'N/A'
                                
                                span_str = format_d_value(
                                    stat_dict.get('span_avg', 'N/A'),
                                    stat_dict.get('span_lower', 'N/A'),
                                    stat_dict.get('span_upper', 'N/A')
                                ) if all(isinstance(stat_dict.get(k), (int, float)) for k in ['span_avg', 'span_lower', 'span_upper']) else 'N/A'
                                
                                stats_list.append({
                                    'Scale': scale.capitalize(),
                                    'Distribution': dist_type.replace('_', ' ').title(),
                                    'D10 (nm)': d10_str,
                                    'D50 (nm)': d50_str,
                                    'D90 (nm)': d90_str,
                                    'Span': span_str,
                                })
                
                if stats_list:
                    stats_df = pd.DataFrame(stats_list)
                    st.dataframe(stats_df, use_container_width=True)
                else:
                    st.info("No statistics data available")
            except Exception as e:
                st.warning(f"Error processing statistics: {str(e)}")
        else:
            st.warning("No statistics available")
    
    # PLOTS TAB
    with tab_plots:
        st.subheader("📊 Sophisticated Data Visualizations")
        st.markdown("Your original sophisticated plotting code with lognormal fits, error bars, D-values, and confidence intervals")
        
        if not PLOTTING_AVAILABLE:
            st.error("⚠️ Plotting module not available. Please ensure nta_plots.py is in the project.")
        else:
            st.markdown("Select plots to generate from your NTA data:")
            
            # Create columns for plot selection
            col1, col2, col3 = st.columns(3)
            
            with col1:
                plot_number = st.checkbox("📊 Number-Weighted Distribution", value=True)
                plot_volume = st.checkbox("📊 Volume-Weighted Distribution", value=False)
            
            with col2:
                plot_surface = st.checkbox("📊 Surface Area-Weighted Distribution", value=False)
                plot_raw = st.checkbox("📊 Raw Particle Counts", value=False)
            
            with col3:
                plot_vol_theo = st.checkbox("📊 Counts vs Volume", value=False)
                plot_surf_theo = st.checkbox("📊 Counts vs Surface Area", value=False)
            
            if st.button("🎨 Generate Selected Plots", use_container_width=True):
                dist_df = st.session_state.distribution_data
                metadata = st.session_state.metadata
                stats = st.session_state.statistics
                
                # Create persistent output directory
                import tempfile
                output_dir = tempfile.mkdtemp(prefix='nta_plots_')
                st.session_state.plot_output_dir = output_dir
                st.session_state.generated_plots = []
                
                try:
                    # Number-weighted distribution
                    if plot_number:
                        with st.expander("📊 Number-Weighted Distribution", expanded=True):
                            try:
                                success, result = nta_plots.generate_number_plots(
                                    dist_df, 
                                    stats_dict=stats,
                                    uniqueID=metadata.get('persistentID', 'NTA_sample'),
                                    metadata=metadata,
                                    output_dir=output_dir,
                                    config=None
                                )
                                if success:
                                    st.success("✅ Number distribution plot generated")
                                    st.session_state.generated_plots.append('number')
                                    # Show generated files
                                    from pathlib import Path
                                    plot_files = list(Path(output_dir).glob('*number*.png'))
                                    if plot_files:
                                        from PIL import Image
                                        for plot_file in sorted(plot_files):
                                            img = Image.open(plot_file)
                                            st.image(img, use_column_width=True)
                                else:
                                    st.warning(f"Could not generate number distribution: {result}")
                            except Exception as e:
                                st.warning(f"Error: {str(e)[:200]}")
                    
                    # Volume-weighted distribution
                    if plot_volume:
                        with st.expander("📊 Volume-Weighted Distribution"):
                            try:
                                success, result = nta_plots.generate_volume_plots(
                                    dist_df,
                                    stats_dict=stats,
                                    uniqueID=metadata.get('persistentID', 'NTA_sample'),
                                    metadata=metadata,
                                    output_dir=output_dir,
                                    config=None
                                )
                                if success:
                                    st.success("✅ Volume distribution plot generated")
                                    st.session_state.generated_plots.append('volume')
                                    from pathlib import Path
                                    plot_files = list(Path(output_dir).glob('*volume*.png'))
                                    if plot_files:
                                        from PIL import Image
                                        for plot_file in sorted(plot_files):
                                            img = Image.open(plot_file)
                                            st.image(img, use_column_width=True)
                                else:
                                    st.warning(f"Could not generate volume distribution: {result}")
                            except Exception as e:
                                st.warning(f"Error: {str(e)[:200]}")
                    
                    # Surface area-weighted distribution
                    if plot_surface:
                        with st.expander("📊 Surface Area-Weighted Distribution"):
                            try:
                                success, result = nta_plots.generate_surface_area_plots(
                                    dist_df,
                                    stats_dict=stats,
                                    uniqueID=metadata.get('persistentID', 'NTA_sample'),
                                    metadata=metadata,
                                    output_dir=output_dir,
                                    config=None
                                )
                                if success:
                                    st.success("✅ Surface area distribution plot generated")
                                    st.session_state.generated_plots.append('surface_area')
                                    from pathlib import Path
                                    plot_files = list(Path(output_dir).glob('*surface*.png'))
                                    if plot_files:
                                        from PIL import Image
                                        for plot_file in sorted(plot_files):
                                            img = Image.open(plot_file)
                                            st.image(img, use_column_width=True)
                                else:
                                    st.warning(f"Could not generate surface area distribution: {result}")
                            except Exception as e:
                                st.warning(f"Error: {str(e)[:200]}")
                    
                    # Raw particle counts
                    if plot_raw:
                        with st.expander("📊 Raw Particle Counts"):
                            try:
                                success, result = nta_plots.generate_raw_particle_plots(
                                    dist_df,
                                    uniqueID=metadata.get('persistentID', 'NTA_sample'),
                                    metadata=metadata,
                                    output_dir=output_dir,
                                    config=None
                                )
                                if success:
                                    st.success("✅ Raw particle plot generated")
                                    st.session_state.generated_plots.append('raw')
                                    from pathlib import Path
                                    plot_files = list(Path(output_dir).glob('*raw*.png'))
                                    if plot_files:
                                        from PIL import Image
                                        for plot_file in sorted(plot_files):
                                            img = Image.open(plot_file)
                                            st.image(img, use_column_width=True)
                                else:
                                    st.warning(f"Could not generate raw particle plot: {result}")
                            except Exception as e:
                                st.warning(f"Error: {str(e)[:200]}")
                    
                    # Raw counts vs volume
                    if plot_vol_theo:
                        with st.expander("📊 Raw Counts vs Theoretical Volume"):
                            try:
                                success, result = nta_plots.generate_volume_theoretical_plots(
                                    dist_df,
                                    stats_dict=stats,
                                    uniqueID=metadata.get('persistentID', 'NTA_sample'),
                                    metadata=metadata,
                                    output_dir=output_dir,
                                    config=None
                                )
                                if success:
                                    st.success("✅ Volume theoretical plot generated")
                                    st.session_state.generated_plots.append('volume_theoretical')
                                    from pathlib import Path
                                    plot_files = list(Path(output_dir).glob('*volume_theo*.png'))
                                    if plot_files:
                                        from PIL import Image
                                        for plot_file in sorted(plot_files):
                                            img = Image.open(plot_file)
                                            st.image(img, use_column_width=True)
                                else:
                                    st.warning(f"Could not generate volume theoretical plot: {result}")
                            except Exception as e:
                                st.warning(f"Error: {str(e)[:200]}")
                    
                    # Raw counts vs surface area
                    if plot_surf_theo:
                        with st.expander("📊 Raw Counts vs Theoretical Surface Area"):
                            try:
                                success, result = nta_plots.generate_surface_area_theoretical_plots(
                                    dist_df,
                                    stats_dict=stats,
                                    uniqueID=metadata.get('persistentID', 'NTA_sample'),
                                    metadata=metadata,
                                    output_dir=output_dir,
                                    config=None
                                )
                                if success:
                                    st.success("✅ Surface area theoretical plot generated")
                                    st.session_state.generated_plots.append('surface_area_theoretical')
                                    from pathlib import Path
                                    plot_files = list(Path(output_dir).glob('*surface_theo*.png'))
                                    if plot_files:
                                        from PIL import Image
                                        for plot_file in sorted(plot_files):
                                            img = Image.open(plot_file)
                                            st.image(img, use_column_width=True)
                                else:
                                    st.warning(f"Could not generate surface area theoretical plot: {result}")
                            except Exception as e:
                                st.warning(f"Error: {str(e)[:200]}")
                    
                    st.success("✅ All selected plots generated successfully! Check the Download tab for PDF files.")
                
                except Exception as e:
                    st.error(f"❌ Error generating plots: {str(e)}")
    
    # DOWNLOAD TAB
    with tab_download:
        st.subheader("📥 Download Analysis Results")
        
        st.markdown("Export your complete analysis results - data, plots, fit parameters, and everything as ZIP:")
        st.info("📝 All files are exported as tab-separated TXT format for easy import into Excel or other programs\n📊 Plots saved as PDF (300 DPI) and PNG\n📈 Fit data includes distribution data and lognormal fit parameters")

        # Get unique ID for file naming
        persistent_id = st.session_state.metadata.get('persistentID', 'NTA_analysis')
        # Extract shorter unique ID (remove _00XX_avg3 suffix)
        import re
        unique_id = re.sub(r'_\d{4}_avg\d+$', '', persistent_id)
        
        # Get number of replicates
        num_replicates = st.session_state.metadata.get('num_replicates', 1)

        # Check if plots have been generated
        has_plots = st.session_state.plot_output_dir is not None and \
                    len(st.session_state.generated_plots) > 0
        
        if has_plots and DOWNLOAD_MANAGER_AVAILABLE:
            # Show summary of available files
            from pathlib import Path
            output_path = Path(st.session_state.plot_output_dir)
            
            # Count files
            pdf_files = list(output_path.glob('*Plot_*.pdf'))
            png_files = list(output_path.glob('*Plot_*.png'))
            fit_files = list(output_path.glob('*Fits_*.json'))
            
            col_summary1, col_summary2, col_summary3, col_summary4 = st.columns(4)
            with col_summary1:
                st.metric("📊 Data Files", 3)  # distribution, metadata, stats
            with col_summary2:
                st.metric("📈 Plot PDFs", len(pdf_files))
            with col_summary3:
                st.metric("🖼 Plot PNGs", len(png_files))
            with col_summary4:
                st.metric("📋 Fit Data", len(fit_files))
            
            st.markdown("---")
        
        # Data Files Section
        st.subheader("📊 Analysis Data")
        col1, col2, col3 = st.columns(3)
        
        with col1:
            dist_df = st.session_state.distribution_data
            tsv_data = dist_df.to_csv(sep='\t', index=False)
            st.download_button(
                label="📥 Distribution Data",
                data=tsv_data,
                file_name=f"Data_{unique_id}_avg{num_replicates}_distribution.txt",
                mime="text/plain",
                key="dist_txt"
            )
        
        with col2:
            if st.session_state.metadata:
                metadata_df = pd.DataFrame(
                    [(k, v) for k, v in st.session_state.metadata.items()],
                    columns=['Field', 'Value']
                )
                tsv_data = metadata_df.to_csv(sep='\t', index=False)
                st.download_button(
                    label="📥 Metadata",
                    data=tsv_data,
                    file_name=f"Data_{unique_id}_avg{num_replicates}_metadata.txt",
                    mime="text/plain",
                    key="meta_txt"
                )
        
        with col3:
            if st.session_state.statistics and isinstance(st.session_state.statistics, dict):
                stats_list = []
                stats = st.session_state.statistics
                
                try:
                    for scale in stats:
                        if isinstance(stats[scale], dict):
                            for dist_type in stats[scale]:
                                stat_dict = stats[scale][dist_type]
                                row = {
                                    'Scale': scale,
                                    'Distribution': dist_type,
                                    **stat_dict
                                }
                                stats_list.append(row)
                    
                    if stats_list:
                        stats_df = pd.DataFrame(stats_list)
                        tsv_data = stats_df.to_csv(sep='\t', index=False)
                        st.download_button(
                            label="📥 Statistics",
                            data=tsv_data,
                            file_name=f"Data_{unique_id}_avg{num_replicates}_statistics.txt",
                            mime="text/plain",
                            key="stats_txt"
                        )
                except Exception as e:
                    st.warning(f"Could not format statistics: {str(e)}")
        
        # Plot Files Section
        if has_plots:
            st.markdown("---")
            st.subheader("📊 Generated Plots")
            
            from pathlib import Path
            output_path = Path(st.session_state.plot_output_dir)
            
            pdf_files = sorted(output_path.glob('*Plot_*.pdf'))
            
            if pdf_files:
                st.markdown(f"**{len(pdf_files)} Plot PDFs** (300 DPI, publication quality)")
                
                # Create grid of download buttons
                cols = st.columns(min(3, len(pdf_files)))
                for idx, pdf_file in enumerate(pdf_files):
                    with cols[idx % 3]:
                        with open(pdf_file, 'rb') as f:
                            pdf_data = f.read()
                        # Extract plot type from filename
                        plot_type = pdf_file.stem.split('_')[-1]
                        st.download_button(
                            label=f"📥 {plot_type.replace('_', ' ').title()}",
                            data=pdf_data,
                            file_name=pdf_file.name,
                            mime="application/pdf",
                            key=f"plot_pdf_{pdf_file.stem}"
                        )
            
            # Fit Data Section
            st.markdown("---")
            st.subheader("📈 Fit Data (Distribution + Parameters)")
            st.markdown("For each plot, get the distribution data and lognormal fit parameters to plot elsewhere")
            
            # Create fit data for each generated plot type
            dist_df = st.session_state.distribution_data
            stats = st.session_state.statistics
            
            plot_fit_configs = [
                ('number', 'linear', 'Number Distribution (Linear)'),
                ('volume', 'linear', 'Volume Distribution (Linear)'),
                ('surface_area', 'linear', 'Surface Area Distribution (Linear)'),
            ]
            
            cols_fit = st.columns(3)
            for idx, (dist_type, scale, label) in enumerate(plot_fit_configs):
                if dist_type in st.session_state.generated_plots:
                    with cols_fit[idx % 3]:
                        try:
                            fit_tsv = nta_download_manager.export_fit_data_tsv(
                                dist_df, stats, scale, dist_type
                            )
                            if fit_tsv:
                                st.download_button(
                                    label=f"📥 {label}",
                                    data=fit_tsv,
                                    file_name=f"FitData_{unique_id}_avg{num_replicates}_{dist_type}_{scale}.txt",
                                    mime="text/plain",
                                    key=f"fit_{dist_type}_{scale}"
                                )
                        except Exception as e:
                            st.warning(f"Could not create fit data: {str(e)}")
        
        # ZIP Download Section
        st.markdown("---")
        st.subheader("📦 Complete Package")
        
        if has_plots and DOWNLOAD_MANAGER_AVAILABLE:
            col_zip1, col_zip2 = st.columns([3, 1])
            
            with col_zip1:
                st.markdown(f"""
                **Download everything as a single ZIP file:**
                - ✅ All data files (distribution, metadata, statistics)
                - ✅ All generated plots (PDFs + PNGs)
                - ✅ All fit data (TSV with parameters)
                """)
            
            with col_zip2:
                try:
                    zip_bytes = nta_download_manager.create_download_zip(
                        st.session_state.plot_output_dir,
                        unique_id,
                        num_replicates
                    )
                    
                    if zip_bytes:
                        st.download_button(
                            label="📦 Download ZIP",
                            data=zip_bytes,
                            file_name=f"Data_{unique_id}_avg{num_replicates}_complete.zip",
                            mime="application/zip",
                            key="complete_zip",
                            use_container_width=True
                        )
                except Exception as e:
                    st.warning(f"Could not create ZIP: {str(e)}")
        else:
            st.markdown("""
            **Generate plots to enable ZIP download:**
            - Go to the **📊 Plots** tab
            - Select plots and click "🎨 Generate Selected Plots"
            - Then return here to download PDFs and fit data
            """)

else:
    # No analysis yet
    st.info("👆 Upload NTA files and click 'RUN ANALYSIS' to begin")
    st.markdown("""
    ## How to use:
    
    1. **Upload files** - Use the sidebar to upload one or more NTA .txt files
    2. **Enter metadata** - Customize experimenter, project, and sample information
    3. **Run analysis** - Click the 'RUN ANALYSIS' button to process your data
    4. **View results** - Check the Summary, Distributions, and Statistics tabs
    5. **Generate plots** - Go to the Plots tab to create visualizations
    6. **Download** - Export results as TSV data files, PDF plots, and fit parameters
    7. **Package** - Download everything as ZIP
    """)
