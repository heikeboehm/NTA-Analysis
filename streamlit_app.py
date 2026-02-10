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
    if 'research_focus' not in st.session_state:
        st.session_state.research_focus = 'Number of particles (size distribution)'

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
st.sidebar.subheader("👤 User Metadata (Optional)")

# Toggle to include user metadata in analysis
include_user_metadata = st.sidebar.checkbox(
    "Add user metadata to analysis?",
    value=False,
    help="Include experimenter, project, location, and PI in the analysis metadata (optional)"
)

if include_user_metadata:
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
else:
    # If not included, clear them
    st.session_state.experimenter = ""
    st.session_state.project = ""
    st.session_state.location = ""
    st.session_state.pi = ""

st.sidebar.markdown("---")
st.sidebar.subheader("🔬 Research Focus")

research_focus = st.sidebar.radio(
    "What are you primarily interested in?",
    options=[
        "Number of particles (size distribution)",
        "Surface area distribution",
        "Internal volume distribution"
    ],
    help="This helps highlight the most relevant metrics",
    index=0
)

st.session_state.research_focus = research_focus

st.sidebar.markdown("---")

# Run analysis button
if st.sidebar.button("RUN ANALYSIS", key="run_button", type="primary", use_container_width=True):
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
                            
                            # AUTO-GENERATE ALL PLOTS in background
                            st.write("🎨 Generating all plots...")
                            try:
                                import tempfile
                                output_dir = tempfile.mkdtemp(prefix='nta_plots_')
                                st.session_state.plot_output_dir = output_dir
                                st.session_state.generated_plots = []
                                
                                # Define all 3 plotting functions
                                plot_functions = [
                                    ('number', nta_plots.generate_number_plots, "Number-Weighted"),
                                    ('volume', nta_plots.generate_volume_plots, "Volume-Weighted"),
                                    ('surface_area', nta_plots.generate_surface_area_plots, "Surface Area-Weighted"),
                                ]
                                
                                for plot_key, plot_func, plot_label in plot_functions:
                                    try:
                                        success, result = plot_func(
                                            distribution_df,
                                            stats_dict=statistics,
                                            uniqueID=metadata.get('persistentID', 'NTA_sample'),
                                            metadata=metadata,
                                            output_dir=output_dir,
                                            config=None
                                        )
                                        if success:
                                            st.session_state.generated_plots.append(plot_key)
                                            st.write(f"✅ {plot_label}")
                                        else:
                                            st.write(f"⚠️ {plot_label}: {result}")
                                    except Exception as e:
                                        st.write(f"⚠️ {plot_label}: {str(e)[:100]}")
                                
                                st.write(f"✅ Generated {len(st.session_state.generated_plots)} plot(s)")
                            except Exception as e:
                                st.warning(f"Plot generation had issues: {str(e)[:200]}")
                            
                            st.success("✅ Analysis complete! View plots in 📊 Plots tab, download in 📥 Download tab")
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
    
    tab_summary, tab_statistics, tab_plots, tab_download = st.tabs(
        ["📊 Summary", "📉 Statistics", "📊 Plots", "📥 Download"]
    )
    
    # SUMMARY TAB
    with tab_summary:
        st.subheader("📊 Analysis Overview")
        
        dist_df = st.session_state.distribution_data
        stats = st.session_state.statistics if st.session_state.statistics else {}
        research_focus = st.session_state.get('research_focus', 'Number of particles (size distribution)')
        
        # Determine which distribution to show based on research focus
        if 'Number of particles' in research_focus:
            dist_type = 'number'
            total_metric_label = "Total Concentration"
            total_metric_col = 'concentration_cm-3_per_mL_avg'
            total_metric_unit = "particles/mL"
            metric_title = "Key metrics for particle number distribution"
            focus_emoji = "📊"
            focus_color = "#4C5B5C"
        elif 'Surface area' in research_focus:
            dist_type = 'surface_area'
            total_metric_label = "Total Surface Area"
            total_metric_col = 'area_nm^2_per_mL_avg'
            total_metric_unit = "nm²/mL"
            metric_title = "Key metrics for surface area distribution"
            focus_emoji = "📈"
            focus_color = "#C45B5B"
        else:  # Internal volume
            dist_type = 'volume'
            total_metric_label = "Total Volume"
            total_metric_col = 'volume_nm^3_per_mL_avg'
            total_metric_unit = "nm³/mL"
            metric_title = "Key metrics for volume distribution"
            focus_emoji = "🔷"
            focus_color = "#2C7F7F"
        
        # Display research focus prominently
        st.markdown(f"### {focus_emoji} {research_focus.split('(')[0].strip()}")
        
        # Display plot - with auto-refresh
        plot_placeholder = st.empty()
        
        if st.session_state.plot_output_dir:
            from pathlib import Path
            from PIL import Image
            
            output_path = Path(st.session_state.plot_output_dir)
            
            # Determine search pattern based on research focus
            if 'Number of particles' in research_focus:
                search_pattern = '*number*linear*.png'
            elif 'Surface area' in research_focus:
                search_pattern = '*surface*linear*.png'
            else:
                search_pattern = '*volume*linear*.png'
            
            # Find the plot file
            plot_files = sorted(output_path.glob(search_pattern))
            
            if plot_files:
                try:
                    img = Image.open(plot_files[0])
                    plot_placeholder.image(img, use_column_width=True)
                except Exception as e:
                    plot_placeholder.warning(f"Could not load plot image: {str(e)}")
            else:
                plot_placeholder.info("📊 Plot image is being generated... The visualization will appear here once plot generation completes. This typically takes a few seconds.\n\n💡 You can also view all plots in the **Plots** tab.")
                # Add a refresh button to manually check for plots
                if st.button("🔄 Check for Plot", key="refresh_plot"):
                    st.rerun()
        else:
            plot_placeholder.info("📊 Plots will be generated after analysis completes.")
        
        st.markdown("---")
        
        # Key metrics section with enhanced styling
        st.markdown(f"### 📈 Key Metrics")
        
        # Create metric cards
        col1, col2, col3 = st.columns(3, gap="large")
        
        # Total metric card
        with col1:
            st.markdown(f"""
            <div style="
                background: linear-gradient(135deg, {focus_color}22 0%, {focus_color}11 100%);
                border-left: 4px solid {focus_color};
                padding: 20px;
                border-radius: 8px;
                text-align: center;
            ">
                <div style="font-size: 14px; color: #666; margin-bottom: 10px; font-weight: 500;">
                    {total_metric_label}
                </div>
                <div style="font-size: 28px; font-weight: bold; color: {focus_color}; margin-bottom: 5px;">
            """, unsafe_allow_html=True)
            
            if total_metric_col in dist_df.columns:
                linear_df = dist_df[dist_df['scale'] == 'linear']
                total_val = linear_df[total_metric_col].sum() if not linear_df.empty else 0
                st.write(f"**{total_val:.2e}**")
            else:
                st.write("**N/A**")
            
            st.markdown(f"""
                </div>
                <div style="font-size: 12px; color: #999;">
                    {total_metric_unit}
                </div>
            </div>
            """, unsafe_allow_html=True)
        
        # D50 card
        with col2:
            d50_avg = "N/A"
            d50_ci = ""
            if 'linear' in stats and dist_type in stats['linear']:
                stat_dict = stats['linear'][dist_type]
                d50_val = stat_dict.get('D50_avg', 'N/A')
                if isinstance(d50_val, (int, float)):
                    d50_avg = f"{d50_val:.2f}"
                    d50_lower = stat_dict.get('D50_lower', 'N/A')
                    d50_upper = stat_dict.get('D50_upper', 'N/A')
                    if isinstance(d50_lower, (int, float)) and isinstance(d50_upper, (int, float)):
                        d50_ci = f"95% CI: {d50_lower:.2f} – {d50_upper:.2f}"
            
            st.markdown(f"""
            <div style="
                background: linear-gradient(135deg, #1f77b422 0%, #1f77b411 100%);
                border-left: 4px solid #1f77b4;
                padding: 20px;
                border-radius: 8px;
                text-align: center;
            ">
                <div style="font-size: 14px; color: #666; margin-bottom: 10px; font-weight: 500;">
                    Median Size (D50)
                </div>
                <div style="font-size: 28px; font-weight: bold; color: #1f77b4; margin-bottom: 5px;">
                    {d50_avg}
                </div>
                <div style="font-size: 11px; color: #999; line-height: 1.6;">
                    nm<br/>{d50_ci}
                </div>
            </div>
            """, unsafe_allow_html=True)
        
        # Span card
        with col3:
            span_avg = "N/A"
            span_ci = ""
            if 'linear' in stats and dist_type in stats['linear']:
                stat_dict = stats['linear'][dist_type]
                span_val = stat_dict.get('span_avg', 'N/A')
                if isinstance(span_val, (int, float)):
                    span_avg = f"{span_val:.2f}"
                    span_lower = stat_dict.get('span_lower', 'N/A')
                    span_upper = stat_dict.get('span_upper', 'N/A')
                    if isinstance(span_lower, (int, float)) and isinstance(span_upper, (int, float)):
                        span_ci = f"95% CI: {span_lower:.2f} – {span_upper:.2f}"
            
            st.markdown(f"""
            <div style="
                background: linear-gradient(135deg, #ff7f0e22 0%, #ff7f0e11 100%);
                border-left: 4px solid #ff7f0e;
                padding: 20px;
                border-radius: 8px;
                text-align: center;
            ">
                <div style="font-size: 14px; color: #666; margin-bottom: 10px; font-weight: 500;">
                    Distribution Width (Span)
                </div>
                <div style="font-size: 28px; font-weight: bold; color: #ff7f0e; margin-bottom: 5px;">
                    {span_avg}
                </div>
                <div style="font-size: 11px; color: #999; line-height: 1.6;">
                    (dimensionless)<br/>{span_ci}
                </div>
            </div>
            """, unsafe_allow_html=True)
        
        st.markdown("---")
        
        # Quality Control Section
        st.markdown("### ✅ Quality Control")
        
        metadata = st.session_state.metadata
        
        # Function to get QC status
        def get_qc_status(value, expected=None, check_type='equals'):
            if value == 'N/A' or value is None:
                return "⚠️ Unknown", "#ff9800"  # Orange
            
            if check_type == 'equals':
                if isinstance(value, str):
                    value_lower = str(value).lower()
                    if 'good' in value_lower or 'ok' in value_lower or 'pass' in value_lower:
                        return "✅ Good", "#4CAF50"  # Green
                    elif 'bad' in value_lower or 'fail' in value_lower or 'error' in value_lower:
                        return "❌ Failed", "#F44336"  # Red
                    else:
                        return f"⚠️ {value}", "#ff9800"  # Orange
            return "⚠️ Unknown", "#ff9800"
        
        # QC Metrics
        qc_col1, qc_col2, qc_col3, qc_col4 = st.columns(4)
        
        # Cell Check
        with qc_col1:
            cell_check = metadata.get('nta_cell_check_result', 'N/A')
            status, color = get_qc_status(cell_check)
            st.markdown(f"""
            <div style="background: {color}11; border-left: 4px solid {color}; padding: 15px; border-radius: 8px; text-align: center;">
                <div style="font-size: 12px; color: #666; margin-bottom: 8px;">Cell Check</div>
                <div style="font-size: 16px; font-weight: bold; color: {color};">{status}</div>
            </div>
            """, unsafe_allow_html=True)
        
        # Particle Drift Check
        with qc_col2:
            drift_check = metadata.get('nta_particle_drift_check_result', 'N/A')
            status, color = get_qc_status(drift_check)
            st.markdown(f"""
            <div style="background: {color}11; border-left: 4px solid {color}; padding: 15px; border-radius: 8px; text-align: center;">
                <div style="font-size: 12px; color: #666; margin-bottom: 8px;">Drift Check</div>
                <div style="font-size: 16px; font-weight: bold; color: {color};">{status}</div>
            </div>
            """, unsafe_allow_html=True)
        
        # Temperature Stability
        with qc_col3:
            temp = metadata.get('nta_temperature', 'N/A')
            if temp != 'N/A' and temp:
                try:
                    temp_val = float(str(temp).split('±')[0].strip())
                    temp_sd = float(str(temp).split('±')[1].strip()) if '±' in str(temp) else 0
                    if temp_sd < 1.0:
                        status, color = "✅ Stable", "#4CAF50"
                    elif temp_sd < 2.0:
                        status, color = "⚠️ Drift", "#ff9800"
                    else:
                        status, color = "❌ Unstable", "#F44336"
                except:
                    status, color = "⚠️ Unknown", "#ff9800"
            else:
                status, color = "⚠️ N/A", "#ff9800"
            
            st.markdown(f"""
            <div style="background: {color}11; border-left: 4px solid {color}; padding: 15px; border-radius: 8px; text-align: center;">
                <div style="font-size: 12px; color: #666; margin-bottom: 8px;">Temperature</div>
                <div style="font-size: 16px; font-weight: bold; color: {color};">{status}</div>
            </div>
            """, unsafe_allow_html=True)
        
        # Total Tracks (threshold: 500 tracks minimum for good quality)
        with qc_col4:
            total_tracks = metadata.get('nta_number_of_traces_sum', 'N/A')
            TRACKS_THRESHOLD = 500  # Minimum for good quality NTA analysis
            
            if total_tracks != 'N/A' and total_tracks:
                try:
                    tracks_val = int(float(total_tracks))
                    if tracks_val >= TRACKS_THRESHOLD:
                        status, color = "✅ Good", "#4CAF50"
                    elif tracks_val >= 200:
                        status, color = "⚠️ Marginal", "#ff9800"
                    else:
                        status, color = "❌ Low", "#F44336"
                except:
                    status, color = "⚠️ Unknown", "#ff9800"
            else:
                status, color = "⚠️ N/A", "#ff9800"
            
            st.markdown(f"""
            <div style="background: {color}11; border-left: 4px solid {color}; padding: 15px; border-radius: 8px; text-align: center;">
                <div style="font-size: 12px; color: #666; margin-bottom: 8px;">Tracks ({TRACKS_THRESHOLD}+ required)</div>
                <div style="font-size: 16px; font-weight: bold; color: {color};">{status}</div>
                <div style="font-size: 11px; color: #999; margin-top: 5px;">
                    {f'{tracks_val:,}' if total_tracks != 'N/A' and total_tracks else 'N/A'}
                </div>
            </div>
            """, unsafe_allow_html=True)
        
        st.markdown("---")
        
        # Sample Information Section
        st.markdown("### 📋 Sample Information")
        
        sample_col1, sample_col2, sample_col3 = st.columns(3)
        
        with sample_col1:
            sample = metadata.get('sample', 'N/A')
            st.markdown(f"""
            <div style="background: #f5f5f5; padding: 15px; border-radius: 8px;">
                <div style="font-size: 12px; color: #666; margin-bottom: 8px; font-weight: 500;">Sample ID</div>
                <div style="font-size: 14px; font-weight: bold; color: #333;">{sample}</div>
            </div>
            """, unsafe_allow_html=True)
        
        with sample_col2:
            electrolyte = metadata.get('nta_electrolyte', 'N/A')
            st.markdown(f"""
            <div style="background: #f5f5f5; padding: 15px; border-radius: 8px;">
                <div style="font-size: 12px; color: #666; margin-bottom: 8px; font-weight: 500;">Buffer/Electrolyte</div>
                <div style="font-size: 14px; font-weight: bold; color: #333;">{electrolyte}</div>
            </div>
            """, unsafe_allow_html=True)
        
        with sample_col3:
            temp_display = metadata.get('nta_temperature', 'N/A')
            st.markdown(f"""
            <div style="background: #f5f5f5; padding: 15px; border-radius: 8px;">
                <div style="font-size: 12px; color: #666; margin-bottom: 8px; font-weight: 500;">Temperature</div>
                <div style="font-size: 14px; font-weight: bold; color: #333;">{temp_display} °C</div>
            </div>
            """, unsafe_allow_html=True)
        
        st.markdown("---")
        
        # Additional Analysis Information
        st.markdown("### 📊 Analysis Details")
        
        detail_col1, detail_col2 = st.columns(2)
        
        with detail_col1:
            num_files = metadata.get('num_replicates', 'N/A')
            st.metric("Replicates Analyzed", num_files)
        
        with detail_col2:
            dilution = metadata.get('nta_dilution', 'N/A')
            st.metric("Dilution Factor", dilution)
    
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
        
        # Get research focus to auto-select appropriate plot
        research_focus = st.session_state.get('research_focus', 'Number of particles (size distribution)')
        if 'Number of particles' in research_focus:
            default_focus_plot = 'number'
        elif 'Surface area' in research_focus:
            default_focus_plot = 'surface_area'
        else:
            default_focus_plot = 'volume'
        
        st.markdown(f"**Currently focused on:** {research_focus.split('(')[0].strip()}")
        st.markdown("Select which plots to view (automatically selected plot is highlighted below):")
        
        if st.session_state.plot_output_dir is None or len(st.session_state.generated_plots) == 0:
            st.info("ℹ️ Plots will be generated automatically after analysis. Check back here after running the analysis!")
        else:
            # Define available plots
            plot_options = {
                'number': ('📊 Number-Weighted Distribution (Your Focus)', ['linear', 'logarithmic']),
                'volume': ('📊 Volume-Weighted Distribution', ['linear', 'logarithmic']),
                'surface_area': ('📊 Surface Area-Weighted Distribution', ['linear', 'logarithmic']),
            }
            
            # Update label for the focused plot
            for key in plot_options:
                if key == 'number':
                    if default_focus_plot == 'number':
                        plot_options[key] = ('📊 Number-Weighted Distribution ⭐ (Your Focus)', ['linear', 'logarithmic'])
                elif key == 'surface_area':
                    if default_focus_plot == 'surface_area':
                        plot_options[key] = ('📊 Surface Area-Weighted Distribution ⭐ (Your Focus)', ['linear', 'logarithmic'])
                elif key == 'volume':
                    if default_focus_plot == 'volume':
                        plot_options[key] = ('📊 Volume-Weighted Distribution ⭐ (Your Focus)', ['linear', 'logarithmic'])
            
            # Create columns for selection
            st.markdown("**Select plots to display:**")
            col1, col2, col3 = st.columns(3)
            selected_plots = {}
            
            for idx, (plot_key, (plot_label, scales)) in enumerate(plot_options.items()):
                col = [col1, col2, col3][idx % 3]
                with col:
                    # Auto-select the focused plot, others default to unchecked
                    is_default = (plot_key == default_focus_plot)
                    selected_plots[plot_key] = st.checkbox(plot_label, value=is_default)
            
            st.markdown("---")
            
            # Display selected plots
            from pathlib import Path
            from PIL import Image
            
            output_path = Path(st.session_state.plot_output_dir)
            
            displayed_count = 0
            for plot_key, (plot_label, scales) in plot_options.items():
                if selected_plots.get(plot_key, False):
                    # For plots with multiple scales, find all PNG files for this plot type
                    if plot_key == 'number':
                        search_pattern = '*number*.png'
                    elif plot_key == 'volume':
                        search_pattern = '*volume*.png'
                    elif plot_key == 'surface_area':
                        search_pattern = '*surface*.png'
                    else:
                        continue
                    
                    plot_files = sorted(output_path.glob(search_pattern))
                    
                    if plot_files:
                        with st.expander(plot_label, expanded=(displayed_count < 1)):
                            for plot_file in plot_files:
                                try:
                                    img = Image.open(plot_file)
                                    st.image(img, use_column_width=True, caption=plot_file.name)
                                except Exception as e:
                                    st.warning(f"Could not load image: {str(e)}")
                        displayed_count += 1
            
            if displayed_count == 0:
                st.info("No plots selected. Check the boxes above to view plots.")
            else:
                st.markdown(f"✅ Displaying {displayed_count} plot(s)")
                st.info("💡 Go to the **Download** tab to get PDFs, fit data, and ZIP package")
    
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
        
        # Download Section - Simplified to ZIP only
        st.subheader("📥 Download Your Results")
        
        if has_plots and DOWNLOAD_MANAGER_AVAILABLE:
            from pathlib import Path
            output_path = Path(st.session_state.plot_output_dir)
            pdf_files = list(output_path.glob('*Plot_*.pdf'))
            fit_files = list(output_path.glob('FitData_*.txt'))
            
            # Simple info box showing what's included
            st.info(f"""
            ### 📦 Complete Analysis Package
            
            Your ZIP file includes everything:
            - **Data Files:** Distribution (PSD) + Metadata with statistics
            - **Plots:** {len(pdf_files)} high-resolution PDFs (300 DPI, publication quality)
            - **Fit Data:** Lognormal fit parameters for {len(fit_files)} distributions
            
            All files use naming convention: `Data_{unique_id}_*`
            
            Organized folders inside ZIP:
            - `data/` - Distribution and metadata
            - `plots/` - All PDF plots
            - `fit_data/` - Fit parameters for re-plotting
            """)
            
            try:
                zip_bytes = nta_download_manager.create_download_zip(
                    st.session_state.plot_output_dir,
                    unique_id,
                    num_replicates,
                    metadata_dict=st.session_state.metadata,
                    distribution_df=st.session_state.distribution_data,
                    statistics_dict=st.session_state.statistics
                )
                
                if zip_bytes:
                    st.download_button(
                        label="📦 Download Complete ZIP",
                        data=zip_bytes,
                        file_name=f"Data_{unique_id}_complete.zip",
                        mime="application/zip",
                        key="complete_zip",
                        use_container_width=True
                    )
            except Exception as e:
                st.warning(f"Could not create ZIP: {str(e)}")
        else:
            st.info("👆 Run analysis and generate plots first, then return here to download results")

else:
    # No analysis yet
    st.info("👆 Upload NTA files and click 'RUN ANALYSIS' to begin")
    st.markdown("""
    ## How to use:
    
    1. **Upload files** - Use the sidebar to upload one or more NTA .txt files
    2. **Select research focus** - Choose whether you're interested in number, surface area, or volume distributions
    3. **Enter optional metadata** - Optionally add experimenter, project, location, and PI information (not required)
    4. **Run analysis** - Click the 'RUN ANALYSIS' button to process your data
    5. **View results** - Summary tab shows your selected focus; check Statistics tab for all detailed data
    6. **Explore plots** - The Plots tab automatically shows your focused plot, but you can view all others too
    7. **Download** - Export results as TSV data files, PDF plots, and fit parameters
    8. **Package** - Download everything as ZIP
    
    **Note:** User metadata is completely optional - just check the box if you want to add it.
    """)
