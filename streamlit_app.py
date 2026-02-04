"""
NTA Analysis Streamlit App
Uses nta_analyzer_cells_01_04.py (Cells 01-06 working code)
"""

import streamlit as st
import os
import tempfile
import pandas as pd
from pathlib import Path

# Import your working analyzer
from nta_analyzer_cells_01_04 import NTAAnalyzer, CONFIG

st.set_page_config(page_title="NTA Analysis", layout="wide")
st.title("NTA Particle Size Analysis - Cells 01-06")

# ============================================================================
# SIDEBAR: Configuration
# ============================================================================
with st.sidebar:
    st.header("Configuration")
    experimenter = st.text_input("Experimenter", value=CONFIG["project_metadata"]["experimenter"])
    location = st.text_input("Location", value=CONFIG["project_metadata"]["location"])
    project = st.text_input("Project", value=CONFIG["project_metadata"]["project"])
    
    # Update config
    CONFIG["project_metadata"]["experimenter"] = experimenter
    CONFIG["project_metadata"]["location"] = location
    CONFIG["project_metadata"]["project"] = project

# ============================================================================
# MAIN: File Upload and Analysis
# ============================================================================
st.header("Upload NTA Files")
uploaded_files = st.file_uploader("Choose NTA .txt files", type="txt", accept_multiple_files=True)

if uploaded_files:
    st.write(f"📁 {len(uploaded_files)} file(s) selected")
    
    if st.button("▶️ Run Analysis", type="primary", use_container_width=True):
        with tempfile.TemporaryDirectory() as temp_dir:
            # Save uploaded files
            temp_files = []
            for uf in uploaded_files:
                temp_path = os.path.join(temp_dir, uf.name)
                with open(temp_path, 'wb') as f:
                    f.write(uf.getbuffer())
                temp_files.append(temp_path)
            
            st.info(f"Processing {len(temp_files)} files...")
            
            # Run analysis
            try:
                analyzer = NTAAnalyzer(config=CONFIG)
                results = analyzer.process(temp_files)
                
                st.success("✅ Analysis Complete!")
                
                # ====================================================================
                # TAB 1: Results Overview
                # ====================================================================
                tab1, tab2, tab3, tab4 = st.tabs(["Overview", "Data", "Metadata", "Download"])
                
                with tab1:
                    st.subheader("Analysis Summary")
                    col1, col2, col3 = st.columns(3)
                    
                    with col1:
                        particles = results['metadata'].get('nta_total_particles_per_mL', 'N/A')
                        st.metric("Total Particles/mL", particles)
                    
                    with col2:
                        volume = results['metadata'].get('nta_total_volume_uL_per_mL', 'N/A')
                        st.metric("Total Volume (µL/mL)", volume)
                    
                    with col3:
                        d50 = results['metadata'].get('nta_linear_number_d50', 'N/A')
                        st.metric("D50 (number-weighted)", d50)
                    
                    st.write(f"Files analyzed: {results['num_replicates']}")
                
                # ====================================================================
                # TAB 2: Distribution Data
                # ====================================================================
                with tab2:
                    st.subheader("Particle Size Distribution")
                    dist = results['distribution']
                    
                    linear = dist[dist['scale'] == 'linear']
                    log = dist[dist['scale'] == 'logarithmic']
                    
                    st.write(f"**Linear scale:** {len(linear)} bins")
                    st.write(f"**Logarithmic scale:** {len(log)} bins")
                    st.write(f"**Total columns:** {len(dist.columns)}")
                    
                    with st.expander("View Linear Scale (first 10 rows)"):
                        st.dataframe(linear.head(10), use_container_width=True)
                    
                    with st.expander("View Logarithmic Scale (first 10 rows)"):
                        st.dataframe(log.head(10), use_container_width=True)
                
                # ====================================================================
                # TAB 3: Metadata
                # ====================================================================
                with tab3:
                    st.subheader("Metadata")
                    metadata_df = pd.DataFrame(
                        [(k, v) for k, v in results['metadata'].items()],
                        columns=['Field', 'Value']
                    )
                    st.dataframe(metadata_df, use_container_width=True, hide_index=True)
                
                # ====================================================================
                # TAB 4: Download
                # ====================================================================
                with tab4:
                    st.subheader("Download Results")
                    
                    # Save outputs
                    output_dir = os.path.join(temp_dir, 'outputs')
                    os.makedirs(output_dir, exist_ok=True)
                    created_files = analyzer.save_outputs(output_dir)
                    
                    st.write(f"Generated {len(created_files)} files:")
                    
                    for filepath in sorted(created_files):
                        filename = os.path.basename(filepath)
                        with open(filepath, 'r') as f:
                            content = f.read()
                        
                        st.download_button(
                            label=f"📥 {filename}",
                            data=content,
                            file_name=filename,
                            mime="text/plain",
                            use_container_width=True
                        )
            
            except Exception as e:
                st.error(f"❌ Error during analysis: {str(e)}")
                import traceback
                st.code(traceback.format_exc())

else:
    st.info("👈 Upload NTA files to get started")

st.markdown("---")
st.markdown("**NTA Analysis** | Cells 01-06 | Using nta_analyzer_cells_01_04.py")
