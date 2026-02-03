"""
NTA Data Analysis - Streamlit Web App

Interactive web interface for particle size distribution analysis.
Users can upload files, adjust parameters, and download results.
"""

import streamlit as st
import os
import tempfile
import shutil
import zipfile
from io import BytesIO

from nta_analyzer_simple import NTAAnalyzer

# ============================================================================
# PAGE CONFIGURATION
# ============================================================================

st.set_page_config(
    page_title="NTA Data Analysis",
    page_icon="🔬",
    layout="wide",
    initial_sidebar_state="expanded"
)

st.title("🔬 NTA Particle Size Distribution Analysis")
st.markdown("""
Analyze Nanoparticle Tracking Analysis (NTA) data from ZetaView instruments.
Upload your data files and get instant analysis with plots and statistics.
""")

# ============================================================================
# SIDEBAR - CONFIGURATION
# ============================================================================

st.sidebar.header("Configuration Parameters")

with st.sidebar.expander("📋 Metadata Fields", expanded=True):
    experimenter = st.text_input(
        "Experimenter Initials",
        value="HB",
        help="Your name or initials"
    )
    
    location = st.text_input(
        "Lab Location",
        value="MPI-CBP",
        help="Where the analysis was performed"
    )
    
    project = st.text_input(
        "Project Name",
        value="NTA_Analysis",
        help="Project or sample name"
    )

with st.sidebar.expander("⚙️ Analysis Parameters", expanded=True):
    dilution = st.number_input(
        "Dilution Factor",
        value=1000.0,
        min_value=1.0,
        help="Sample dilution factor applied during measurement"
    )

# ============================================================================
# MAIN CONTENT - FILE UPLOAD
# ============================================================================

col1, col2 = st.columns([2, 1])

with col1:
    st.subheader("📁 Upload NTA Data Files")
    st.info("ℹ️ You can upload one or multiple .txt files from ZetaView for replicate analysis")
    
    uploaded_files = st.file_uploader(
        "Choose NTA data files",
        type="txt",
        accept_multiple_files=True,
        key="file_uploader"
    )

with col2:
    st.subheader("📊 Analysis Type")
    st.metric("Files Selected", len(uploaded_files) if uploaded_files else 0)

# ============================================================================
# PROCESSING & ANALYSIS
# ============================================================================

if uploaded_files:
    
    # Create temporary directory for uploaded files
    temp_dir = tempfile.mkdtemp()
    temp_files = []
    
    try:
        # Save uploaded files
        for uploaded_file in uploaded_files:
            file_path = os.path.join(temp_dir, uploaded_file.name)
            with open(file_path, 'wb') as f:
                f.write(uploaded_file.getbuffer())
            temp_files.append(file_path)
        
        # Show processing status
        with st.spinner("🔄 Processing files..."):
            
            # Initialize analyzer
            analyzer = NTAAnalyzer()
            
            # Prepare metadata overrides
            metadata_overrides = {
                'experimenter': experimenter,
                'location': location,
                'project': project,
                'dilution_factor': dilution
            }
            
            # Process files
            results = analyzer.process(temp_files, metadata_overrides)
            
            # Generate plot to get D-values
            fig = analyzer.get_plot()
            
            st.success("✅ Analysis complete!")
        
        # ====================================================================
        # RESULTS DISPLAY
        # ====================================================================
        
        st.divider()
        st.subheader("📈 Results")
        
        # Display metadata
        with st.expander("📋 Analysis Information", expanded=False):
            col1, col2, col3, col4 = st.columns(4)
            
            with col1:
                st.metric("Sample", results['metadata'].get('sample_name', 'N/A'))
            
            with col2:
                st.metric("Files Analyzed", results['metadata'].get('num_files', 1))
            
            with col3:
                st.metric("Dilution Factor", f"{results['metadata'].get('dilution_factor', 1):.1f}")
            
            with col4:
                st.metric("Temperature (°C)", f"{results['metadata'].get('temperature', 'N/A')}")
            
            st.write("**Full Metadata:**")
            for key, value in results['metadata'].items():
                st.write(f"• **{key}:** {value}")
        
        # Generate and display plot
        st.subheader("📊 Linear Number-Weighted Distribution")
        
        try:
            fig = analyzer.get_plot()
            st.pyplot(fig)
            
            # Display D-values from analyzer results
            if 'd_values' in analyzer.results:
                st.divider()
                col1, col2, col3 = st.columns(3)
                
                d_values = analyzer.results['d_values']
                
                with col1:
                    st.metric("D10 (nm)", f"{d_values.get('D10', 0):.1f}")
                
                with col2:
                    st.metric("D50 (nm)", f"{d_values.get('D50', 0):.1f}")
                
                with col3:
                    st.metric("D90 (nm)", f"{d_values.get('D90', 0):.1f}")
                
                st.caption("D-values represent the size where 10%, 50%, and 90% of particles are smaller")
        
        except Exception as e:
            st.error(f"Error generating plot: {str(e)}")
        
        # Display distribution data table
        with st.expander("📊 Distribution Data Table", expanded=False):
            st.dataframe(
                results['distribution'],
                use_container_width=True,
                height=400
            )
        
        # ====================================================================
        # DOWNLOAD RESULTS
        # ====================================================================
        
        st.divider()
        st.subheader("⬇️ Download Results")
        
        # Create output directory
        output_dir = tempfile.mkdtemp()
        
        try:
            # Save all outputs
            created_files = analyzer.save_outputs(output_dir)
            
            # Create zip file
            zip_buffer = BytesIO()
            unique_id = results['metadata'].get('persistentID', 'analysis')
            
            with zipfile.ZipFile(zip_buffer, 'w', zipfile.ZIP_DEFLATED) as zip_file:
                for file_path in created_files:
                    arcname = os.path.basename(file_path)
                    zip_file.write(file_path, arcname=arcname)
            
            zip_buffer.seek(0)
            
            # Download button
            st.download_button(
                label="📦 Download All Results (ZIP)",
                data=zip_buffer.getvalue(),
                file_name=f"NTA_analysis_{unique_id}.zip",
                mime="application/zip",
                help="Contains metadata, distribution data, plots, and statistics"
            )
            
            # Individual file downloads
            st.write("**Or download individual files:**")
            
            col1, col2, col3 = st.columns(3)
            
            for file_path in created_files:
                filename = os.path.basename(file_path)
                
                with open(file_path, 'rb') as f:
                    file_data = f.read()
                
                if 'metadata' in filename:
                    with col1:
                        st.download_button(
                            label="📄 Metadata",
                            data=file_data,
                            file_name=filename,
                            mime="text/plain"
                        )
                elif 'Stats' in filename:
                    with col2:
                        st.download_button(
                            label="📊 Statistics",
                            data=file_data,
                            file_name=filename,
                            mime="text/plain"
                        )
                elif 'png' in filename or 'pdf' in filename:
                    with col3:
                        mime_type = "image/png" if 'png' in filename else "application/pdf"
                        st.download_button(
                            label="📈 Plot",
                            data=file_data,
                            file_name=filename,
                            mime=mime_type
                        )
        
        except Exception as e:
            st.error(f"Error preparing downloads: {str(e)}")
        
        finally:
            # Cleanup
            if os.path.exists(output_dir):
                shutil.rmtree(output_dir)
    
    except Exception as e:
        st.error(f"❌ Analysis failed: {str(e)}")
        st.info("Please check your files and try again.")
    
    finally:
        # Cleanup temporary directory
        if os.path.exists(temp_dir):
            shutil.rmtree(temp_dir)

else:
    st.info("👆 Start by uploading one or more NTA data files above")
    
    # Show example/instructions
    with st.expander("📖 How to use this app", expanded=False):
        st.markdown("""
        ### Quick Start
        
        1. **Upload Files:** Click the file uploader to select your ZetaView .txt files
           - Single file: Straightforward analysis
           - Multiple files: Automatic averaging with uncertainty calculation
        
        2. **Configure Parameters:** Adjust metadata and dilution factor in the sidebar
           - **Dilution Factor:** Must match your sample dilution
           - **Metadata:** For record keeping (experimenter, location, project)
        
        3. **Review Results:** 
           - Check the linear number-weighted distribution plot
           - Examine the data table
           - Read the statistics
        
        4. **Download:** Get all results as a ZIP or individual files
           - Metadata (TXT)
           - Distribution data (TXT) - ready for plotting in Origin/GraphPad
           - Plot (PNG)
           - Statistics (TXT)
        
        ### File Format Requirements
        - Format: ZetaView .txt output files
        - Encoding: Latin-1 (standard for ZetaView)
        - Each file must contain a "Size Distribution" section
        """)

# ============================================================================
# FOOTER
# ============================================================================

st.divider()

col1, col2, col3 = st.columns(3)

with col1:
    st.caption("🔬 NTA Analysis Tool")

with col2:
    st.caption("Powered by Streamlit")

with col3:
    st.caption("Version 1.0 MVP")

st.caption("Built for analyzing particle size distributions from ZetaView NTA instruments")
