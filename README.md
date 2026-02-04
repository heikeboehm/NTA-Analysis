# NTA Data Analysis - MVP Version 1.0

A web-based particle size distribution analyzer for ZetaView NTA instrument data.

## Features

âœ… **Single & Multi-file Analysis** - Upload one or multiple NTA files  
âœ… **Automatic Averaging** - Replicate analysis with uncertainty calculation  
âœ… **Linear Number-Weighted Plots** - Publication-ready visualizations  
âœ… **Full Data Export** - Metadata, statistics, distribution data, and plots  
âœ… **Web Interface** - Streamlit cloud deployment (publicly accessible)  

## Files Included

- **`streamlit_app.py`** - Main web application
- **`nta_analyzer_simple.py`** - Analysis engine (no dependencies on notebook)
- **`requirements.txt`** - Python package dependencies

## Installation

### Local Testing (Windows/Mac/Linux)

```bash
# Clone/download files to a directory
cd your_nta_app_directory

# Create virtual environment (optional but recommended)
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate

# Install dependencies
pip install -r requirements.txt

# Run the app
streamlit run streamlit_app.py
```

The app will open at `http://localhost:8501`

### Cloud Deployment (Streamlit Cloud - GitHub)

**This is the easiest way to deploy publicly:**

1. **Create a GitHub Repository**
   - New repo: `nta-analysis` (or similar)
   - Upload the 3 files above to the repo

2. **Deploy to Streamlit Cloud**
   - Go to https://streamlit.io/cloud
   - Click "New app"
   - Select your GitHub repo
   - Main file: `streamlit_app.py`
   - Click "Deploy"

3. **That's it!** Your app is now public
   - Share the URL with your lab: `https://share.streamlit.io/YourUsername/nta-analysis/streamlit_app.py`
   - Anyone can use it (no login needed)

## How to Use

### Basic Workflow

1. **Upload Files**
   - Single file: Standard analysis
   - Multiple files: Automatic averaging with uncertainty

2. **Configure Parameters**
   - **Dilution Factor**: Must match your sample dilution
   - **Metadata**: Experimenter, location, project (for record-keeping)

3. **Review Results**
   - Linear number-weighted distribution plot (main & cumulative)
   - Distribution data table
   - D50 (median size)

4. **Download**
   - ZIP archive with all outputs
   - OR individual file downloads

### Output Files

When you download results, you get:

- **`Data_[ID]_metadata.txt`** - Measurement parameters and metadata
- **`Data_[ID]_PSD.txt`** - Particle size distribution table (for Origin/GraphPad)
- **`Plot_[ID]_linear_number.png`** - Distribution plot (high-quality)
- **`Stats_[ID].txt`** - Key statistics (D50, etc.)

## File Format

**Requires**: ZetaView `.txt` output files  
**Format**: Standard ZetaView measurement export  
**Encoding**: Latin-1 (handled automatically)

## Current Limitations (MVP)

- Linear scale number-weighted plot only (other plots can be added later)
- Single plot type displayed (others available in next version)
- No advanced statistics filtering

## What's Next (Future Versions)

- Volume-weighted plots
- Surface area-weighted plots
- Log-scale visualizations
- Statistics filtering (D10, D90, span, etc.)
- Batch processing
- Data comparison tools

## Troubleshooting

### "File format not recognized"
- Ensure it's a direct ZetaView `.txt` export
- Check file is not corrupted

### "No data rows found"
- File may not have "Size Distribution" section
- Try a different file to confirm

### Plot looks strange
- Check dilution factor is correct
- May need to filter data range in future versions

## Technical Details

- **Language**: Python 3.8+
- **Web Framework**: Streamlit 1.31
- **Data Processing**: Pandas, NumPy
- **Visualization**: Matplotlib
- **Statistics**: SciPy

## Support

For issues or feature requests, contact the author or create an issue in the GitHub repo.

---

**Version**: 1.0 MVP  
**Author**: Adapted from Heike Boehm's NTA analysis workflow  
**Last Updated**: February 2026  
**License**: Open use for research
