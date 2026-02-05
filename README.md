# 🧪 NTA Data Analysis Tool

[![Streamlit App](https://static.streamlit.io/badges/streamlit_badge_black_white.svg)](https://nta-analysis.streamlit.app)

A modern, interactive web application for analyzing Nanoparticle Tracking Analysis (NTA) data with real-time processing, visualization, and comprehensive data export.

## 📊 Features

✅ **Easy File Upload** - Drag & drop NTA .txt files directly  
✅ **Multi-Tab Interface** - Organized views for different analyses  
✅ **Persistent Metadata** - Pre-filled fields that stay between sessions  
✅ **Real-Time Processing** - Automatic averaging and statistics  
✅ **Quality Control** - Built-in QC alerts and validation  
✅ **Comprehensive Statistics** - D10, D50, D90, span with uncertainty bounds  
✅ **Publication-Quality Plots** - Number, volume, and surface area weighted distributions  
✅ **Complete Data Export** - Download as CSV and PNG  

## 🚀 Quick Start

### Option 1: Use the Deployed App (Easiest)

Visit: **[https://nta-analysis.streamlit.app](https://nta-analysis.streamlit.app)**

No installation needed! Just upload your files and analyze.

### Option 2: Run Locally

#### Prerequisites
- Python 3.8+ (3.11 recommended)
- pip (Python package manager)

#### Install & Run

```bash
# Clone the repository
git clone https://github.com/yourusername/NTA-Analysis.git
cd NTA-Analysis

# Install dependencies
pip install -r requirements.txt

# Run the app
streamlit run streamlit_app.py
```

The app will open at `http://localhost:8501`

## 📋 Application Tabs

### 📊 Summary
Quick overview with key metrics:
- Number of files processed
- Total detected particles/tracks
- Concentration (particles/mL)
- Median size (D50)
- Quality control status
- Sample metadata
- Quick visualization

### 📈 Distributions
Detailed distribution data tables:
- All size bins with measurements
- Uncertainties (standard deviations)
- Linear and logarithmic scales
- Filter and download as CSV

### 📉 Statistics
Statistical summary:
- D10, D50, D90 percentile sizes
- Distribution span
- Separate rows for each distribution type
- Uncertainty bounds
- Download as CSV

### 🎨 Plots
Generate publication-quality visualizations:
- Number-weighted distribution
- Volume-weighted distribution
- Surface area-weighted distribution
- Linear and logarithmic scales
- Error bars and D-value markers
- Download as PNG (300 DPI)

### 📋 Metadata
Complete measurement parameters:
- All extracted metadata fields
- Measurement conditions
- Instrument settings
- Quality control notes
- Download as CSV

### ⬇️ Download
Export all results:
- Distribution data (CSV)
- Statistics (CSV)
- Metadata (CSV)
- All with timestamps and traceability

## 📁 Project Structure

```
NTA-Analysis/
├── streamlit_app.py          # Main Streamlit application
├── nta_analysis.py           # Core analysis module (57 functions)
├── requirements.txt          # Python dependencies
├── .streamlit/
│   └── config.toml          # Streamlit Cloud configuration
├── .gitignore               # Git ignore rules
└── README.md                # This file
```

## 📖 Documentation

See the documentation files in the outputs folder:

- **PROJECT_SUMMARY.md** - Complete project overview
- **SETUP_GUIDE.md** - Detailed installation and usage guide
- **TRANSFORMATION_SUMMARY.md** - Technical details

## 🔧 Configuration

### Default Metadata

Edit `nta_analysis.py` to customize default values:

```python
CONFIG['project_metadata'] = {
    "experimenter": "Your_Initials",
    "location": "Your_Lab_Location",
    "project": "Your_Project_Name",
    "pi": "Principal_Investigator_Initials",
}
```

These values will be pre-filled in the app and persist between analyses.

### Streamlit Theme

Customize the app appearance by editing `.streamlit/config.toml`:

```toml
[theme]
primaryColor = "#1f77b4"
backgroundColor = "#ffffff"
secondaryBackgroundColor = "#f0f2f6"
textColor = "#262730"
font = "sans serif"
```

## 📦 Dependencies

All dependencies are in `requirements.txt`:

```
streamlit>=1.28.0
pandas>=1.5.0
numpy>=1.23.0
matplotlib>=3.6.0
scipy>=1.9.0
scikit-learn>=1.2.0
Pillow>=9.0.0
```

Install all at once:
```bash
pip install -r requirements.txt
```

## 🌐 Deployment to Streamlit Cloud

### Step 1: Prepare Repository

Make sure your repo has:
- ✅ `streamlit_app.py` (entry point)
- ✅ `requirements.txt` (dependencies)
- ✅ `nta_analysis.py` (core module)
- ✅ `.streamlit/config.toml` (optional but recommended)
- ✅ `.gitignore` (to keep repo clean)

### Step 2: Push to GitHub

```bash
git add .
git commit -m "Initial Streamlit app setup"
git push origin main
```

### Step 3: Deploy to Streamlit Cloud

1. Go to [https://share.streamlit.io](https://share.streamlit.io)
2. Click "Deploy an app"
3. Fill in:
   - GitHub repo: `yourusername/NTA-Analysis`
   - Branch: `main`
   - File path: `streamlit_app.py`
4. Click "Deploy!"

Your app will be live at: `https://nta-analysis.streamlit.app`

### Step 4: Automatic Updates

Every time you push to GitHub, Streamlit Cloud automatically redeploys your app! 🚀

```bash
# Make a change locally
# Commit and push
git push origin main

# Your app updates automatically within seconds!
```

## 🔄 Workflow

### Local Development
```bash
# Make changes to code
# Test locally
streamlit run streamlit_app.py

# Push to GitHub when ready
git add .
git commit -m "Your message"
git push origin main

# App auto-deploys to streamlit.io
```

### For Collaborators
```bash
# Clone the repo
git clone https://github.com/yourusername/NTA-Analysis.git

# Create your branch
git checkout -b feature/your-feature

# Make changes, test, commit
streamlit run streamlit_app.py
git add .
git commit -m "Add your feature"

# Push and create a Pull Request
git push origin feature/your-feature
```

## 📊 Data Privacy & Security

✅ **All processing is local** - Data never leaves your machine (when running locally)  
✅ **Streamlit Cloud** - Data transmitted securely, never stored on Streamlit servers  
✅ **No account required** - Works immediately upon upload  
✅ **Open source** - All code is visible and auditable  
✅ **Reproducible** - Same input always produces same output  

## 🐛 Troubleshooting

### App won't start
```bash
# Verify all dependencies installed
pip install -r requirements.txt --upgrade

# Check Python version (3.8+)
python --version

# Try running with verbose output
streamlit run streamlit_app.py --logger.level=debug
```

### Files won't upload
- Ensure files have `.txt` extension
- Check file size (should be > 100 KB)
- Verify file is valid ZetaView output

### Streamlit Cloud deployment issues
- Check repo visibility (should be public)
- Verify `streamlit_app.py` is in root directory
- Ensure `requirements.txt` is up-to-date
- Check the Streamlit Cloud logs for errors

## 📚 Documentation

For detailed information, see:

- **Installation**: SETUP_GUIDE.md
- **Features**: PROJECT_SUMMARY.md
- **Technical Details**: TRANSFORMATION_SUMMARY.md

## 🤝 Contributing

To contribute improvements:

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/improvement`)
3. Make your changes
4. Commit your work (`git commit -am 'Add improvement'`)
5. Push to GitHub (`git push origin feature/improvement`)
6. Open a Pull Request

## 📝 Development Notes

### Original Framework
- Based on analysis code by: Heike Boehm (MPI for Medical Research)
- Reference: ZetaView Nanoparticle Tracking Analyzer documentation

### Transformation
- Original format: Jupyter Notebook (NTA_analysis_average.ipynb)
- Transformed to: Python module + Streamlit web app
- Date: February 2026
- Status: Production Ready

### Future Enhancements
Planned features for future versions:
- [ ] Batch sample comparison
- [ ] PDF report generation
- [ ] Historical data storage
- [ ] Parameter templates/presets
- [ ] Advanced QC dashboards

## 📜 License

This application is derived from the original NTA analysis framework.
Use according to your institution's policies.

## 📞 Support

### Getting Help
1. Check the Troubleshooting section above
2. Review documentation files
3. Check browser console for errors (F12)
4. Verify dependencies: `pip list | grep -E "streamlit|pandas|numpy"`

### Common Issues

**Q: How long does analysis take?**  
A: Usually 10-30 seconds for typical samples, depending on file count and size.

**Q: Can I use this offline?**  
A: Yes! When running locally with `streamlit run`, everything is offline.

**Q: How do I modify the analysis?**  
A: Edit `nta_analysis.py` - all functions are documented and modular.

**Q: Can I add new features?**  
A: Absolutely! Edit `streamlit_app.py` to change the UI or add tabs.

## 🎉 Quick Links

- 🌐 **Live App**: https://nta-analysis.streamlit.app
- 📦 **GitHub**: https://github.com/yourusername/NTA-Analysis
- 📚 **Streamlit Docs**: https://docs.streamlit.io
- 🐍 **Python**: https://www.python.org

## 📊 Statistics

- **Code**: 4,900+ lines
- **Functions**: 57 analysis functions
- **Documentation**: 3 comprehensive guides
- **Browser Support**: Chrome, Firefox, Safari, Edge
- **Deployment**: Fully automated via GitHub

---

**Made with ❤️ for better NTA data analysis**

Version: 1.0.0 | Status: ✅ Production Ready | Updated: February 2026
