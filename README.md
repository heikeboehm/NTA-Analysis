# 🧪 NTA Data Analysis Tool

Fast, reliable analysis of Nanoparticle Tracking Analysis data. Upload your ZetaView files and get instant results with quality control checks.

## 🚀 Quick Start

### Online (No Installation)
Go to: **https://nta-analysis.streamlit.app** and start analyzing.

### Local Installation
```bash
git clone https://github.com/heikeboehm/NTA-Analysis.git
cd NTA-Analysis
pip install -r requirements.txt
streamlit run streamlit_app.py
```

Opens at `http://localhost:8501`

---

## 📊 How to Use

**1. Upload Files**  
Drag & drop your NTA `.txt` files from ZetaView.

**2. Select Research Focus**  
Choose what you want to analyze: particle number, surface area, or internal volume.

**3. Add Sample Info (Optional)**  
Enter sample ID, buffer, temperature, dilution if you want it saved.

**4. Run Analysis**  
Click "RUN ANALYSIS" and wait ~30 seconds.

**5. View Results**

- **Summary**: Quick overview with quality control checks
- **Statistics**: D10, D50, D90, span with confidence intervals  
- **Plots**: Publication-quality distributions (auto-focused on your selection)
- **Download**: Export everything as CSV and PNG

---

## ✅ Quality Control

The Summary tab shows you:
- ✅ Cell check result
- ✅ Particle drift check  
- ✅ Temperature stability (< 0.5°C is good)
- ✅ Number of tracks detected (500+ is good)
- ✅ Number of replicates (3+ is good)

Green checkmark = all good. Warning symbol = check the details.

---

## 📋 What You Get

- Particle size distribution analysis
- Number, volume, and surface area weighted data
- Uncertainty estimates for all metrics
- Automated quality control
- Publication-ready plots
- CSV exports for further analysis

---

## 🛠️ Settings

### Customize Default Metadata
Edit `nta_analysis.py`:
```python
CONFIG['project_metadata'] = {
    "experimenter": "Your_Initials",
    "location": "Your_Lab",
    "project": "Your_Project",
    "pi": "Your_Initials",
}
```

### Change Theme
Edit `.streamlit/config.toml` to customize colors and appearance.

---

## ❓ FAQ

**How long does it take?**  
Usually 20-30 seconds per analysis.

**Can I use it offline?**  
Yes - run it locally with `streamlit run streamlit_app.py`.

**What file format do I need?**  
ZetaView `.txt` export files. One or multiple replicates.

**Can I change the analysis?**  
Yes - `nta_analysis.py` has all the functions. Everything is documented.

**How do I report bugs?**  
Open an issue on GitHub or contact the maintainer.

---

## 📦 Files You Need

- `streamlit_app.py` - The web interface
- `nta_analysis.py` - Analysis functions
- `requirements.txt` - Python dependencies

---

## 📚 More Info

- **Detailed Setup**: SETUP_GUIDE.md
- **How It Works**: PROJECT_SUMMARY.md  
- **Technical Details**: TRANSFORMATION_SUMMARY.md

---

**Questions?** Check the documentation or review your data file format.

**Version**: 1.0.0 | **Status**: Production Ready | **Updated**: February 2026
