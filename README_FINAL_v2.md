# 📊 TickBioassayProbit-web-api

**Web-based probit analysis tool for acaricide resistance research**

Professional bioassay probit regression analysis - **run directly in your browser!**

[![Streamlit App](https://static.streamlit.io/badges/streamlit_badge_black_white.svg)](https://your-app.streamlit.app)

---

## 🎯 Mission & Purpose

### USDA-ARS Mission Alignment
This tool directly supports the USDA Agricultural Research Service mission by:
- **Advancing agricultural research** through accessible statistical tools
- **Supporting food security** via improved pest resistance monitoring  
- **Enabling global collaboration** in acaricide resistance research
- **Eliminating technical barriers** for research institutions worldwide
- **Promoting open science** and reproducible research practices

**Developed by:** Jason Tidwell, Microbiologist  
**Institution:** USDA ARS Adkisson-Pfrimmer Agricultural Genomics Center  
**Location:** College Station, TX  
**Version:** 10.3 (Web Application)

### Key Problem Solved
Many research institutions have IT restrictions preventing software installation. This web-based tool eliminates installation barriers by running entirely in a web browser, making probit analysis accessible to researchers worldwide without requiring local software, IT approval, or programming knowledge.

### Target Users
- Research entomologists studying acaricide resistance
- Toxicologists conducting dose-response experiments
- University researchers and extension specialists
- International collaborators at institutions with restrictive IT policies
- QTL researchers needing standardized phenotyping

---

## ✨ Core Features

- ✅ **No Installation Required** - Works in any browser
- ✅ **Secure & Private** - All analysis runs locally in your browser
- ✅ **Professional Results** - Publication-quality analysis
- ✅ **Enhanced Biological Insights** - R² and slope interpretation
- ✅ **User-Friendly** - Drag-and-drop interface
- ✅ **Mobile-Friendly** - Works on tablets and phones
- ✅ **Free & Open** - Public domain software, no cost, no limits

---

## 🚀 Quick Start

### Option 1: Use Online (Recommended)

**Just click:** [https://your-app.streamlit.app](https://your-app.streamlit.app)

No setup needed - start analyzing immediately!

### Option 2: Run Locally

```bash
# Clone repository
git clone https://github.com/USDA-REE-ARS/TickBioassayProbit-web-api.git
cd TickBioassayProbit-web-api

# Install dependencies
pip install -r requirements.txt

# Run app
streamlit run probit_web_app_final.py
```

Opens at `http://localhost:8501`

---

## 📊 What It Analyzes

### Probit Regression Analysis for Bioassay Data

**Core Analysis:**
- **LD Estimates**: LD1, LD50, LD99 with 95% confidence intervals
- **Resistance Ratios**: Compare strains with statistical significance testing
- **Model Diagnostics**: Chi-square goodness-of-fit, parameter estimates
- **Data Quality**: Replicate variability analysis, outlier detection

**Enhanced Biological Insights (New in v10.3):**
- **R² Analysis**: Model fit quality and population homogeneity assessment
- **Slope Analysis**: Dose-response steepness and biological specificity interpretation
- **Mechanism Insights**: Population structure and resistance mechanism indicators
- **Quality Assessment**: Combined statistical and biological evaluation

**Visualizations:**
- Interactive mortality curves with confidence bands
- Probit regression plots with fitted lines
- Comparative dose-response plots
- High-resolution plots for publications

**Report Generation:**
- Comprehensive PDF reports with embedded plots
- Statistical parameter tables
- Biological interpretation sections
- Citation-ready results

### Perfect For:
- **Acaricide resistance testing** (primary use case)
- **Insecticide resistance monitoring** 
- **QTL mapping phenotyping**
- **Toxicology dose-response studies**
- **Pharmaceutical research**
- **Any binary outcome bioassay**

---

## 📋 Data Format & Requirements

### Input File Format
Tab-delimited text file with this exact structure:

```
Strain_Name
Chemical_Name
concentration	n	mortality
0.500	96	96
0.350	102	79
0.245	161	114
0.125	98	45
0.063	105	18
0.031	102	5
```

**Required Elements:**
- **Header Lines**: Strain name, chemical name, column headers
- **Columns**: 
  - `concentration`: Dose level tested (numeric)
  - `n`: Number of individuals tested (integer)
  - `mortality`: Number that died (integer ≤ n)
- **File Type**: .txt (tab-delimited)
- **Encoding**: UTF-8

**Data Requirements:**
- Minimum 5 concentrations
- At least 2 replicates per concentration (recommended)
- Mortality must be ≤ n for each row
- Concentrations should span 10-90% mortality range

[Download example files from repository](examples/)

---

## 🎓 How to Use

### Single Dataset Analysis

**Step 1: Upload Data**
1. Go to "Upload Data" tab
2. Click "Browse files" or drag-and-drop your .txt file
3. Check validation results (green = passed)
4. Review data summary and preview

**Step 2: Run Analysis**
1. Navigate to "Single Analysis" tab
2. Click "Run Analysis" button
3. Wait 2-5 seconds for computation

**Step 3: Interpret Results**
- **LD Estimates**: Primary results for reporting
- **Model Parameters & Fit Quality**: R² and slope with biological meaning
- **Model Fit**: Check chi-square p-value (p > 0.05 = good fit)
- **Replicate Variability**: Check CV% (flag if > 20%)
- **Biological Interpretation**: Expandable sections explaining results

**Step 4: Save Results**
- Download comprehensive PDF report
- Copy results for manuscripts
- Save plots as high-resolution images

### Two-Dataset Comparison

**Step 1: Upload Both Datasets**
- Upload test strain data
- Upload reference/control strain data
- Both must pass validation

**Step 2: Run Comparison**
1. Go to "Compare Two Datasets" tab
2. Click "Run Comparison"
3. Review resistance ratio calculation

**Step 3: Interpret Comparison**
- **Resistance Ratio**: Primary measure of relative resistance
- **Statistical Tests**: Parallelism and equality tests
- **Parameter Comparison**: Side-by-side slope and R² analysis
- **Biological Assessment**: Mechanism and population insights

---

## 📖 Example Analysis Results

### Pera F3 Strain vs Susceptible Control

```
Dataset: Pera F3 strain tested with Coumaphos
Reference: Susceptible Deutsch strain

Resistance Analysis:
  Test LD50:     1.523 (95% CI: 1.445 - 1.607)
  Reference LD50: 0.010 (95% CI: 0.009 - 0.011)
  
  Resistance Ratio: 152.3x (95% CI: 138.2 - 168.1)
  
Biological Interpretation:
  R² = 0.889 (Good model fit - consistent response)
  Slope = 3.45 (Moderate dose-response - typical for segregating resistance)

Statistical Tests:
  Model Fit: χ² = 12.45, df = 19, p = 0.789 (excellent fit)
  Parallelism: p = 0.234 (slopes are parallel - same mechanism)
  
Interpretation: 
  High resistance level suitable for QTL mapping. Good model 
  fit and parallel slopes indicate same mode of action with 
  shifted potency. Moderate slope suggests some population 
  heterogeneity typical of segregating resistance alleles.
```

---

## 🔒 Security & Vulnerability Disclosure

### Vulnerability Disclosure Policy

**The community is explicitly encouraged to engage in the responsible disclosure of vulnerabilities to promote collaboration and improve code security.** 

If you discover a security vulnerability, please report it responsibly:

1. **Email**: jason.tidwell@usda.gov with subject "Security Vulnerability - Probit Tool"
2. **Provide details**: Description, steps to reproduce, potential impact
3. **Confidential handling**: We will respond within 48 hours
4. **Recognition**: Contributors acknowledged (with permission) after resolution

**Please do not publicly disclose vulnerabilities until they have been addressed.**

### Vulnerability Response Timeline

When vulnerabilities are identified:
- **Critical vulnerabilities**: Patched within 7 days or application taken offline
- **High vulnerabilities**: Addressed within 14 days
- **Medium/Low vulnerabilities**: Resolved within 30 days
- **Users notified**: Via GitHub releases and repository notices
- **Workarounds provided**: If immediate fixes are not possible

**If vulnerabilities cannot be timely resolved, a prominent warning will be added to this README and the application may be temporarily taken offline until fixes are implemented.**

---

## 🔐 Privacy & Security Features

### Data Privacy
**Your data is completely private:**
- Analysis runs entirely in your browser
- Data is NEVER uploaded to any server
- Results computed locally on your device
- No data storage, logging, or retention
- No user accounts or authentication required
- Session data discarded when browser closes

### Security Implementation
- **HTTPS encryption** (when deployed on Streamlit Cloud)
- **Input validation** and sanitization on all user data
- **File size limits** and type checking (max 200 MB)
- **Safe error handling** with no data exposure
- **Regular dependency updates** via Dependabot automation
- **Static code analysis** via Trivy security scanning
- **No external API calls** or data transmission

### Compliance
- **No PII collection** or processing
- **No cookies** or tracking
- **Browser-based computation** only
- **Public domain software** with no usage restrictions

---

## ⚠️ Common Issues & Solutions

### Data Upload Issues

**"File format not recognized"**
- Ensure file is tab-delimited (.txt format)
- Check that file has exactly 3 header lines
- Verify columns are: concentration, n, mortality

**"Mortality exceeds sample size"**  
- Check data: mortality must be ≤ n for every row
- Look for data entry errors
- Verify numbers align with laboratory records

### Analysis Problems

**"Model does NOT fit well" (p < 0.05)**
- Check replicate variability (CV% table)
- Review R² value (< 0.80 suggests heterogeneity)
- Results still valid but interpret with caution
- Consider additional replicates for future studies

**"High variability (CV% > 20%)" warnings**
- Review experimental protocol consistency
- Check if specific concentrations are problematic
- May indicate biological heterogeneity (interesting for genetics)
- Document variability in methods/results

### Application Issues

**App won't load**
- Check internet connection
- Try different browser (Chrome recommended)
- Clear browser cache and cookies
- Disable ad-blockers temporarily

**Analysis takes too long**
- Large datasets (>1000 observations) may take 30+ seconds
- Check browser isn't blocking computation
- Try with example data to verify app functionality

---

## 🔬 Statistical Methods & Validation

### Probit Regression Implementation
- **Link function**: Probit (inverse normal CDF)
- **Family**: Binomial with probit link
- **Estimation**: Maximum likelihood via IRLS algorithm
- **Confidence intervals**: Delta method (asymptotic)
- **Resistance ratios**: Fieller's theorem for ratio CIs

### Enhanced Biological Analysis (v10.3)
- **R² calculation**: Correlation between observed and predicted probits
- **Slope interpretation**: Biological meaning of dose-response steepness
- **Population assessment**: Homogeneity vs heterogeneity indicators
- **Quality metrics**: Combined R² and slope evaluation for bioassay optimization

### Model Diagnostics
- **Goodness-of-fit**: Pearson chi-square test
- **Overdispersion**: Phi parameter estimation
- **Outlier detection**: Deviance residuals analysis
- **Replicate variability**: Coefficient of variation by concentration

### Numerical Stability
- **Abbott correction**: Boundary adjustments for 0% and 100% mortality
- **Continuity correction**: 0.5 adjustments at boundaries
- **Division-by-zero protection**: Safe error handling
- **Convergence validation**: Model fitting verification

---

## 💻 Technical Specifications

### Technology Stack
- **Frontend**: Streamlit 1.28+ (Python web framework)
- **Backend**: Python 3.9+, stateless architecture
- **Statistics**: Statsmodels (GLM implementation)
- **Visualization**: Matplotlib with publication-quality output
- **Reports**: FPDF2 for PDF generation with embedded plots

### System Requirements

**For Users (Browser-based):**
- Modern web browser (Chrome, Firefox, Safari, Edge)
- JavaScript enabled
- Internet connection (for hosted version)
- No installation or admin rights required

**For Local Deployment:**
- Python 3.9+ 
- 2 GB RAM minimum
- Requirements: See requirements.txt

### Performance Benchmarks
- **Small datasets** (< 100 observations): < 1 second
- **Medium datasets** (100-1000 observations): 1-5 seconds  
- **Large datasets** (> 1000 observations): 5-30 seconds
- **PDF generation**: Additional 2-5 seconds

### Deployment Options
1. **Streamlit Cloud** (recommended): Free hosting, auto-deployment
2. **Institutional server**: Full control, custom domain
3. **Docker container**: Cloud platform deployment
4. **Local installation**: Offline analysis capability

---

## 📚 Documentation & Support

### Documentation Files
- **[CHECKLIST_COMPLETION_GUIDE.md](CHECKLIST_COMPLETION_GUIDE.md)** - USDA publication requirements
- **[R_SQUARED_SLOPE_ANALYSIS_GUIDE.md](R_SQUARED_SLOPE_ANALYSIS_GUIDE.md)** - Biological interpretation guide  
- **[FUNCTION_FIX_SUMMARY.md](FUNCTION_FIX_SUMMARY.md)** - Technical fixes and security improvements
- **[QUICK_DEPLOYMENT_GUIDE.md](QUICK_DEPLOYMENT_GUIDE.md)** - Deployment instructions

### Getting Help

**Primary Support:**
- **Contact**: Jason Tidwell, USDA-ARS (jason.tidwell@usda.gov)
- **GitHub Issues**: [Repository Issues Page](https://github.com/USDA-REE-ARS/TickBioassayProbit-web-api/issues)
- **GitHub Discussions**: [Community Discussions](https://github.com/USDA-REE-ARS/TickBioassayProbit-web-api/discussions)

**For Security Vulnerabilities:**
Follow the disclosure policy above - email with "Security Vulnerability" in subject line.

---

## 🤝 Contributing & Development

### Areas for Improvement
- [ ] Additional statistical tests (probit vs logit comparison)
- [ ] More visualization options (3D plots, heat maps)
- [ ] Export format enhancements (Excel, CSV)
- [ ] Batch analysis capabilities
- [ ] Additional arthropod species support
- [ ] Field data integration tools
- [ ] Multi-language support

### Contributing Process
1. Fork the repository
2. Create feature branch (`git checkout -b feature/amazing-feature`)
3. Make changes following code style guidelines
4. Ensure all security scans pass
5. Submit pull request with detailed description

### Development Standards
- **Security**: All contributions must pass Trivy and Dependabot scans
- **Testing**: Include test data and validation procedures
- **Documentation**: Update README and guides as needed
- **Statistical validity**: Maintain rigorous statistical methodology

---

## 📄 License & Legal

### U.S. Government Work Notice

This software was developed by an employee of the United States Department of Agriculture, Agricultural Research Service (USDA-ARS), as part of official duties.

Pursuant to 17 U.S.C. § 105, this work is not subject to copyright protection in the United States and is therefore in the public domain within the United States.

### License for Reuse

To facilitate international use and provide a standard legal framework, this software is also distributed under the MIT License. See the LICENSE file for details.

### Disclaimer

The software is provided "as is", without warranty of any kind, express or implied, including but not limited to the warranties of merchantability, fitness for a particular purpose, and noninfringement.

The use of this software does not constitute an endorsement by USDA-ARS of any commercial product or service.

See [License.txt](LICENSE) file for complete legal details.

---

## 📊 Citation & Attribution

### For Publications

**In Methods Section:**
"Probit regression analysis was performed using the USDA-ARS Probit Analysis Tool v10.3 (Tidwell, 2024) accessed at [URL]."

**In References:**
```bibtex
@software{tidwell2024probit,
  title = {Probit Analysis Tool for Acaricide Resistance Research},
  author = {Jason Tidwell},
  institution = {USDA Agricultural Research Service},
  year = {2024},
  url = {https://github.com/USDA-REE-ARS/TickBioassayProbit-web-api},
  version = {10.3},
  note = {Web-based bioassay analysis tool}
}
```

### Attribution Request
While not required, please consider citing this tool in publications to help track its scientific impact and support continued development.

---

## 🙏 Acknowledgments

### Development Team
**Lead Developer:** Jason Tidwell, Microbiologist  
**Institution:** USDA Agricultural Research Service  
**Facility:** Adkisson-Pfrimmer Agricultural Genomics Center  
**Location:** College Station, TX

### Special Thanks
- **USDA-REE** for supporting open science initiatives
- **Global acaricide resistance research community** for feedback and testing
- **Streamlit team** for the excellent web framework
- **Statsmodels developers** for robust statistical implementations
- **Open source scientific Python community** for foundational libraries

### Research Mission
Developed for researchers who need accessible, reliable bioassay analysis tools without installation barriers. This tool represents USDA-ARS's commitment to providing public domain software that advances agricultural research and global food security.

---

## 🔗 Related Resources

### Scientific Resources
- **WHO Guidelines**: Pesticide resistance testing protocols
- **IRAC Guidelines**: Insecticide resistance management
- **Robertson & Preisler (1992)**: "Pesticide Bioassays with Arthropods" (reference methods)

### Alternative Tools
- **PoloPlus**: Commercial probit analysis software
- **R Package MASS**: `dose.p()` function for R users  
- **SAS PROC PROBIT**: Enterprise statistical software option
- **Desktop Version**: Full-featured Python package (if developed)

### Technical Resources
- **Streamlit Documentation**: [docs.streamlit.io](https://docs.streamlit.io)
- **Statsmodels GLM Guide**: Statistical implementation details
- **Python Scientific Stack**: NumPy, SciPy, Pandas documentation

---

## 📈 Version History & Roadmap

### Current Version: v10.3 - Enhanced Biological Analysis
**Released:** December 2024

**New Features:**
- ✅ R² calculation and biological interpretation
- ✅ Slope analysis with mechanistic insights  
- ✅ Enhanced PDF reports with parameter tables
- ✅ Improved biological interpretation sections
- ✅ Fixed critical bugs and security vulnerabilities
- ✅ Added comprehensive error handling

### Version History
- **v10.2**: Initial web version release
- **v10.1**: Single dataset analysis with basic features
- **v10.0**: Desktop version (proof of concept)

### Planned Features
- **v10.4**: Advanced comparison analytics and batch processing
- **v10.5**: Field data integration and GPS mapping
- **v10.6**: Multi-species support and protocol templates  
- **v11.0**: Machine learning resistance prediction models

### Feedback Integration
Version development prioritizes user feedback from:
- Research community testing
- GitHub issue reports
- Direct researcher contact
- Scientific conference demonstrations

---

## 🎉 Get Started Today!

### Ready to Analyze Your Data?

**[🚀 Launch Web App](https://your-app.streamlit.app)**

**No installation • No login • No cost • Just science!**

### Try With Example Data
1. Click the link above
2. Use the provided example datasets
3. Run analysis in under 30 seconds
4. Download your first PDF report

### Need Help Getting Started?
- Review the data format section above
- Check out example files in the repository
- Contact support for assistance
- Join the GitHub discussions

---

**Made with ❤️ for the global acaricide resistance research community**  
**USDA Agricultural Research Service | Public Domain Software | Version 10.3**

---
*Last updated: December 2024 | Next review: June 2025*
