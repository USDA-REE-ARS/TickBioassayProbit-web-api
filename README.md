# 📊 TickBioassayProbit-web-api

**Web-based probit analysis tool for acaricide resistance research**

Professional bioassay probit regression analysis - **run directly in your browser!**

[![Streamlit App](https://static.streamlit.io/badges/streamlit_badge_black_white.svg)](https://tickbioprobit.streamlit.app)

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
**Institution:** USDA ARS Cattle Fever Tick Research Unit  
**Location:** Edinburg, TX  
**Version:** 11.9 (Web Application)

### Key Problem Solved
Many research institutions have IT restrictions preventing software installation. This web-based tool eliminates installation barriers by providing a browser-based interface to the hosted Streamlit application, making probit analysis accessible to researchers without requiring local software installation or programming knowledge.

### Target Users
- Research entomologists studying acaricide resistance
- Toxicologists conducting dose-response experiments
- University researchers and extension specialists
- International collaborators at institutions with restrictive IT policies
- QTL researchers needing standardized phenotyping

---

## ✨ Core Features

- ✅ **No Installation Required** - Works in any browser
- ✅ **Flexible Data Upload** - Supports legacy text files plus CSV/TSV formats with column mapping
- ✅ **Professional Results** - Publication-quality analysis and downloadable PDF reports
- ✅ **Single & Multi-Dataset Analysis** - Analyze one population or compare multiple populations
- ✅ **Statistical Diagnostics** - Pearson chi-square, dispersion, descriptive R², slope, and LC confidence intervals
- ✅ **Resistance Comparisons** - Resistance Ratios for designated susceptible references and LC50 Fold-Differences for pairwise comparisons
- ✅ **User-Friendly** - Drag-and-drop interface with shared concentration-unit entry and validation
- ✅ **Optional AI Assistant** - Available only when an approved administrator-configured endpoint is enabled
- ✅ **Free & Open** - Public domain software in the United States with MIT licensing for reuse

---

## 🚀 Quick Start

### Option 1: Use Online (Recommended)

**Just click:** [https://tickbioprobit.streamlit.app](https://tickbioprobit.streamlit.app)

No setup needed - start analyzing immediately!

### Option 2: Run Locally

```bash
# Clone repository
git clone https://github.com/USDA-REE-ARS/TickBioassayProbit-web-api.git
cd TickBioassayProbit-web-api

# Install dependencies
pip install -r requirements.txt

# Run app
streamlit run probit_web_app.py
```

#Opens at `http://localhost:8501`

---

## 📊 What It Analyzes

### Probit Regression Analysis for Bioassay Data

**Core Analysis:**
- **LC Estimates**: User-selected lethal concentration estimates (default LC1, LC50, LC99) with 95% confidence intervals
- **Resistance Ratios**: Test population LC50 divided by a user-designated susceptible reference LC50
- **LC50 Fold-Differences**: Pairwise comparison using the larger LC50 divided by the smaller LC50
- **Model Diagnostics**: Pearson chi-square goodness-of-fit, residual degrees of freedom, dispersion, slope, intercept, and descriptive probit-scale R²
- **Data Quality**: Replicate variability, concentration-range coverage, control mortality, and extrapolation warnings

**Multi-Dataset Analysis:**
- **Multiple Populations**: Upload and analyze multiple datasets in one session
- **Reference vs All**: Compare selected populations with one reference dataset
- **All Pairwise**: Compare every selected dataset pair
- **Global Curve Tests**: Dataset × concentration interaction and common-slope dataset-shift likelihood-ratio tests
- **Multiple Testing**: Holm (default), Bonferroni, or no p-value adjustment

**Visualizations:**
- Mortality concentration-response curves
- Probit regression plots with fitted lines
- Multi-dataset comparison plots
- LC50 forest plots
- Ratio/fold-difference forest plots

**Report Generation:**
- Individual and multi-dataset PDF reports with embedded plots
- Statistical parameter tables
- LC estimates and confidence intervals
- Comparison statistics and ratio/fold-difference results

### Perfect For:
- **Acaricide resistance testing** (primary use case)
- **Insecticide resistance monitoring**
- **QTL mapping phenotyping**
- **Toxicology concentration-response studies**
- **Other grouped binary-outcome bioassays**

---

## 📋 Data Format & Requirements

### Input File Format

The application accepts legacy metadata-first text files and conventional header-first CSV/TSV/text files.

**Legacy format:**

```text
Strain_Name
Chemical_Name
concentration    n    mortality
0.500            96   96
0.350            102  79
0.245            161  114
0.125            98   45
0.063            105  18
0.031            102  5
```

**Header-first format:**

```csv
population,chemical,units,concentration,n,mortality
Strain_Name,Chemical_Name,ppm,0.500,96,96
Strain_Name,Chemical_Name,ppm,0.350,102,79
Strain_Name,Chemical_Name,ppm,0.245,161,114
```

**Required Analytical Fields:**
- `concentration`: Concentration tested (numeric)
- `n`: Number of individuals tested (whole-number count)
- `mortality`: Number that died (whole-number count ≤ n)

Common alternative column names are detected automatically, and columns can be mapped manually in the application.

**Optional Metadata:**
- Strain/population
- Chemical/acaricide
- Concentration units

If concentration units are not included in the uploaded files, the user enters them once for the analysis. For multi-dataset comparison, the user confirms that all selected datasets use the same directly comparable concentration units.

**Data Requirements:**
- At least **3 distinct positive treatment concentrations**
- Mortality must be ≤ n for each row
- `n` and mortality must contain whole-number raw counts
- Replicates are recommended when practical
- Concentrations should span enough of the response range to support the LC estimates of interest

A row with `concentration = 0` is treated as an optional untreated control. If control mortality is greater than 0%, Abbott's correction is applied to treatment mortality before model fitting.

The **Help** tab in Version 11.9 includes small working examples of both supported dataset formats.

[Download example files from repository](examples/)

---

## 🎓 How to Use

### Single Dataset Analysis

**Step 1: Upload Data**
1. Go to "Upload Data" tab
2. Click "Browse files" or drag-and-drop your .txt, .tsv, or .csv file
3. Review the detected strain/population, chemical, and column mapping
4. Check validation results and data preview
5. Enter shared concentration units if they are not available in the uploaded data

**Step 2: Run Analysis**
1. Navigate to "Individual Analyses" tab
2. Select the dataset
3. Click "Run Individual Analysis"

**Step 3: Interpret Results**
- **LC Estimates**: Review the requested LC values and 95% confidence intervals
- **LC Coverage**: Note whether an estimate is interpolated or extrapolated
- **Model Parameters**: Review slope and intercept
- **Model Fit**: Review Pearson chi-square, residual df, p-value, and dispersion
- **Descriptive R²**: Use as an auxiliary description of the fitted probit line, not as the formal goodness-of-fit test
- **Replicate Variability**: Review concentration-specific variation

**Step 4: Save Results**
- Download the individual PDF report
- Copy results for manuscripts
- Save figures as needed

### Multi-Dataset Comparison

**Step 1: Upload Datasets**
- Upload two or more datasets
- Confirm that each dataset passes validation
- Confirm that the chemical and concentration units are compatible

**Step 2: Select Comparison Mode**
1. Go to "Multi-Dataset Comparison"
2. Choose **Reference vs all** or **All pairwise**
3. Select Holm, Bonferroni, or no multiple-comparison adjustment
4. If using a susceptible reference, identify it explicitly before interpreting the ratio as a Resistance Ratio

**Step 3: Interpret Comparison**
- **Resistance Ratio**: Used when the denominator is explicitly designated as a susceptible reference
- **LC50 Fold-Difference**: Used for all-pairwise comparisons; the smaller LC50 is placed in the denominator
- **Slope Interaction Test**: Tests for non-parallel concentration-response slopes
- **Dataset Shift Test**: Tests for a population shift under a common-slope model
- **Global Tests**: Summarize slope heterogeneity and common-slope dataset shifts across all selected datasets

---

## 📖 Example Analysis Results

### Pera F3 Strain vs Susceptible Control

```text
Dataset: Pera F3 strain tested with Coumaphos
Reference: Susceptible Deutsch strain

Resistance Analysis:
  Test LC50:       1.523 (95% CI: 1.445 - 1.607)
  Reference LC50:  0.010 (95% CI: 0.009 - 0.011)

  Resistance Ratio: 152.3x (95% CI: 138.2 - 168.1)

Model Summary:
  Descriptive probit-scale R² = 0.889
  Slope = 3.45

Statistical Tests:
  Pearson goodness-of-fit: χ² = 12.45, df = 19, p = 0.789
  Slope interaction: p = 0.234

Interpretation:
  The test population has a substantially higher LC50 than the designated
  susceptible reference in this illustrative example. The slope-interaction
  test does not detect a significant difference in slopes. Slope and R²
  describe features of the fitted response but do not by themselves identify
  a molecular resistance mechanism.
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
- Uploaded assay files are processed by the Streamlit server hosting the application
- The application does not intentionally persist uploaded raw assay files to permanent storage
- No user account is required by the application itself
- Session handling, temporary storage, logs, retention, and transport security depend on the deployment environment
- Users should follow applicable organizational requirements before submitting sensitive data

### Optional AI Assistant
- The AI assistant is disabled unless an administrator configures an approved endpoint
- When enabled, the application sends structured analysis summaries and assay metadata rather than the uploaded raw observation table
- AI credentials are supplied through server-side environment variables and are not embedded in the source code
- Endpoint authorization, retention, and data-handling requirements remain the responsibility of the deployment administrator

### Security Implementation
- **HTTPS encryption** when provided by the hosting environment
- **Input validation** and file-type checking
- **Safe error handling**
- **Regular dependency updates** via Dependabot automation
- **Static code analysis** via Trivy security scanning
- **Server-side credential handling** for optional AI configuration

### Compliance
- No PII collection is required by the application workflow
- Public domain U.S. Government work with MIT licensing for reuse
- Deployment administrators are responsible for applicable organizational security and data-handling requirements

---

## ⚠️ Common Issues & Solutions

### Data Upload Issues

**"File format not recognized"**
- Confirm that the file is a supported .txt, .tsv, or .csv file
- Verify that the table contains concentration, n tested, and mortality/deaths fields
- Use the Column Mapping controls if the headings are unusual
- See the working examples in the Help tab

**"Mortality exceeds sample size"**
- Check data: mortality must be ≤ n for every row
- Look for data entry errors
- Verify numbers against laboratory records

**"Not enough treatment concentrations"**
- At least 3 distinct positive concentrations are required
- Additional concentrations are recommended when practical to improve response-range coverage

### Analysis Problems

**"Evidence of lack of fit" (Pearson p < 0.05)**
- Review replicate variability
- Check concentration spacing and assay consistency
- Review whether the concentration range adequately captures the response
- Interpret LC estimates and comparisons cautiously when lack of fit is substantial

**"High variability (CV% > 20%)" warnings**
- Review experimental protocol consistency
- Check whether specific concentrations are unusually variable
- Consider whether additional replication is warranted
- Document important variability in methods/results

**"Non-positive slope"**
- Check concentration coding and data transcription
- Review whether mortality generally increases with concentration
- LC estimates and LC50 ratios/fold-differences are suppressed when the fitted slope is non-positive

**"Extrapolated LC estimate"**
- The estimated LC lies outside the tested positive-concentration range
- Consider additional concentrations closer to the target response

### Application Issues

**App won't load**
- Check internet connection
- Try another modern browser
- Refresh the application
- Verify the Streamlit deployment is online

**PDF report won't download correctly**
- Retry after the analysis has completed
- Confirm the current deployed application version
- Report persistent problems through the repository issue tracker

---

## 🔬 Statistical Methods & Validation

### Probit Regression Implementation
- **Link function**: Probit (inverse normal CDF)
- **Family**: Grouped binomial
- **Predictor**: log10(concentration)
- **Estimation**: Maximum likelihood via Statsmodels GLM
- **LC confidence intervals**: Delta method on the log10 concentration scale
- **Resistance ratio / fold-difference confidence intervals**: Delta method on the log10 ratio scale

### Lethal Concentration Estimates
- User-selected LC levels are supported (default LC1, LC50, LC99)
- LC estimates are suppressed when the fitted slope is non-positive
- Estimates outside the tested concentration range are labeled **Extrapolated**
- The application separately reports whether the target response was represented in the observed mortality range

### Model Diagnostics
- **Goodness-of-fit**: Pearson chi-square test
- **Dispersion**: Pearson χ² / residual df
- **Descriptive R²**: Auxiliary probit-scale summary; not the formal GLM goodness-of-fit test
- **Replicate variability**: Coefficient of variation by concentration
- A non-significant Pearson test does not prove that the model is correct

### Multi-Dataset Comparisons
- **Global slope interaction**: Likelihood-ratio comparison of full and parallel grouped-binomial probit models
- **Global dataset shift**: Likelihood-ratio comparison of parallel and common models
- **Pairwise slope and shift tests**: Combined-model likelihood-ratio tests
- **Multiple testing**: Holm (default), Bonferroni, or none
- P-value adjustment does not convert the reported ratio confidence intervals into simultaneous confidence intervals

### Resistance Ratios and LC50 Fold-Differences
- **Resistance Ratio** = LC50(test) / LC50(susceptible reference) when a susceptible reference is explicitly designated
- **LC50 Fold-Difference** = larger LC50 / smaller LC50 for all-pairwise comparisons
- For an inverted pairwise ratio, the confidence interval is inverted as `(1 / upper, 1 / lower)`

### Numerical Stability and Controls
- Raw 0% and 100% treatment responses are retained for model fitting
- Small boundary adjustments are used only for empirical probit display calculations
- Untreated controls (`concentration = 0`) are excluded from the dose-response fit
- Abbott's correction is applied when untreated-control mortality is greater than 0%
- Model convergence and finite-parameter checks are performed before LC interpretation

---

## 💻 Technical Specifications

### Technology Stack
- **Frontend / Web Framework**: Streamlit
- **Backend**: Python
- **Statistics**: Statsmodels grouped-binomial GLM
- **Scientific Computing**: NumPy, pandas, SciPy
- **Visualization**: Matplotlib with publication-quality output
- **Reports**: FPDF-compatible PDF generation with embedded plots
- **Optional AI**: Administrator-configured compatible endpoint using server-side environment variables

### System Requirements

**For Users (Browser-based):**
- Modern web browser (Chrome, Firefox, Safari, Edge)
- JavaScript enabled
- Internet connection for the hosted version
- No local installation or admin rights required

**For Local Deployment:**
- Python environment compatible with `requirements.txt`
- Requirements: See requirements.txt

### Performance
Performance depends on dataset size, number of populations, report options, and the hosting environment.

### Deployment Options
1. **Streamlit Cloud**
2. **Institutional server**
3. **Docker/container deployment**
4. **Local installation**

---

## 📚 Documentation & Support

This README provides the primary documentation for installation, data formatting, analysis workflow, interpretation, security, and citation of the application.

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

See [License](License.txt) file for complete legal details.

---

## 📊 Citation & Attribution

### For Publications

**In Methods Section:**
"Probit regression analysis was performed using the USDA-ARS Probit Analysis Tool v11.9 (Tidwell, 2026) accessed at https://tickbioprobit.streamlit.app."

**In References:**
```bibtex
@software{tidwell2026probit,
  title = {Probit Analysis Tool for Acaricide Resistance Research},
  author = {Jason Tidwell},
  institution = {USDA Agricultural Research Service},
  year = {2026},
  url = {https://github.com/USDA-REE-ARS/TickBioassayProbit-web-api},
  version = {11.9},
  note = {Web-based bioassay probit analysis tool}
}
```

### Attribution Request
While not required, please consider citing this tool in publications to help track its scientific impact and support continued development.

---

## 🙏 Acknowledgments

### Development Team
**Lead Developer:** Jason Tidwell, Microbiologist  
**Institution:** USDA Agricultural Research Service  
**Facility:** Cattle Fever Tick Research Unit  
**Location:** Edinburg, TX

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

### Current Version: v11.9 - Multi-Dataset Probit Analysis
**Released:** September 2026

**Current Features:**
- ✅ Flexible legacy and header-first data import with column mapping
- ✅ User-selectable LC estimates with delta-method confidence intervals
- ✅ Untreated-control detection and Abbott correction
- ✅ Multi-dataset reference-vs-all and all-pairwise analysis
- ✅ Directional Resistance Ratios for designated susceptible references
- ✅ LC50 Fold-Differences for pairwise comparisons
- ✅ Global and pairwise combined-model curve tests
- ✅ Holm and Bonferroni multiple-testing adjustment
- ✅ Shared concentration-unit workflow
- ✅ Alphabetical dataset organization
- ✅ Pearson dispersion and extrapolation safeguards
- ✅ Scientific notation for very small p-values and extreme slope displays
- ✅ Individual and multi-dataset PDF reports
- ✅ Optional administrator-configured AI results assistant
- ✅ Working dataset-format examples in the Help tab

### Recent Version History
- **v11.9**: Added working dataset-format examples to the Help tab
- **v11.8**: Corrected LC50 ratio/fold-difference calculation compatibility
- **v11.7**: Added shared concentration-unit entry and confirmation
- **v11.6**: Improved scientific notation for small p-values
- **v11.5**: Standardized slope display formatting
- **v11.4**: Alphabetical dataset ordering
- **v11.3**: Corrected PDF generation/download path
- **v11.2**: Added LC50 Fold-Difference reporting for all-pairwise comparisons
- **v11.0**: Added multi-dataset architecture and optional AI assistant
- **v10.x**: Earlier single- and two-dataset web application development

### Feedback Integration
Version development prioritizes user feedback from:
- Research community testing
- GitHub issue reports
- Direct researcher contact
- Scientific conference demonstrations

---

## 🎉 Get Started Today!

### Ready to Analyze Your Data?

**[🚀 Launch Web App](https://tickbioprobit.streamlit.app)**

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
**USDA Agricultural Research Service | Public Domain Software | Version 11.9**

---
*Last updated: September 2026*
