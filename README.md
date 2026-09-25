# 📊 TickBioassayProbit-web-api

**Web-based probit analysis tool for acaricide and bioassay research**

Professional grouped-binomial probit regression analysis delivered through a Streamlit web application.

[![Streamlit App](https://static.streamlit.io/badges/streamlit_badge_black_white.svg)](https://tickbioprobit.streamlit.app)

---

## 🎯 Mission & Purpose

### USDA-ARS Mission Alignment

This tool supports the USDA Agricultural Research Service mission by:

- **Advancing agricultural research** through accessible statistical tools
- **Supporting pest and acaricide-resistance monitoring**
- **Enabling collaboration** among research institutions
- **Reducing software-installation barriers** through a web interface
- **Promoting reproducible and transparent analysis**

**Developed by:** Jason Tidwell, Microbiologist  
**Institution:** USDA Agricultural Research Service, Cattle Fever Tick Research Unit  
**Location:** Edinburg, Texas, USA  
**Current version:** **11.8 Web**

### Why the Tool Exists

Bioassay researchers often need probit analysis without installing specialized commercial or desktop software. TickBioassayProbit provides a browser-accessible workflow for fitting concentration-response models, estimating lethal concentrations, comparing populations, and generating PDF reports.

The application runs on the server hosting the Streamlit deployment. Users interact with it through a web browser.

### Target Users

- Research entomologists studying acaricide or insecticide resistance
- Toxicologists conducting concentration-response experiments
- USDA, university, and extension researchers
- International collaborators
- Researchers generating resistance phenotypes for genetic or QTL studies
- Other investigators using grouped binary-response bioassays

---

## ✨ Version 11.8 Highlights

### Single-Dataset Analysis

- Grouped-binomial GLM with a **probit link**
- User-selectable lethal concentration levels, with defaults of **LC1, LC50, and LC99**
- Delta-method **95% confidence intervals**
- Interpolation/extrapolation assessment based on the tested concentration range
- Pearson chi-square goodness-of-fit
- Pearson dispersion statistic
- Descriptive probit-scale R²
- Dose-response slope and intercept
- Replicate variability summaries
- Abbott correction when untreated-control mortality is present
- Mortality and probit regression plots
- Downloadable PDF report

### Multi-Dataset Analysis

- Upload and analyze **multiple populations/strains at once**
- Datasets organized alphabetically by population/strain name
- Reference-vs-all comparisons
- All-pairwise comparisons
- Global dataset × log10(concentration) interaction test
- Global dataset-shift test under a common-slope model
- Holm, Bonferroni, or no multiple-testing adjustment
- LC50 forest plots
- Multi-dataset mortality and probit plots
- Consolidated multi-dataset PDF reports

### Ratio Reporting

The application distinguishes between two different comparison concepts:

**Resistance Ratio (RR)**  
Used when the user explicitly designates a dataset as a susceptible reference:

```text
RR = LC50(test population) / LC50(susceptible reference)
```

The susceptible reference remains the denominator, so the direction of the ratio is preserved.

**LC50 Fold-Difference**  
Used for all-pairwise comparisons where no susceptible-reference direction is assumed:

```text
LC50 fold-difference = larger LC50 / smaller LC50
```

The fold-difference is therefore always ≥ 1. The application separately identifies the dataset with the higher and lower LC50 values.

Ratio and fold-difference confidence intervals are individual 95% confidence intervals calculated with the delta method on the log10 ratio scale.

### Flexible Data Import

Version 11 accepts:

- Legacy metadata-first text files
- Header-first CSV files
- Header-first TSV files
- Other simple delimited text tables
- Common alternative column names
- Manual column mapping when automatic detection is insufficient

Strain/population and chemical metadata are automatically populated when present in the uploaded file and remain editable.

### Shared Concentration Units

- Units embedded in uploaded data are detected automatically.
- If units are missing, the user enters the concentration units **once for the analysis**.
- The shared value is applied to datasets that lack embedded unit metadata.
- Conflicting embedded units are flagged.
- Multi-dataset analysis requires user confirmation that selected datasets use the same concentration scale.

### Optional AI Results Assistant

Version 11 includes an optional AI interpretation interface.

The statistical engine remains authoritative. When enabled, the assistant receives **derived analysis summaries and metadata rather than the uploaded raw bioassay table** and is instructed not to recalculate or replace statistical results.

The AI feature is disabled unless an administrator configures an approved endpoint using server-side environment variables:

```text
PROBIT_AI_ENDPOINT
PROBIT_AI_MODEL
PROBIT_AI_API_KEY   # optional, depending on the approved endpoint
```

No AI credential is embedded in the repository.

---

## 🚀 Quick Start

### Option 1: Use the Online Application

**Launch:** [https://tickbioprobit.streamlit.app](https://tickbioprobit.streamlit.app)

No local installation is required for end users.

### Option 2: Run Locally

```bash
# Clone repository
git clone https://github.com/USDA-REE-ARS/TickBioassayProbit-web-api.git
cd TickBioassayProbit-web-api

# Install dependencies
pip install -r requirements.txt

# Run the current application
streamlit run probit_web_app_v11.8.py
```

Streamlit normally opens the local application at:

```text
http://localhost:8501
```

If the repository uses a different deployment entry-point filename, run that file instead.

---

## 📋 Input Data

### Required Analytical Fields

Each dataset requires:

| Field | Description |
|---|---|
| `concentration` | Tested concentration or dose |
| `n` | Number of individuals tested |
| `mortality` | Number of individuals dead |

The program recognizes common alternatives such as `conc`, `dose`, `total`, `n tested`, `dead`, and `deaths`. Unusual headings can be assigned through the Column Mapping panel.

### Legacy Format

```text
Yucatan
Amitraz
concentration    n    mortality
0.250            100  95
0.125            80   60
0.063            75   20
0.031            80   5
```

Line 1 contains the strain/population and line 2 contains the chemical.

### Conventional Header-First Format

```csv
population,acaricide,concentration,units,total,deaths
Yucatan,Amitraz,0.250,ppm,100,95
Yucatan,Amitraz,0.125,ppm,80,60
Yucatan,Amitraz,0.063,ppm,75,20
Yucatan,Amitraz,0.031,ppm,80,5
```

Constant metadata columns may be used to auto-populate population, chemical, and units.

### Untreated Controls

A row with:

```text
concentration = 0
```

is interpreted as an untreated control.

When pooled control mortality is greater than 0%, treatment mortality is adjusted using Abbott's correction before model fitting. Control observations themselves are excluded from the concentration-response regression.

### Validation Rules

The application checks for:

- Numeric concentration, `n`, and mortality values
- Positive sample sizes
- Whole-number raw counts for `n` and mortality
- Mortality ≤ number tested
- At least **3 distinct positive treatment concentrations**
- Mixed strain/population metadata within a single uploaded dataset
- Mixed chemical metadata within a single uploaded dataset
- Incompatible concentration units
- Replicate variability
- Unusual or non-positive fitted slopes

Four or more treatment concentrations are recommended when practical, and concentration selection should provide useful coverage of the response range.

---

## 🎓 How to Use

### 1. Upload Data

1. Open **Upload Data**.
2. Upload one or more `.txt`, `.tsv`, or `.csv` files.
3. Review the auto-detected strain/population and chemical.
4. Confirm or correct column assignments.
5. Review validation messages.
6. If units are not embedded in the files, enter the concentration units once in the shared units field.
7. For multiple datasets, confirm that the concentration units are directly comparable.

### 2. Individual Analysis

1. Open **Individual Analyses**.
2. Select the desired population.
3. Run the analysis.
4. Review:
   - LC estimates and confidence intervals
   - interpolation/extrapolation status
   - slope and intercept
   - Pearson chi-square
   - Pearson p-value
   - dispersion
   - replicate variability
   - fitted plots
5. Download the individual PDF report.

### 3. Multi-Dataset Comparison

Select two or more validated datasets and choose:

#### Reference vs All

Use when every test population should be compared with a selected reference.

If the reference is explicitly identified as a **susceptible reference population**, comparison values are reported as resistance ratios.

#### All Pairwise

Every selected population is compared with every other selected population.

Pairwise magnitude is reported as the **LC50 fold-difference**, with the lower LC50 placed in the denominator.

Example:

```text
Population A LC50 = 0.01
Population B LC50 = 0.50

LC50 fold-difference = 0.50 / 0.01 = 50-fold
```

The output identifies Population B as the higher-LC50 dataset.

### 4. Multiple-Comparison Adjustment

Available options:

- **Holm** — recommended default
- **Bonferroni**
- **None**

Adjusted p-values are provided for pairwise curve-comparison tests. LC50 ratio/fold-difference confidence intervals remain individual 95% confidence intervals and are not converted into simultaneous confidence intervals by the p-value adjustment.

---

## 📊 Statistical Methods

### Probit Regression

Each dataset is modeled using a grouped-binomial generalized linear model:

```text
Probit(mortality probability) = β0 + β1 × log10(concentration)
```

- **Family:** Binomial
- **Link:** Probit
- **Estimation:** Maximum likelihood through Statsmodels GLM

Raw 0% and 100% treatment responses are retained in model fitting. Small boundary adjustments are used only where needed to display empirical values on the probit scale.

### Lethal Concentration Estimation

For a requested mortality level x:

```text
log10(LCx) = [Probit(x) - β0] / β1
```

LC confidence intervals are calculated using the delta method and the fitted parameter covariance matrix.

LC estimates are suppressed when the fitted concentration-response slope is non-positive.

### LC Coverage

Each LC estimate is classified as:

- **Interpolated** — estimated concentration falls within the tested positive-concentration range
- **Extrapolated** — estimated concentration lies outside the tested positive-concentration range

The program also reports whether the requested target response lies within the observed treatment-mortality range.

### Pearson Goodness-of-Fit

The application reports the Pearson chi-square statistic from the grouped-binomial GLM together with residual degrees of freedom and a p-value.

- **p < 0.05:** evidence of lack of fit and/or extra-binomial variation
- **p ≥ 0.05:** no statistically significant lack of fit is detected

A non-significant result does not prove that the probit model is correct.

### Pearson Dispersion

```text
φ = Pearson χ² / residual df
```

Values above 1 indicate residual variability greater than the nominal binomial expectation.

### Descriptive Probit-Scale R²

The application also provides an auxiliary probit-scale R² describing how closely empirical probits track the fitted line.

This quantity is **descriptive only** and is not used as the formal GLM goodness-of-fit test.

### Global Multi-Dataset Comparison

For multiple datasets, three nested grouped-binomial probit models are fitted:

```text
Full model:
log10(concentration) + dataset + dataset × log10(concentration)

Parallel model:
log10(concentration) + dataset

Common model:
log10(concentration)
```

The application uses likelihood-ratio tests to evaluate:

1. **Global slope interaction**  
   Full model vs parallel model

2. **Global dataset shift under a common slope**  
   Parallel model vs common model

A significant slope interaction indicates heterogeneous concentration-response slopes and means that a single LC50 ratio or fold-difference does not fully describe the difference among curves.

### Pairwise Curve Tests

For each pair, combined grouped-binomial probit models are used to test:

- slope interaction/non-parallelism
- dataset shift under a common-slope model

Raw and multiple-testing-adjusted p-values are reported.

### Resistance Ratio

When a susceptible reference is explicitly designated:

```text
RR = LC50(test) / LC50(susceptible reference)
```

The 95% CI is calculated by delta-method propagation of the two independent log10(LC50) variances.

### LC50 Fold-Difference

For all-pairwise comparisons:

```text
Fold-difference = larger LC50 / smaller LC50
```

When the directional ratio is inverted, its confidence interval is also transformed:

```text
1 / ratio
CI = (1 / original upper, 1 / original lower)
```

The fold-difference point estimate is always ≥ 1, while the confidence interval can cross 1 when the direction of the true difference is uncertain.

---

## 📈 Output and Visualization

### Individual Analysis

- LC table
- 95% confidence intervals
- LC coverage flags
- slope and intercept
- Pearson chi-square and p-value
- dispersion
- descriptive probit-scale R²
- replicate variability
- mortality-response curve
- probit regression plot
- PDF report

### Multi-Dataset Analysis

- Alphabetically ordered dataset summary
- Global model-comparison tests
- Reference-based or all-pairwise comparison table
- Raw and adjusted p-values
- LC50 forest plot
- Multi-dataset mortality curves
- Multi-dataset probit regression plot
- Ratio/fold-difference forest plot
- Consolidated PDF report

Very small p-values are displayed in scientific notation for readability. Ordinary slope values are displayed to three decimal places, with scientific notation used for very small or very large values.

---

## 📄 PDF Reports

Both individual and multi-dataset analyses can generate downloadable PDF reports.

Reports may contain:

- assay metadata
- raw data tables
- control mortality and Abbott-correction information
- LC estimates and confidence intervals
- interpolation/extrapolation information
- model parameters
- goodness-of-fit statistics
- dispersion
- comparison statistics
- resistance ratios or LC50 fold-differences
- fitted plots

The PDF export path writes and verifies a completed PDF before making it available through the Streamlit download button.

---

## 🔒 Security & Vulnerability Disclosure

### Vulnerability Disclosure Policy

The community is encouraged to responsibly disclose vulnerabilities so they can be addressed without exposing users unnecessarily.

If you discover a security vulnerability:

1. **Email:** jason.tidwell@usda.gov  
   Subject: `Security Vulnerability - Probit Tool`
2. **Provide details:** description, steps to reproduce, and potential impact.
3. **Confidential handling:** reports will be evaluated before public disclosure.
4. **Recognition:** contributors may be acknowledged, with permission, after resolution.

**Please do not publicly disclose vulnerabilities until they have been addressed.**

### Vulnerability Response Timeline

When vulnerabilities are identified:

- **Critical vulnerabilities:** patched within 7 days or the application may be taken offline
- **High vulnerabilities:** addressed within 14 days
- **Medium/Low vulnerabilities:** resolved within 30 days
- **Users notified:** through GitHub releases and repository notices when appropriate
- **Workarounds provided:** when an immediate fix is not possible

If a vulnerability cannot be resolved in a timely manner, a prominent warning may be added to the README and the application may be temporarily taken offline until the issue is addressed.

### Repository Security Practices

The project uses repository security and dependency-management practices that may include:

- Dependabot
- Trivy scanning
- code review
- dependency updates
- input validation
- file-type restrictions
- HTTPS when provided by the hosting environment
- server-side handling of optional AI credentials

Repository dependency versions may be maintained at newer minimum releases to satisfy security or organizational requirements.

---

## 🔐 Data Handling and Privacy

Uploaded files are processed by the Python process hosting the Streamlit application.

The application itself does not intentionally persist uploaded assay files to permanent storage; however:

- transport security
- server logging
- temporary storage
- retention
- access controls

depend on the environment in which the application is deployed.

Users and administrators should follow applicable USDA or organizational requirements before submitting sensitive information.

### AI Data Handling

When the optional AI assistant is enabled:

- the application sends structured derived statistics and assay metadata to the configured endpoint;
- the uploaded raw observation table is not included in the AI context by the application;
- endpoint authorization, retention, security, and approval remain the responsibility of the deployment administrator.

---

## ⚠️ Common Issues

### Model Does Not Converge

Possible causes include:

- concentrations cover too narrow a range
- little change in mortality among concentrations
- extreme separation
- high replicate variability
- insufficient concentration levels

### Non-Positive Slope

The application will fit and display the diagnostic model, but LC estimates and LC50 ratio/fold-difference calculations are suppressed.

Review:

- concentration coding
- assay direction
- data transcription
- treatment range
- biological response

### Extrapolated LC Estimate

An LC estimate is flagged as extrapolated when the estimated concentration lies outside the tested positive-concentration range.

Consider testing additional concentrations closer to the target response.

### Significant Pearson Test

A significant Pearson goodness-of-fit result indicates evidence of lack of fit and/or extra-binomial variation.

Review:

- replicate variability
- concentration spacing
- experimental consistency
- possible heterogeneity
- additional model diagnostics

Do not assume that a significant lack-of-fit test leaves all estimates unaffected.

### Multi-Dataset Comparison Disabled

Check for:

- different chemicals
- incompatible concentration units
- missing unit confirmation
- validation errors
- fewer than two successfully fitted datasets

---

## 💻 Technical Specifications

### Technology Stack

- **Web framework:** Streamlit
- **Language:** Python
- **Scientific computing:** NumPy, pandas, SciPy
- **Statistical modeling:** Statsmodels
- **Visualization:** Matplotlib
- **PDF generation:** FPDF2-compatible `FPDF` interface
- **Optional AI connection:** Python standard-library HTTP client to an administrator-configured compatible endpoint

### Dependencies

See:

```text
requirements.txt
```

Repository dependency constraints may be stricter than the minimum required by the source code because of security-scanning and organizational requirements.

### Deployment Options

- Streamlit Community Cloud
- USDA/institutional server
- containerized deployment
- local Streamlit installation

---

## 📚 Statistical References

1. **Finney, D.J. (1971).** *Probit Analysis*, 3rd ed. Cambridge University Press.
2. **Robertson, J.L., Russell, R.M., Preisler, H.K., and Savin, N.E. (2007).** *Bioassays with Arthropods*, 2nd ed. CRC Press.
3. **Ver Hoef, J.M. (2012).** Who invented the delta method? *The American Statistician* 66:124-127.
4. **Abbott, W.S. (1925).** A method of computing the effectiveness of an insecticide. *Journal of Economic Entomology* 18:265-267.

---

## 📚 Documentation & Support

### Repository

[USDA-REE-ARS/TickBioassayProbit-web-api](https://github.com/USDA-REE-ARS/TickBioassayProbit-web-api)

### Support

**Contact:** Jason Tidwell, USDA-ARS  
**Email:** jason.tidwell@usda.gov

For bugs and feature requests, use the GitHub Issues page associated with the repository.

---

## 🤝 Contributing

Contributions should:

1. use a feature branch;
2. preserve statistical validity;
3. include appropriate validation or test data;
4. maintain clear documentation;
5. pass applicable repository security/dependency checks;
6. avoid embedding credentials or sensitive deployment information.

---

## 📄 License & Legal

### U.S. Government Work Notice

This software was developed by an employee of the United States Department of Agriculture, Agricultural Research Service (USDA-ARS), as part of official duties.

Pursuant to 17 U.S.C. § 105, U.S. Government works are not subject to copyright protection in the United States.

### License for Reuse

To facilitate international use and provide a standard legal framework, this software is also distributed under the **MIT License**. See the repository license file for the controlling terms.

### Disclaimer

The software is provided "as is", without warranty of any kind, express or implied, including but not limited to warranties of merchantability, fitness for a particular purpose, and noninfringement.

Use of this software does not constitute USDA-ARS endorsement of any commercial product or service.

---

## 📊 Citation & Attribution

### Suggested Methods Wording

> Probit regression analysis was performed using the USDA-ARS Probit Analysis Tool v11.8. The application fits grouped-binomial generalized linear models with a probit link to log10-transformed concentration data. Lethal concentration estimates and individual 95% confidence intervals were calculated using the delta method. Multi-population comparisons used nested grouped-binomial probit models and likelihood-ratio tests, with multiple-comparison adjustment where specified.

### Suggested Software Citation

```bibtex
@software{tidwell2026probit,
  title       = {Probit Analysis Tool for Bioassay and Acaricide Resistance Research},
  author      = {Tidwell, Jason},
  institution = {USDA Agricultural Research Service},
  year        = {2026},
  url         = {https://github.com/USDA-REE-ARS/TickBioassayProbit-web-api},
  version     = {11.8},
  note        = {Web-based grouped-binomial probit analysis tool}
}
```

---

## 📈 Version History

### v11.8 — Current

Key capabilities include:

- multi-file dataset upload
- individual and multi-dataset analysis
- reference-vs-all and all-pairwise comparisons
- LC50 fold-differences for pairwise comparisons
- directional resistance ratios for designated susceptible references
- Holm/Bonferroni multiplicity adjustment
- global and pairwise combined-model tests
- shared concentration-unit workflow
- alphabetically ordered datasets
- Pearson dispersion
- extrapolation/interpolation flags
- scientific notation for very small p-values
- verified individual and multi-dataset PDF generation
- optional administrator-configured AI results assistant

### Recent Development Milestones

- **v11.7:** shared concentration-units workflow
- **v11.6:** improved scientific-notation formatting for small p-values
- **v11.5:** standardized slope display formatting
- **v11.4:** alphabetical dataset ordering
- **v11.3:** corrected PDF byte-generation path
- **v11.2:** LC50 fold-difference reporting for all-pairwise comparisons
- **v11.0:** multi-dataset architecture and optional AI assistant
- **v10.7:** major statistical and validation revision
- **v10.x:** earlier single/two-dataset web application development

---

## 🎉 Get Started

**[🚀 Launch the Web App](https://tickbioprobit.streamlit.app)**

No local installation is required for users of the hosted application.

For development or offline use, clone the repository and run the Streamlit entry point locally.

---

**USDA Agricultural Research Service | Cattle Fever Tick Research Unit | Version 11.8**

*README last updated: September 2026*
