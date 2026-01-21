# Probit Analysis Tool - Deployment Package
## USDA ARS Software Submission

**Application Name:** Web-Based Probit Analysis Tool for Bioassay Research

**Developer:** Jason Tidwell, Research Entomologist  
**Institution:** USDA ARS Adkisson-Pfrimmer Agricultural Genomics Center  
**Location:** College Station, TX  
**Date:** December 2024  
**Version:** 10.3 (Web)

---

## 1. APPLICATION OVERVIEW

### Purpose
Web-based statistical tool for analyzing dose-response bioassay data, specifically designed for acaricide resistance studies in ticks and other arthropods. Implements probit regression analysis with comprehensive diagnostics.

### Target Users
- Research entomologists
- Toxicologists  
- Acaricide resistance researchers
- University researchers
- Extension specialists
- International collaborators at institutions with restrictive IT policies

### Key Problem Solved
Many research institutions have IT restrictions preventing software installation. This web-based tool eliminates installation barriers by running entirely in a web browser, making probit analysis accessible to researchers worldwide without requiring local software, IT approval, or programming knowledge.

### Core Features
- Single dataset probit analysis
- Two-dataset comparison with resistance ratio calculation
- Real-time data validation
- Interactive mortality curves and probit regression plots
- PDF report generation with embedded high-resolution plots
- No installation required (browser-based)
- No user data storage (privacy-preserving)

---

## 2. TECHNICAL SPECIFICATIONS

### Technology Stack

**Frontend Framework:**
- Streamlit 1.28+ (Python web framework)
- Interactive web interface
- Real-time computation

**Backend:**
- Python 3.8 or higher
- Stateless architecture (no persistent storage)
- Session-based processing

**Core Dependencies:**
```
streamlit>=1.28.0
pandas>=1.3.0
numpy>=1.21.0
scipy>=1.7.0
statsmodels>=0.13.0
matplotlib>=3.4.0
fpdf>=1.7.2
```

**Statistical Methods:**
- Probit regression (GLM with probit link function)
- Maximum likelihood estimation (IRLS algorithm)
- Delta method for confidence intervals
- Fieller's method for resistance ratios
- Chi-square goodness-of-fit test

**File Processing:**
- Input: Tab-delimited text files (.txt)
- Output: PDF reports with embedded plots
- No database required
- No persistent storage

**Security Considerations:**
- No user authentication (public access)
- No PII collection
- No data retention
- File size limits (200 MB)
- Input validation on all user data
- HTTPS encryption (when deployed on Streamlit Cloud)

### Application Architecture

```
User Browser
    ↓ [Upload data file]
Streamlit Web Interface
    ↓ [Validate data]
Python Analysis Engine
    ├─ Data preprocessing
    ├─ Probit regression
    ├─ Statistical tests
    ├─ Plot generation
    └─ PDF creation
    ↓ [Return results]
User Browser
    ↓ [Download PDF]
Session Ends (no data retained)
```

**Data Flow:**
1. User uploads bioassay data file
2. Data validated in browser session
3. Analysis performed in memory
4. Results displayed interactively
5. Optional PDF download
6. Session ends, data discarded
7. No server-side storage

### System Requirements

**For Users (Browser-based):**
- Modern web browser (Chrome, Firefox, Safari, Edge)
- Internet connection
- No installation required
- No admin rights needed

**For Deployment:**
- Option 1: Streamlit Cloud (recommended - managed hosting)
- Option 2: Linux server with Python 3.8+
- Option 3: Docker container
- Option 4: Windows/Mac server

---

## 3. DEPLOYMENT OPTIONS

### Option 1: Streamlit Cloud (RECOMMENDED)

**Advantages:**
- Free for public applications
- Managed hosting (no server maintenance)
- Automatic HTTPS
- Auto-deployment from GitHub
- Zero downtime updates
- Professional URL

**Deployment Steps:**
1. Create GitHub repository
2. Push code to GitHub
3. Connect to Streamlit Cloud (share.streamlit.io)
4. Deploy application (automatic)
5. Get public URL: https://[app-name].streamlit.app

**Time to Deploy:** 10-15 minutes  
**Cost:** Free  
**Maintenance:** Automatic updates via GitHub

**Technical Requirements:**
- GitHub account
- Repository with:
  - probit_web_app.py
  - requirements_web.txt
  - README.md

---

### Option 2: Institutional Server

**Advantages:**
- Full control
- No external dependencies
- Can run on internal network
- Custom domain possible

**Deployment Steps:**
1. Install Python 3.8+ on server
2. Install dependencies: `pip install -r requirements_web.txt`
3. Run: `streamlit run probit_web_app.py --server.port 8501`
4. Configure reverse proxy (Nginx/Apache) for HTTPS
5. Set up systemd service for persistence

**Server Requirements:**
- Ubuntu 20.04+ or similar
- 2 GB RAM minimum
- Python 3.8+
- Port 8501 available

**Time to Deploy:** 1-2 hours  
**Maintenance:** Manual updates, server monitoring

---

### Option 3: Docker Container

**Advantages:**
- Reproducible environment
- Easy to deploy to cloud platforms
- Isolated from host system

**Dockerfile provided in deployment guide**

**Cloud Platforms:**
- Google Cloud Run
- AWS ECS
- Azure Container Instances
- DigitalOcean App Platform

**Time to Deploy:** 30-60 minutes  
**Cost:** Varies by platform (~$10-50/month)

---

### Option 4: Local Desktop Use

Users can run locally for offline analysis:

```bash
pip install -r requirements_web.txt
streamlit run probit_web_app.py
```

Opens in browser at http://localhost:8501

---

## 4. INSTALLATION INSTRUCTIONS

### Quick Start (Streamlit Cloud)

**Prerequisites:**
- GitHub account
- Git installed locally

**Step 1: Prepare Repository**
```bash
# Create new directory
mkdir probit-analysis-tool
cd probit-analysis-tool

# Initialize git
git init

# Add files
cp /path/to/probit_web_app.py .
cp /path/to/requirements_web.txt .
cp /path/to/README.md .

# Commit
git add .
git commit -m "Initial commit"

# Push to GitHub
git remote add origin https://github.com/[username]/probit-tool.git
git push -u origin main
```

**Step 2: Deploy to Streamlit Cloud**
1. Go to https://share.streamlit.io
2. Sign in with GitHub
3. Click "New app"
4. Select repository: [username]/probit-tool
5. Main file: probit_web_app.py
6. Click "Deploy"

**Step 3: Application Ready**
- URL: https://[app-name].streamlit.app
- Updates automatically when you push to GitHub
- No further maintenance required

---

### Installation Verification

**Test Checklist:**
1. ✓ Application loads without errors
2. ✓ Upload example dataset
3. ✓ Run single analysis
4. ✓ Verify plots appear
5. ✓ Download PDF report
6. ✓ Upload two datasets
7. ✓ Run comparison
8. ✓ Verify resistance ratio calculated
9. ✓ Download comparison PDF

**Example test files provided in repository**

---

## 5. USER DOCUMENTATION

### Input Data Format

**File Requirements:**
- Tab-delimited text file (.txt)
- Three header lines:
  - Line 1: Strain/population name
  - Line 2: Chemical/acaricide name
  - Line 3: Column headers (concentration, n, mortality)
- Data columns:
  - concentration: numeric values
  - n: number of individuals tested
  - mortality: number of dead individuals

**Example:**
```
Yucatan
Coumaphos
concentration	n	mortality
0.500	100	95
0.250	100	70
0.125	100	40
0.063	100	15
0.031	100	5
```

### Using the Application

**Single Dataset Analysis:**
1. Navigate to "Upload Data" tab
2. Upload data file for Dataset 1
3. Verify validation passes
4. Go to "Single Analysis" tab
5. Click "Run Analysis"
6. Review results:
   - Data summary
   - LD estimates (LD1, LD50, LD99)
   - Model fit statistics
   - Model parameters
   - Replicate variability
   - Mortality curve
   - Probit regression plot
7. Download PDF report

**Two-Dataset Comparison:**
1. Upload both datasets
2. Go to "Compare Two Datasets" tab
3. Click "Run Comparison"
4. Review results:
   - Resistance ratio with 95% CI
   - LD comparison table
   - Statistical tests (parallelism, equality)
   - Combined plots
   - Interpretation text
5. Download comparison PDF

### Interpreting Results

**LD Values:**
- LD1: Lethal dose for 1% of population
- LD50: Lethal dose for 50% of population (median)
- LD99: Lethal dose for 99% of population
- 95% CI: Confidence interval around estimate

**Model Fit:**
- Chi-square test assesses goodness-of-fit
- p > 0.05: Good fit (model appropriate)
- p < 0.05: Poor fit (check data quality/variability)

**Resistance Ratio:**
- Ratio of LD50 values
- RR > 1: Test strain more resistant
- RR < 1: Test strain more susceptible
- 95% CI not including 1.0: Statistically significant

**Replicate Variability:**
- CV% (coefficient of variation) per concentration
- CV% > 20%: High variability (highlighted in red)
- High CV% may indicate experimental issues

---

## 6. DATA PRIVACY & SECURITY

### No PII Collection
- Application does not collect, store, or transmit any personally identifiable information
- No user accounts or authentication
- No user tracking or analytics

### Data Handling
- All data processing occurs in user's browser session
- No server-side data storage
- Files uploaded by user are:
  - Processed in memory
  - Deleted when session ends
  - Never written to disk
  - Never transmitted to third parties

### Session Management
- Each user gets isolated session
- Sessions expire after inactivity
- No cross-session data access
- No persistent storage

### Network Security
- HTTPS encryption (Streamlit Cloud)
- No sensitive data transmission
- Input validation prevents injection attacks
- File size limits prevent abuse

---

## 7. MAINTENANCE & SUPPORT

### Update Procedure
**If deployed on Streamlit Cloud:**
1. Make changes to code locally
2. Commit to GitHub
3. Push to repository
4. Streamlit Cloud auto-deploys
5. Users see updates immediately

**No downtime required**

### Monitoring
- Streamlit Cloud provides basic analytics
- No additional monitoring needed for public deployment

### Known Limitations
- Maximum file size: 200 MB
- Browser-based (requires internet for cloud version)
- Not optimized for very large datasets (>10,000 observations)
- PDF plots are static (not interactive)

### Future Enhancements (Optional)
- Batch processing multiple files
- Additional plot types
- Export results to CSV/Excel
- More statistical tests
- Dose-response model comparisons

---

## 8. TESTING GUIDE

### Pre-Deployment Testing

**Test 1: Single Dataset**
1. Upload example_dataset1.txt
2. Verify validation passes
3. Run analysis
4. Check all outputs appear:
   - ✓ Data summary (4 metrics)
   - ✓ LD table (LD1, LD50, LD99)
   - ✓ Model fit (chi-square, p-value)
   - ✓ Parameters (intercept, slope)
   - ✓ Variability table
   - ✓ Mortality curve
   - ✓ Probit plot
   - ✓ PDF download works

**Test 2: Two-Dataset Comparison**
1. Upload both example datasets
2. Run comparison
3. Check outputs:
   - ✓ Resistance ratio calculated
   - ✓ Statistical tests complete
   - ✓ Combined plots appear
   - ✓ Interpretation text generated
   - ✓ PDF download works

**Test 3: Data Validation**
1. Upload file with errors
2. Verify error messages appear
3. Check warnings for high variability

**Test 4: Cross-Browser**
- ✓ Chrome
- ✓ Firefox
- ✓ Safari
- ✓ Edge

**Test 5: Mobile (Optional)**
- ✓ Loads on mobile browser
- ✓ Basic functionality works

---

## 9. CITATION & ATTRIBUTION

### Public Domain Notice
```
USDA Public Domain Notice

This software was developed by the U.S. Department of Agriculture, 
Agricultural Research Service (USDA ARS). As a work of the United 
States Government, this software is in the public domain within 
the United States.

Additionally, USDA ARS waives copyright and related rights in 
the work worldwide through the CC0 1.0 Universal Public Domain 
Dedication.

Developed by: Jason Tidwell, Research Entomologist
USDA ARS Adkisson-Pfrimmer Agricultural Genomics Center
College Station, TX
2024
```

### Citation for Publications
```
Tidwell, J. (2024). Probit Analysis Tool for Bioassay Research 
(Version 10.3). USDA Agricultural Research Service. 
Available at: [GitHub URL or Web App URL]
```

**In Methods Section:**
"Probit regression analysis was performed using the Probit Analysis 
Tool v10.3 (Tidwell, 2024; USDA ARS)."

---

## 10. SUPPORT & CONTACT

### Primary Contact
Jason Tidwell  
Research Entomologist  
USDA ARS Adkisson-Pfrimmer Agricultural Genomics Center  
College Station, TX  
Email: [your.email@usda.gov]

### Documentation
- GitHub Repository: [URL]
- User Guide: README.md in repository
- Deployment Guide: WEB_DEPLOYMENT_GUIDE.md
- Testing Guide: WEB_APP_TESTING_GUIDE.md

### Reporting Issues
Users can report issues via:
1. GitHub Issues (if repository public)
2. Direct email to contact above
3. Through Area Office

---

## 11. COMPLIANCE & APPROVALS

### ARS Review Status
- Submitted for review: [Date]
- Public release approved: [Pending]
- Public domain status: Confirmed (U.S. government work)

### Technology Transfer
- OTT Review: Not required (educational software, no commercial value)
- Patent Review: Not applicable (implements standard methods)
- Export Control: Not applicable (no restrictions)

### Open Source
- License: Public Domain (U.S. government work)
- No proprietary dependencies
- All code freely redistributable
- No usage restrictions

---

## APPENDIX A: FILE MANIFEST

**Core Application:**
- probit_web_app.py (26 KB) - Main application
- requirements_web.txt (1 KB) - Dependencies

**Documentation:**
- README_WEB.md (8 KB) - User guide
- WEB_DEPLOYMENT_GUIDE.md (15 KB) - Deployment instructions
- WEB_APP_TESTING_GUIDE.md (9 KB) - Testing procedures
- PDF_REPORTS_COMPLETE_GUIDE.md (7 KB) - PDF documentation
- CITATION_GUIDE.md (12 KB) - Citation instructions

**Example Data:**
- example_dataset1.txt (1 KB) - Yucatan strain
- example_dataset2.txt (1 KB) - Reference strain

**Total Package Size:** ~80 KB (excluding images)

---

## APPENDIX B: DEPLOYMENT CHECKLIST

### Pre-Deployment
- [ ] Code tested locally
- [ ] All documentation complete
- [ ] Example files prepared
- [ ] GitHub repository created
- [ ] Public domain notice added
- [ ] Citation information included
- [ ] ARS approval obtained

### Deployment
- [ ] Repository pushed to GitHub
- [ ] Streamlit Cloud account created
- [ ] Application deployed
- [ ] URL obtained and tested
- [ ] Cross-browser testing completed
- [ ] PDF downloads verified

### Post-Deployment
- [ ] URL shared with stakeholders
- [ ] Listed on ARS software page
- [ ] Announced to target user community
- [ ] Monitoring configured (optional)
- [ ] Feedback mechanism established

---

**END OF DEPLOYMENT PACKAGE**

**Document Version:** 1.0  
**Last Updated:** December 2024  
**Status:** Ready for Deployment
