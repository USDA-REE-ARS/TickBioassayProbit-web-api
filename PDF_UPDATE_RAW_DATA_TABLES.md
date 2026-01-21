# 📄 PDF Report Update - Raw Data Tables Added

## ✅ **What Was Updated**

The PDF reports now include **raw data tables** showing the actual data that was analyzed.

---

## 📊 **Single Dataset PDF - Now Includes:**

### **Page 1: Complete Analysis**
1. Header (strain, chemical, date)
2. Data Summary (4 metrics)
3. **NEW: Raw Data Table** ✨
   - Concentration
   - N Tested
   - Mortality  
   - Mortality %
4. LD Estimates (LD1, LD50, LD99 with CIs)
5. Model Fit Statistics
6. Model Parameters

### **Page 2: Mortality Curve**
- Full-page plot with observed vs. fitted

### **Page 3: Probit Regression Line**
- Full-page probit plot

---

## ⚖️ **Comparison PDF - Now Includes:**

### **Page 1: Analysis & Data**
1. Header (both datasets, date)
2. **NEW: Raw Data Table - Dataset 1** ✨
   - All observations with mortality percentages
3. **NEW: Raw Data Table - Dataset 2** ✨
   - All observations with mortality percentages
4. Resistance Ratio with CI
5. LD Comparison Table
6. Statistical Tests
7. Model Fit Comparison

### **Page 2: Combined Mortality Curves**
- Both strains overlaid

### **Page 3: Combined Probit Lines**
- Both regression lines

---

## 📋 **Raw Data Table Format**

### **Columns Included:**

| Concentration | N Tested | Mortality | Mortality % |
|---------------|----------|-----------|-------------|
| 0.500000      | 100      | 95        | 95.00%      |
| 0.250000      | 100      | 70        | 70.00%      |
| 0.125000      | 100      | 40        | 40.00%      |
| ...           | ...      | ...       | ...         |

### **Features:**
- ✅ All concentrations from input file
- ✅ Sample sizes (n)
- ✅ Raw mortality counts
- ✅ Calculated mortality percentages
- ✅ Professional table formatting
- ✅ Clear borders and headers

---

## 🎯 **Benefits of Including Raw Data**

### **1. Complete Documentation**
- Readers can see exactly what data was analyzed
- No need to refer back to original files
- Self-contained report

### **2. Transparency**
- Full reproducibility
- Data quality assessment
- Easy verification

### **3. Publication Ready**
- Meets journal requirements for data availability
- Supplementary material ready
- Reviewers can see raw data

### **4. Archival Quality**
- Complete record in single PDF
- Years later, you know exactly what was run
- No lost data files

---

## 📄 **Example: What the PDF Now Shows**

### **Single Dataset Report:**

```
═══════════════════════════════════════════════
         Probit Analysis Report
═══════════════════════════════════════════════

Strain: Yucatan
Chemical: Coumaphos  
Date: 2024-12-11 14:30

Data Summary
───────────────────────────────────────────────
Number of observations: 21
Total individuals tested: 2,466
Total mortality: 988
Overall mortality: 40.07%

Raw Data                              ← NEW!
───────────────────────────────────────────────
┌─────────────┬──────────┬───────────┬──────────────┐
│Concentration│ N Tested │ Mortality │ Mortality %  │
├─────────────┼──────────┼───────────┼──────────────┤
│  0.500000   │   100    │    95     │   95.00%     │
│  0.250000   │   100    │    70     │   70.00%     │
│  0.125000   │   100    │    40     │   40.00%     │
│  0.063000   │   100    │    15     │   15.00%     │
│  0.031000   │   100    │     5     │    5.00%     │
│     ...     │   ...    │    ...    │     ...      │
└─────────────┴──────────┴───────────┴──────────────┘

Lethal Dose Estimates (95% CI)
───────────────────────────────────────────────
LD1:  0.066374 (0.061360 - 0.071796)
LD50: 0.209622 (0.202819 - 0.216652)
LD99: 0.662029 (0.607538 - 0.721408)

[... rest of analysis ...]
```

---

### **Comparison Report:**

```
═══════════════════════════════════════════════
    Probit Analysis Comparison Report
═══════════════════════════════════════════════

Dataset 1: cora_f3 - coumaphos
Dataset 2: deutch - coumaphos
Date: 2024-12-11 14:35

Raw Data - Dataset 1                  ← NEW!
───────────────────────────────────────────────
┌─────────────┬──────────┬───────────┬──────────────┐
│Concentration│ N Tested │ Mortality │ Mortality %  │
├─────────────┼──────────┼───────────┼──────────────┤
│  0.500000   │   100    │    95     │   95.00%     │
│  0.250000   │   100    │    70     │   70.00%     │
│     ...     │   ...    │    ...    │     ...      │
└─────────────┴──────────┴───────────┴──────────────┘

Raw Data - Dataset 2                  ← NEW!
───────────────────────────────────────────────
┌─────────────┬──────────┬───────────┬──────────────┐
│Concentration│ N Tested │ Mortality │ Mortality %  │
├─────────────┼──────────┼───────────┼──────────────┤
│  0.500000   │   100    │    99     │   99.00%     │
│  0.250000   │   100    │    85     │   85.00%     │
│     ...     │   ...    │    ...    │     ...      │
└─────────────┴──────────┴───────────┴──────────────┘

Resistance Ratio (LD50 Basis)
───────────────────────────────────────────────
Ratio: 9.23x
95% CI: (6.89 - 12.36)

[... rest of comparison ...]
```

---

## 🔍 **Technical Details**

### **Table Formatting:**
- **Font:** Arial 9-10pt for data tables (fits more data)
- **Headers:** Bold 9pt
- **Borders:** All cells bordered for clarity
- **Alignment:** Numeric values right-aligned
- **Decimals:** Concentration (6 places), Mortality % (2 places)

### **Space Management:**
- Tables use smaller font to fit on page
- Comparison uses 9pt to accommodate two tables
- Auto-adjusts for number of rows
- Page breaks handled automatically

### **PDF File Size:**
Minimal impact on file size:
- Each data row: ~50 bytes
- 20-row table: ~1 KB
- Total increase: < 2-3 KB
- **Still well under 500 KB total**

---

## ✅ **What This Means for Users**

### **Before:**
- PDF had statistics and plots only
- Had to keep original data files
- Couldn't verify data from PDF alone

### **After:**
- PDF is self-contained ✅
- Can verify all results ✅
- Data visible for reviewers ✅
- Complete documentation ✅

---

## 📝 **For Publications**

The PDF now serves as:

1. **Methods Documentation**
   - Complete data in report
   - No separate data file needed

2. **Supplementary Material**
   - Can submit PDF as-is
   - Meets data availability requirements

3. **Reviewer Response**
   - "Raw data is included in Supplemental PDF"
   - Full transparency

4. **Archival Record**
   - Years later, you know exactly what was done
   - No missing files

---

## 🎯 **Testing the Update**

### **To Verify:**

1. Run single dataset analysis
2. Download PDF
3. **Check Page 1** → Should see raw data table after "Data Summary"
4. Verify all your concentrations are listed

For comparison:
1. Run two-dataset comparison
2. Download PDF  
3. **Check Page 1** → Should see TWO raw data tables
4. Both datasets visible before "Resistance Ratio"

---

## 📊 **Example Use Case**

**Scenario:** Submitting manuscript on acaricide resistance

**Before:**
- Manuscript: "See supplementary file data.xlsx"
- Reviewer: "Can't open Excel file, what were the actual values?"

**After:**
- Manuscript: "See supplementary PDF"
- Reviewer: Opens PDF → Sees data, plots, and statistics all in one place ✅

---

## 💡 **Additional Notes**

### **Data Preservation:**
- Original data is unchanged
- Table shows data "as uploaded"
- Mortality % calculated for convenience

### **Quality Check:**
- Users can visually verify their data
- Catch typos before sharing
- Confirm correct file was analyzed

### **Professional Appearance:**
- Clean table formatting
- Consistent with rest of PDF
- Publication quality

---

## 🎉 **Summary**

**What Changed:**
- ✅ Added raw data tables to single dataset PDF (Page 1)
- ✅ Added raw data tables to comparison PDF (Page 1, both datasets)
- ✅ Tables show: Concentration, N, Mortality, Mortality %
- ✅ Professional formatting with borders

**Benefits:**
- ✅ Self-contained reports
- ✅ Full transparency
- ✅ Publication ready
- ✅ Better documentation

**File Size Impact:**
- Minimal (< 3 KB increase)
- Still easily emailable

**User Experience:**
- ✅ Download PDF → See everything
- ✅ No need for separate data files
- ✅ Complete archival record

---

**Your PDF reports are now even more comprehensive and publication-ready!** 📄✨
