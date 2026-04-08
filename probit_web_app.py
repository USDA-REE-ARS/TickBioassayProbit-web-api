"""
Probit Analysis Tool - Web Application
Streamlit-based web interface for bioassay probit regression analysis
Version: 10.3 Web

Run with: streamlit run probit_web_app.py
"""

import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import statsmodels.api as sm
from scipy.stats import norm, chi2
from io import BytesIO, StringIO
import tempfile
import os
from datetime import datetime
from fpdf import FPDF
import base64

# Set page configuration
st.set_page_config(
    page_title="Probit Analysis Tool",
    page_icon="📊",
    layout="wide",
    initial_sidebar_state="expanded"
)

# Constants
CONTINUITY_CORRECTION = 0.5
Z_SCORE_95_CI = 1.96
CHI_SQUARE_EPSILON = 1e-10
ALPHA_LEVEL = 0.05
LD_LEVELS = [1, 50, 99]

# Custom CSS for better appearance
st.markdown("""
<style>
    .main-header {
        font-size: 2.5rem;
        font-weight: bold;
        color: #1f77b4;
        text-align: center;
        padding: 1rem 0;
    }
    .sub-header {
        font-size: 1.5rem;
        color: #2c3e50;
        margin-top: 1rem;
    }
    .success-box {
        background-color: #d4edda;
        border: 1px solid #c3e6cb;
        border-radius: 5px;
        padding: 1rem;
        margin: 1rem 0;
    }
    .warning-box {
        background-color: #fff3cd;
        border: 1px solid #ffc107;
        border-radius: 5px;
        padding: 1rem;
        margin: 1rem 0;
    }
    .error-box {
        background-color: #f8d7da;
        border: 1px solid #f5c6cb;
        border-radius: 5px;
        padding: 1rem;
        margin: 1rem 0;
    }
    .metric-card {
        background-color: #f8f9fa;
        border-radius: 10px;
        padding: 1rem;
        box-shadow: 0 2px 4px rgba(0,0,0,0.1);
    }
</style>
""", unsafe_allow_html=True)

# Helper Functions
def validate_bioassay_data(df, strain, chemical):
    """Validate bioassay data for common issues"""
    errors = []
    warnings = []
    
    # Check for required columns
    required_cols = ['concentration', 'n', 'mortality']
    missing = [col for col in required_cols if col not in df.columns]
    if missing:
        errors.append(f"Missing required columns: {', '.join(missing)}")
        return errors, warnings
    
    # Check for negative values
    if (df['n'] < 0).any() or (df['mortality'] < 0).any() or (df['concentration'] < 0).any():
        errors.append("Data contains negative values")
    
    # Check mortality > n
    invalid = df[df['mortality'] > df['n']]
    if len(invalid) > 0:
        errors.append(f"Mortality exceeds sample size in {len(invalid)} row(s)")
    
    # Check for zero concentration
    if (df['concentration'] == 0).any():
        warnings.append("Data contains zero concentrations (will be filtered)")
    
    # Check for minimum data points
    if len(df) < 3:
        errors.append("Need at least 3 data points for analysis")
    
    # Check for high mortality variability
    if len(df) > 1:
        df_check = df.copy()
        df_check['mort_pct'] = df_check['mortality'] / df_check['n'] * 100
        by_conc = df_check.groupby('concentration')['mort_pct'].agg(['std', 'count'])
        high_var = by_conc[(by_conc['std'] > 20) & (by_conc['count'] > 1)]
        if len(high_var) > 0:
            warnings.append(f"High variability (CV% > 20%) at {len(high_var)} concentration(s)")
    
    return errors, warnings

def preprocess_data(df):
    """Preprocess bioassay data"""
    # Remove rows with n=0 or concentration=0
    df = df[(df['n'] > 0) & (df['concentration'] > 0)].copy()
    
    # Apply Abbott correction for boundary values
    df['mortality_adj'] = df['mortality'].copy()
    df['mortality_adj'] = df['mortality_adj'].clip(lower=CONTINUITY_CORRECTION)
    df['mortality_adj'] = df['mortality_adj'].clip(upper=df['n'] - CONTINUITY_CORRECTION)
    
    # Calculate adjusted alive count
    df['alive_adj'] = df['n'] - df['mortality_adj']
    
    # Log transform concentration
    df['log_concentration'] = np.log10(df['concentration'])
    
    return df

def fit_probit_model(df):
    """Fit probit regression model with robust convergence checking"""
    # Prepare data
    X = sm.add_constant(df['log_concentration'])
    y = df[['mortality_adj', 'alive_adj']].values
    
    # Fit model with error handling
    try:
        model = sm.GLM(y, X, family=sm.families.Binomial(link=sm.families.links.Probit()))
        result = model.fit()
        
        # CRITICAL FIX: Check model convergence
        if not result.converged:
            raise RuntimeError("Model did not converge. This may indicate:\n" +
                             "• Insufficient data points\n" +
                             "• Poor concentration selection\n" +
                             "• High variability between replicates\n" +
                             "Try using more data points or better concentration spacing.")
        
        # CRITICAL FIX: Validate model parameters
        if not np.isfinite(result.params).all():
            raise RuntimeError("Model produced non-finite parameters. This may indicate:\n" +
                             "• Extreme concentration values\n" +
                             "• Complete separation (all 0% or 100% mortality)\n" +
                             "• Numerical instability\n" +
                             "Check data for extreme values.")
        
        # CRITICAL FIX: Check parameter standard errors
        if not np.isfinite(result.bse).all() or (result.bse == 0).any():
            raise RuntimeError("Invalid parameter standard errors. This may indicate:\n" +
                             "• Model fitting problems\n" +
                             "• Insufficient variability in data\n" +
                             "• Numerical issues with covariance matrix")
    
    except Exception as e:
        if isinstance(e, RuntimeError):
            raise  # Re-raise our custom error messages
        else:
            raise RuntimeError(f"Model fitting failed: {str(e)}\n" +
                             "This may be due to data quality issues. Please check:\n" +
                             "• Data format and values\n" +
                             "• Concentration range selection\n" +
                             "• Mortality counts vs. sample sizes")
    
    # Goodness of fit test with validation
    observed = df['mortality'].values
    n_values = df['n'].values
    predicted_prob = result.predict(X)
    
    # CRITICAL FIX: Validate predictions
    if not np.isfinite(predicted_prob).all() or (predicted_prob < 0).any() or (predicted_prob > 1).any():
        raise RuntimeError("Model produced invalid probability predictions. Check model parameters.")
    
    expected = predicted_prob * n_values
    
    # Avoid division by zero in chi-square calculation
    chi_square = np.sum(((observed - expected) ** 2) / (expected + CHI_SQUARE_EPSILON))
    df_degrees = len(df) - 2
    
    if df_degrees <= 0:
        raise RuntimeError("Insufficient degrees of freedom for goodness-of-fit test. Need at least 3 data points.")
    
    p_value = 1 - chi2.cdf(chi_square, df_degrees)
    
    return result, chi_square, df_degrees, p_value

def compute_ldx_with_ci(result, x, alpha=ALPHA_LEVEL):
    """Compute LDx with confidence intervals with robust error checking"""
    intercept, slope = result.params
    
    # CRITICAL FIX: Validate model parameters
    if not np.isfinite(intercept) or not np.isfinite(slope):
        raise ValueError("Model parameters are not finite. Check data quality and model convergence.")
    
    # CRITICAL FIX: Check for problematic slopes that would cause division by zero
    if abs(slope) < 1e-10:
        raise ValueError("Slope is too close to zero (flat dose-response curve). " +
                        "Check if concentrations span an appropriate range for the mortality response.")
    
    target_quantile = norm.ppf(x / 100.0)
    
    # Point estimate with validation
    log_conc_ld = (target_quantile - intercept) / slope
    ld = 10 ** log_conc_ld
    
    # CRITICAL FIX: Validate LD estimate
    if not np.isfinite(ld) or ld <= 0:
        raise ValueError(f"Invalid LD{x} estimate: {ld}. " +
                        "This may indicate model fitting problems or extreme parameter values.")
    
    # Confidence interval using delta method with validation
    cov_matrix = result.cov_params()
    var_intercept = cov_matrix.iloc[0, 0]
    var_slope = cov_matrix.iloc[1, 1]
    cov_intercept_slope = cov_matrix.iloc[0, 1]
    
    # CRITICAL FIX: Validate covariance matrix values
    if not (np.isfinite(var_intercept) and np.isfinite(var_slope) and np.isfinite(cov_intercept_slope)):
        raise ValueError("Invalid covariance matrix values. Check model fitting.")
    
    grad_intercept = -1 / slope
    grad_slope = -(target_quantile - intercept) / (slope ** 2)
    
    var_log_ld = (grad_intercept ** 2) * var_intercept + \
                  (grad_slope ** 2) * var_slope + \
                  2 * grad_intercept * grad_slope * cov_intercept_slope
    
    # CRITICAL FIX: Check for negative variance
    if var_log_ld < 0:
        raise ValueError("Negative variance calculated in LD estimation. Check model covariance matrix.")
    
    se_log_ld = np.sqrt(var_log_ld)
    z_critical = norm.ppf(1 - alpha / 2)
    
    log_ld_lower = log_conc_ld - z_critical * se_log_ld
    log_ld_upper = log_conc_ld + z_critical * se_log_ld
    
    ld_lower = 10 ** log_ld_lower
    ld_upper = 10 ** log_ld_upper
    
    # CRITICAL FIX: Validate confidence interval bounds
    if not (np.isfinite(ld_lower) and np.isfinite(ld_upper) and ld_lower > 0 and ld_upper > 0):
        raise ValueError(f"Invalid confidence interval bounds: [{ld_lower}, {ld_upper}]")
    
    return ld, ld_lower, ld_upper

def calculate_replicate_variability(df):
    """Calculate variability statistics for replicates"""
    df_with_pct = df.copy()
    df_with_pct['mortality_pct'] = (df_with_pct['mortality'] / df_with_pct['n']) * 100
    
    variability = df_with_pct.groupby('concentration').agg({
        'mortality_pct': ['mean', 'std', 'count'],
        'n': 'sum',
        'mortality': 'sum'
    }).reset_index()
    
    variability.columns = ['concentration', 'mean_mortality_pct', 'std_mortality_pct', 
                          'n_replicates', 'total_n', 'total_mortality']
    
    variability['cv_pct'] = (variability['std_mortality_pct'] / variability['mean_mortality_pct']) * 100
    variability['cv_pct'] = variability['cv_pct'].fillna(0)
    variability['std_mortality_pct'] = variability['std_mortality_pct'].fillna(0)
    
    return variability

def create_mortality_plot(df, result, strain, chemical):
    """Create mortality curve plot"""
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # Observed data
    mortality_pct = (df['mortality'] / df['n']) * 100
    ax.scatter(df['concentration'], mortality_pct, 
              color='red', s=100, alpha=0.6, label=f'{strain} (Observed)', zorder=3)
    
    # Fitted curve
    conc_range = np.logspace(np.log10(df['concentration'].min()), 
                             np.log10(df['concentration'].max()), 100)
    log_conc_range = np.log10(conc_range)
    X_pred = sm.add_constant(log_conc_range)
    pred_prob = result.predict(X_pred)
    pred_mortality_pct = pred_prob * 100
    
    ax.plot(conc_range, pred_mortality_pct, 
           color='blue', linewidth=2, label=f'{strain} (Fitted)', zorder=2)
    
    ax.set_xlabel(f'{chemical} Concentration', fontsize=12, fontweight='bold')
    ax.set_ylabel('Mortality (%)', fontsize=12, fontweight='bold')
    ax.set_title(f'Observed vs. Fitted Mortality Curve\n{strain} - {chemical}', 
                fontsize=14, fontweight='bold')
    ax.set_xscale('log')
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=10)
    ax.set_ylim(-5, 105)
    
    plt.tight_layout()
    return fig

def create_probit_plot(df, result, strain, chemical):
    """Create probit transformation plot"""
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # Calculate empirical probits
    mortality_pct = (df['mortality'] / df['n']) * 100
    # Adjust for 0 and 100% using Abbott correction
    mortality_pct_adj = mortality_pct.copy()
    mortality_pct_adj[mortality_pct_adj <= 0] = 0.1
    mortality_pct_adj[mortality_pct_adj >= 100] = 99.9
    empirical_probits = norm.ppf(mortality_pct_adj / 100)
    
    # Observed data (probit scale)
    ax.scatter(df['log_concentration'], empirical_probits,
              color='red', s=100, alpha=0.6, label=f'{strain} (Observed)', zorder=3)
    
    # Fitted line
    log_conc_range = np.linspace(df['log_concentration'].min(),
                                 df['log_concentration'].max(), 100)
    fitted_probits = result.params[0] + result.params[1] * log_conc_range
    
    ax.plot(log_conc_range, fitted_probits,
           color='blue', linewidth=2, label=f'{strain} (Fitted)', zorder=2)
    
    # Reference lines for LD levels
    for ld_level, label in [(0.01, 'LD1'), (0.5, 'LD50'), (0.99, 'LD99')]:
        probit_value = norm.ppf(ld_level)
        ax.axhline(y=probit_value, color='gray', linestyle='--', alpha=0.3, linewidth=1)
        ax.text(df['log_concentration'].max(), probit_value, f' {label}',
               verticalalignment='center', fontsize=9, alpha=0.7)
    
    ax.set_xlabel(f'Log10({chemical} Concentration)', fontsize=12, fontweight='bold')
    ax.set_ylabel('Probit (Mortality)', fontsize=12, fontweight='bold')
    ax.set_title(f'Probit Regression Line\n{strain} - {chemical}',
                fontsize=14, fontweight='bold')
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=10)
    
    plt.tight_layout()
    return fig

def create_combined_probit_plot(df1, result1, strain1, df2, result2, strain2, chemical):
    """Create combined probit plot for two datasets"""
    fig, ax = plt.subplots(figsize=(12, 7))
    
    # Dataset 1 - empirical probits
    mortality_pct1 = (df1['mortality'] / df1['n']) * 100
    mortality_pct1_adj = mortality_pct1.copy()
    mortality_pct1_adj[mortality_pct1_adj <= 0] = 0.1
    mortality_pct1_adj[mortality_pct1_adj >= 100] = 99.9
    empirical_probits1 = norm.ppf(mortality_pct1_adj / 100)
    
    ax.scatter(df1['log_concentration'], empirical_probits1,
              color='red', s=100, alpha=0.6, label=f'{strain1} (Observed)', zorder=3)
    
    # Dataset 1 - fitted line
    log_conc_range1 = np.linspace(df1['log_concentration'].min(),
                                  df1['log_concentration'].max(), 100)
    fitted_probits1 = result1.params[0] + result1.params[1] * log_conc_range1
    ax.plot(log_conc_range1, fitted_probits1,
           color='darkred', linewidth=2, label=f'{strain1} (Fitted)', zorder=2)
    
    # Dataset 2 - empirical probits
    mortality_pct2 = (df2['mortality'] / df2['n']) * 100
    mortality_pct2_adj = mortality_pct2.copy()
    mortality_pct2_adj[mortality_pct2_adj <= 0] = 0.1
    mortality_pct2_adj[mortality_pct2_adj >= 100] = 99.9
    empirical_probits2 = norm.ppf(mortality_pct2_adj / 100)
    
    ax.scatter(df2['log_concentration'], empirical_probits2,
              color='blue', s=100, alpha=0.6, label=f'{strain2} (Observed)', zorder=3)
    
    # Dataset 2 - fitted line
    log_conc_range2 = np.linspace(df2['log_concentration'].min(),
                                  df2['log_concentration'].max(), 100)
    fitted_probits2 = result2.params[0] + result2.params[1] * log_conc_range2
    ax.plot(log_conc_range2, fitted_probits2,
           color='darkblue', linewidth=2, label=f'{strain2} (Fitted)', zorder=2)
    
    # Reference lines for LD50
    probit_50 = norm.ppf(0.5)
    ax.axhline(y=probit_50, color='gray', linestyle='--', alpha=0.5, linewidth=1)
    
    # Calculate LD50 positions
    ld50_1 = (probit_50 - result1.params[0]) / result1.params[1]
    ld50_2 = (probit_50 - result2.params[0]) / result2.params[1]
    ax.axvline(x=ld50_1, color='red', linestyle=':', alpha=0.5, linewidth=1)
    ax.axvline(x=ld50_2, color='blue', linestyle=':', alpha=0.5, linewidth=1)
    
    ax.set_xlabel(f'Log10({chemical} Concentration)', fontsize=12, fontweight='bold')
    ax.set_ylabel('Probit (Mortality)', fontsize=12, fontweight='bold')
    ax.set_title(f'Probit Regression Comparison\n{strain1} vs {strain2}',
                fontsize=14, fontweight='bold')
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=10, loc='best')
    
    plt.tight_layout()
    return fig

def create_pdf_single_dataset(df, result, chi_square, df_degrees, p_value, strain, chemical):
    """Create PDF report for single dataset analysis"""
    pdf = FPDF()
    pdf.set_auto_page_break(auto=True, margin=15)
    pdf.add_page()
    
    # Title
    pdf.set_font('Arial', 'B', 16)
    pdf.cell(0, 10, 'Probit Analysis Report', 0, 1, 'C')
    pdf.ln(5)
    
    # Dataset info
    pdf.set_font('Arial', 'B', 12)
    pdf.cell(0, 8, f'Strain: {strain}', 0, 1)
    pdf.cell(0, 8, f'Chemical: {chemical}', 0, 1)
    pdf.cell(0, 8, f'Date: {datetime.now().strftime("%Y-%m-%d %H:%M")}', 0, 1)
    pdf.ln(5)
    
    # Data Summary
    pdf.set_font('Arial', 'B', 14)
    pdf.cell(0, 10, 'Data Summary', 0, 1)
    pdf.set_font('Arial', '', 11)
    pdf.cell(0, 6, f'Number of observations: {len(df)}', 0, 1)
    pdf.cell(0, 6, f'Total individuals tested: {int(df["n"].sum())}', 0, 1)
    pdf.cell(0, 6, f'Total mortality: {int(df["mortality"].sum())}', 0, 1)
    overall_mort = (df['mortality'].sum() / df['n'].sum() * 100)
    pdf.cell(0, 6, f'Overall mortality: {overall_mort:.2f}%', 0, 1)
    pdf.ln(5)
    
    # LD Estimates
    pdf.set_font('Arial', 'B', 14)
    pdf.cell(0, 10, 'Lethal Dose Estimates (95% CI)', 0, 1)
    pdf.set_font('Arial', '', 11)
    for ld_level in LD_LEVELS:
        ld, lower, upper = compute_ldx_with_ci(result, ld_level)
        pdf.cell(0, 6, f'LD{ld_level}: {ld:.6f} ({lower:.6f} - {upper:.6f})', 0, 1)
    pdf.ln(5)
    
    # Model Fit
    pdf.set_font('Arial', 'B', 14)
    pdf.cell(0, 10, 'Model Fit Statistics', 0, 1)
    pdf.set_font('Arial', '', 11)
    pdf.cell(0, 6, f'Chi-Square: {chi_square:.4f}', 0, 1)
    pdf.cell(0, 6, f'Degrees of Freedom: {df_degrees}', 0, 1)
    pdf.cell(0, 6, f'p-value: {p_value:.4f}', 0, 1)
    fit_text = 'Good fit (p > 0.05)' if p_value >= 0.05 else 'Poor fit (p < 0.05)'
    pdf.cell(0, 6, f'Interpretation: {fit_text}', 0, 1)
    pdf.ln(5)
    
    # Model Parameters
    pdf.set_font('Arial', 'B', 14)
    pdf.cell(0, 10, 'Model Parameters', 0, 1)
    pdf.set_font('Arial', '', 11)
    pdf.cell(0, 6, f'Intercept: {result.params[0]:.4f} (SE: {result.bse[0]:.4f})', 0, 1)
    pdf.cell(0, 6, f'Slope: {result.params[1]:.4f} (SE: {result.bse[1]:.4f})', 0, 1)
    
    return pdf

def create_pdf_single_with_plots(df, df_processed, result, chi_square, df_degrees, p_value, strain, chemical, options=None):
    """Create PDF report with plots for single dataset analysis
    
    Args:
        options: Dict with boolean flags for what to include in PDF
                Keys: include_data_summary, include_raw_data, include_ld_estimates,
                      include_model_fit, include_parameters, include_mortality_plot, include_probit_plot
    """
    # Default options - include everything
    if options is None:
        options = {
            'include_data_summary': True,
            'include_raw_data': True,
            'include_ld_estimates': True,
            'include_model_fit': True,
            'include_parameters': True,
            'include_mortality_plot': True,
            'include_probit_plot': True
        }
    
    pdf = FPDF()
    pdf.set_auto_page_break(auto=True, margin=15)
    pdf.add_page()
    
    # Title
    pdf.set_font('Arial', 'B', 16)
    pdf.cell(0, 10, 'Probit Analysis Report', 0, 1, 'C')
    pdf.ln(5)
    
    # Dataset info
    pdf.set_font('Arial', 'B', 12)
    pdf.cell(0, 8, f'Strain: {strain}', 0, 1)
    pdf.cell(0, 8, f'Chemical: {chemical}', 0, 1)
    pdf.cell(0, 8, f'Date: {datetime.now().strftime("%Y-%m-%d %H:%M")}', 0, 1)
    pdf.ln(5)
    
    # Data Summary
    if options.get('include_data_summary', True):
        pdf.set_font('Arial', 'B', 14)
        pdf.cell(0, 10, 'Data Summary', 0, 1)
        pdf.set_font('Arial', '', 11)
        pdf.cell(0, 6, f'Number of observations: {len(df)}', 0, 1)
        pdf.cell(0, 6, f'Total individuals tested: {int(df["n"].sum())}', 0, 1)
        pdf.cell(0, 6, f'Total mortality: {int(df["mortality"].sum())}', 0, 1)
        overall_mort = (df['mortality'].sum() / df['n'].sum() * 100)
        pdf.cell(0, 6, f'Overall mortality: {overall_mort:.2f}%', 0, 1)
        pdf.ln(5)
    
    # Raw Data Table
    if options.get('include_raw_data', True):
        pdf.set_font('Arial', 'B', 14)
        pdf.cell(0, 10, 'Raw Data', 0, 1)
        pdf.set_font('Arial', 'B', 10)
        pdf.cell(45, 6, 'Concentration', 1)
        pdf.cell(35, 6, 'N Tested', 1)
        pdf.cell(35, 6, 'Mortality', 1)
        pdf.cell(45, 6, 'Mortality %', 1)
        pdf.ln()
        
        pdf.set_font('Arial', '', 10)
        for idx, row in df.iterrows():
            mort_pct = (row['mortality'] / row['n']) * 100
            pdf.cell(45, 6, f"{row['concentration']:.6f}", 1)
            pdf.cell(35, 6, f"{int(row['n'])}", 1)
            pdf.cell(35, 6, f"{int(row['mortality'])}", 1)
            pdf.cell(45, 6, f"{mort_pct:.2f}%", 1)
            pdf.ln()
        pdf.ln(5)
    
    # LD Estimates
    if options.get('include_ld_estimates', True):
        pdf.set_font('Arial', 'B', 14)
        pdf.cell(0, 10, 'Lethal Dose Estimates (95% CI)', 0, 1)
        pdf.set_font('Arial', '', 11)
        for ld_level in LD_LEVELS:
            ld, lower, upper = compute_ldx_with_ci(result, ld_level)
            pdf.cell(0, 6, f'LD{ld_level}: {ld:.6f} ({lower:.6f} - {upper:.6f})', 0, 1)
        pdf.ln(5)
    
    # Model Fit
    if options.get('include_model_fit', True):
        pdf.set_font('Arial', 'B', 14)
        pdf.cell(0, 10, 'Model Fit Statistics', 0, 1)
        pdf.set_font('Arial', '', 11)
        pdf.cell(0, 6, f'Chi-Square: {chi_square:.4f}', 0, 1)
        pdf.cell(0, 6, f'Degrees of Freedom: {df_degrees}', 0, 1)
        pdf.cell(0, 6, f'p-value: {p_value:.4f}', 0, 1)
        fit_text = 'Good fit (p > 0.05)' if p_value >= 0.05 else 'Poor fit (p < 0.05)'
        pdf.cell(0, 6, f'Interpretation: {fit_text}', 0, 1)
        pdf.ln(5)
    
    # Model Parameters
    if options.get('include_parameters', True):
        pdf.set_font('Arial', 'B', 14)
        pdf.cell(0, 10, 'Model Parameters', 0, 1)
        pdf.set_font('Arial', '', 11)
        pdf.cell(0, 6, f'Intercept: {result.params[0]:.4f} (SE: {result.bse[0]:.4f})', 0, 1)
        pdf.cell(0, 6, f'Slope: {result.params[1]:.4f} (SE: {result.bse[1]:.4f})', 0, 1)
        pdf.ln(10)
    
    # Add plots on separate pages (if any are selected)
    if options.get('include_mortality_plot', True) or options.get('include_probit_plot', True):
        pdf.add_page()
    
    # Mortality curve
    if options.get('include_mortality_plot', True):
        pdf.set_font('Arial', 'B', 14)
        pdf.cell(0, 10, 'Mortality Curve', 0, 1)
        pdf.ln(2)
        
        with tempfile.NamedTemporaryFile(delete=False, suffix='.png') as tmpfile:
            fig = create_mortality_plot(df_processed, result, strain, chemical)
            fig.savefig(tmpfile.name, format='png', dpi=150, bbox_inches='tight')
            plt.close(fig)
            pdf.image(tmpfile.name, x=10, w=190)
            os.unlink(tmpfile.name)
        
        pdf.ln(5)
    
    # Probit plot
    if options.get('include_probit_plot', True):
        # Add new page if mortality plot was included, otherwise use current page
        if options.get('include_mortality_plot', True):
            pdf.add_page()
        
        pdf.set_font('Arial', 'B', 14)
        pdf.cell(0, 10, 'Probit Regression Line', 0, 1)
        pdf.ln(2)
        
        with tempfile.NamedTemporaryFile(delete=False, suffix='.png') as tmpfile:
            fig = create_probit_plot(df_processed, result, strain, chemical)
            fig.savefig(tmpfile.name, format='png', dpi=150, bbox_inches='tight')
            plt.close(fig)
            pdf.image(tmpfile.name, x=10, w=190)
            os.unlink(tmpfile.name)
    
    return pdf

def create_pdf_comparison(df1, result1, chi_sq1, df_deg1, p_val1, strain1, chemical1,
                          df2, result2, chi_sq2, df_deg2, p_val2, strain2, chemical2,
                          resistance_ratio, rr_lower, rr_upper, p_parallel, p_equality):
    """Create PDF report for two-dataset comparison"""
    pdf = FPDF()
    pdf.set_auto_page_break(auto=True, margin=15)
    pdf.add_page()
    
    # Title
    pdf.set_font('Arial', 'B', 16)
    pdf.cell(0, 10, 'Probit Analysis Comparison Report', 0, 1, 'C')
    pdf.ln(5)
    
    # Dataset info
    pdf.set_font('Arial', 'B', 12)
    pdf.cell(0, 8, f'Dataset 1: {strain1} - {chemical1}', 0, 1)
    pdf.cell(0, 8, f'Dataset 2: {strain2} - {chemical2}', 0, 1)
    pdf.cell(0, 8, f'Date: {datetime.now().strftime("%Y-%m-%d %H:%M")}', 0, 1)
    pdf.ln(5)
    
    # Resistance Ratio
    pdf.set_font('Arial', 'B', 14)
    pdf.cell(0, 10, 'Resistance Ratio (LD50 Basis)', 0, 1)
    pdf.set_font('Arial', '', 11)
    pdf.cell(0, 6, f'Ratio: {resistance_ratio:.2f}x', 0, 1)
    pdf.cell(0, 6, f'95% Confidence Interval: ({rr_lower:.2f} - {rr_upper:.2f})', 0, 1)
    interpretation = f'{strain1} is {resistance_ratio:.2f}x {"more resistant" if resistance_ratio > 1 else "more susceptible"} than {strain2}'
    pdf.multi_cell(0, 6, f'Interpretation: {interpretation}')
    pdf.ln(5)
    
    # LD Estimates Comparison
    pdf.set_font('Arial', 'B', 14)
    pdf.cell(0, 10, 'Lethal Dose Estimates', 0, 1)
    pdf.set_font('Arial', 'B', 11)
    pdf.cell(30, 6, 'LD Level', 1)
    pdf.cell(45, 6, strain1, 1)
    pdf.cell(45, 6, strain2, 1)
    pdf.cell(30, 6, 'Ratio', 1)
    pdf.ln()
    
    pdf.set_font('Arial', '', 11)
    for ld_level in LD_LEVELS:
        ld1, _, _ = compute_ldx_with_ci(result1, ld_level)
        ld2, _, _ = compute_ldx_with_ci(result2, ld_level)
        ratio = ld1 / ld2
        pdf.cell(30, 6, f'LD{ld_level}', 1)
        pdf.cell(45, 6, f'{ld1:.6f}', 1)
        pdf.cell(45, 6, f'{ld2:.6f}', 1)
        pdf.cell(30, 6, f'{ratio:.2f}', 1)
        pdf.ln()
    pdf.ln(5)
    
    # Statistical Tests
    pdf.set_font('Arial', 'B', 14)
    pdf.cell(0, 10, 'Statistical Tests', 0, 1)
    pdf.set_font('Arial', '', 11)
    pdf.cell(0, 6, f'Parallelism Test p-value: {p_parallel:.4f}', 0, 1)
    parallel_text = 'Lines are parallel (p > 0.05)' if p_parallel > 0.05 else 'Lines differ (p < 0.05)'
    pdf.cell(0, 6, f'  {parallel_text}', 0, 1)
    pdf.cell(0, 6, f'Equality Test p-value: {p_equality:.4f}', 0, 1)
    equality_text = 'No significant difference (p > 0.05)' if p_equality > 0.05 else 'Significant difference (p < 0.05)'
    pdf.cell(0, 6, f'  {equality_text}', 0, 1)
    pdf.ln(5)
    
    # Model Fit Comparison
    pdf.set_font('Arial', 'B', 14)
    pdf.cell(0, 10, 'Model Fit Statistics', 0, 1)
    pdf.set_font('Arial', 'B', 11)
    pdf.cell(60, 6, 'Dataset', 1)
    pdf.cell(40, 6, 'Chi-Square', 1)
    pdf.cell(20, 6, 'df', 1)
    pdf.cell(40, 6, 'p-value', 1)
    pdf.ln()
    
    pdf.set_font('Arial', '', 11)
    pdf.cell(60, 6, strain1, 1)
    pdf.cell(40, 6, f'{chi_sq1:.4f}', 1)
    pdf.cell(20, 6, f'{df_deg1}', 1)
    pdf.cell(40, 6, f'{p_val1:.4f}', 1)
    pdf.ln()
    
    pdf.cell(60, 6, strain2, 1)
    pdf.cell(40, 6, f'{chi_sq2:.4f}', 1)
    pdf.cell(20, 6, f'{df_deg2}', 1)
    pdf.cell(40, 6, f'{p_val2:.4f}', 1)
    pdf.ln()
    
    return pdf

def create_pdf_comparison_with_plots(df1, df1_processed, result1, chi_sq1, df_deg1, p_val1, strain1, chemical1,
                                     df2, df2_processed, result2, chi_sq2, df_deg2, p_val2, strain2, chemical2,
                                     resistance_ratio, rr_lower, rr_upper, p_parallel, p_equality, options=None):
    """Create PDF report with plots for two-dataset comparison
    
    Args:
        options: Dict with boolean flags for what to include in PDF
                Keys: include_raw_data, include_rr, include_ld_comparison,
                      include_statistical_tests, include_model_fit,
                      include_mortality_plot, include_probit_plot
    """
    # Default options - include everything
    if options is None:
        options = {
            'include_raw_data': True,
            'include_rr': True,
            'include_ld_comparison': True,
            'include_statistical_tests': True,
            'include_model_fit': True,
            'include_mortality_plot': True,
            'include_probit_plot': True
        }
    
    pdf = FPDF()
    pdf.set_auto_page_break(auto=True, margin=15)
    pdf.add_page()
    
    # Title
    pdf.set_font('Arial', 'B', 16)
    pdf.cell(0, 10, 'Probit Analysis Comparison Report', 0, 1, 'C')
    pdf.ln(5)
    
    # Dataset info
    pdf.set_font('Arial', 'B', 12)
    pdf.cell(0, 8, f'Dataset 1: {strain1} - {chemical1}', 0, 1)
    pdf.cell(0, 8, f'Dataset 2: {strain2} - {chemical2}', 0, 1)
    pdf.cell(0, 8, f'Date: {datetime.now().strftime("%Y-%m-%d %H:%M")}', 0, 1)
    pdf.ln(5)
    
    # Raw Data Tables
    if options.get('include_raw_data', True):
        pdf.set_font('Arial', 'B', 14)
        pdf.cell(0, 10, 'Raw Data - Dataset 1', 0, 1)
        pdf.set_font('Arial', 'B', 9)
        pdf.cell(40, 6, 'Concentration', 1)
        pdf.cell(30, 6, 'N Tested', 1)
        pdf.cell(30, 6, 'Mortality', 1)
        pdf.cell(40, 6, 'Mortality %', 1)
        pdf.ln()
        
        pdf.set_font('Arial', '', 9)
        for idx, row in df1.iterrows():
            mort_pct = (row['mortality'] / row['n']) * 100
            pdf.cell(40, 5, f"{row['concentration']:.6f}", 1)
            pdf.cell(30, 5, f"{int(row['n'])}", 1)
            pdf.cell(30, 5, f"{int(row['mortality'])}", 1)
            pdf.cell(40, 5, f"{mort_pct:.2f}%", 1)
            pdf.ln()
        pdf.ln(5)
        
        pdf.set_font('Arial', 'B', 14)
        pdf.cell(0, 10, 'Raw Data - Dataset 2', 0, 1)
        pdf.set_font('Arial', 'B', 9)
        pdf.cell(40, 6, 'Concentration', 1)
        pdf.cell(30, 6, 'N Tested', 1)
        pdf.cell(30, 6, 'Mortality', 1)
        pdf.cell(40, 6, 'Mortality %', 1)
        pdf.ln()
        
        pdf.set_font('Arial', '', 9)
        for idx, row in df2.iterrows():
            mort_pct = (row['mortality'] / row['n']) * 100
            pdf.cell(40, 5, f"{row['concentration']:.6f}", 1)
            pdf.cell(30, 5, f"{int(row['n'])}", 1)
            pdf.cell(30, 5, f"{int(row['mortality'])}", 1)
            pdf.cell(40, 5, f"{mort_pct:.2f}%", 1)
            pdf.ln()
        pdf.ln(5)
    
    # Resistance Ratio
    if options.get('include_rr', True):
        pdf.set_font('Arial', 'B', 14)
        pdf.cell(0, 10, 'Resistance Ratio (LD50 Basis)', 0, 1)
        pdf.set_font('Arial', '', 11)
        pdf.cell(0, 6, f'Ratio: {resistance_ratio:.2f}x', 0, 1)
        pdf.cell(0, 6, f'95% Confidence Interval: ({rr_lower:.2f} - {rr_upper:.2f})', 0, 1)
        interpretation = f'{strain1} is {resistance_ratio:.2f}x {"more resistant" if resistance_ratio > 1 else "more susceptible"} than {strain2}'
        pdf.multi_cell(0, 6, f'Interpretation: {interpretation}')
        pdf.ln(5)
    
    # LD Estimates Comparison
    pdf.set_font('Arial', 'B', 14)
    pdf.cell(0, 10, 'Lethal Dose Estimates', 0, 1)
    pdf.set_font('Arial', 'B', 11)
    pdf.cell(30, 6, 'LD Level', 1)
    pdf.cell(45, 6, strain1, 1)
    pdf.cell(45, 6, strain2, 1)
    pdf.cell(30, 6, 'Ratio', 1)
    pdf.ln()
    
    pdf.set_font('Arial', '', 11)
    for ld_level in LD_LEVELS:
        ld1, _, _ = compute_ldx_with_ci(result1, ld_level)
        ld2, _, _ = compute_ldx_with_ci(result2, ld_level)
        ratio = ld1 / ld2
        pdf.cell(30, 6, f'LD{ld_level}', 1)
        pdf.cell(45, 6, f'{ld1:.6f}', 1)
        pdf.cell(45, 6, f'{ld2:.6f}', 1)
        pdf.cell(30, 6, f'{ratio:.2f}', 1)
        pdf.ln()
    pdf.ln(5)
    
    # Statistical Tests
    pdf.set_font('Arial', 'B', 14)
    pdf.cell(0, 10, 'Statistical Tests', 0, 1)
    pdf.set_font('Arial', '', 11)
    pdf.cell(0, 6, f'Parallelism Test p-value: {p_parallel:.4f}', 0, 1)
    parallel_text = 'Lines are parallel (p > 0.05)' if p_parallel > 0.05 else 'Lines differ (p < 0.05)'
    pdf.cell(0, 6, f'  {parallel_text}', 0, 1)
    pdf.cell(0, 6, f'Equality Test p-value: {p_equality:.4f}', 0, 1)
    equality_text = 'No significant difference (p > 0.05)' if p_equality > 0.05 else 'Significant difference (p < 0.05)'
    pdf.cell(0, 6, f'  {equality_text}', 0, 1)
    pdf.ln(5)
    
    # Model Fit Comparison
    pdf.set_font('Arial', 'B', 14)
    pdf.cell(0, 10, 'Model Fit Statistics', 0, 1)
    pdf.set_font('Arial', 'B', 11)
    pdf.cell(60, 6, 'Dataset', 1)
    pdf.cell(40, 6, 'Chi-Square', 1)
    pdf.cell(20, 6, 'df', 1)
    pdf.cell(40, 6, 'p-value', 1)
    pdf.ln()
    
    pdf.set_font('Arial', '', 11)
    pdf.cell(60, 6, strain1, 1)
    pdf.cell(40, 6, f'{chi_sq1:.4f}', 1)
    pdf.cell(20, 6, f'{df_deg1}', 1)
    pdf.cell(40, 6, f'{p_val1:.4f}', 1)
    pdf.ln()
    
    pdf.cell(60, 6, strain2, 1)
    pdf.cell(40, 6, f'{chi_sq2:.4f}', 1)
    pdf.cell(20, 6, f'{df_deg2}', 1)
    pdf.cell(40, 6, f'{p_val2:.4f}', 1)
    pdf.ln()
    pdf.ln(10)
    
    # Add plots
    pdf.add_page()
    
    # Combined mortality curve
    pdf.set_font('Arial', 'B', 14)
    pdf.cell(0, 10, 'Combined Mortality Curves', 0, 1)
    pdf.ln(2)
    
    with tempfile.NamedTemporaryFile(delete=False, suffix='.png') as tmpfile:
        # Create combined mortality plot
        fig, ax = plt.subplots(figsize=(12, 7))
        
        # Dataset 1
        mortality_pct1 = (df1_processed['mortality'] / df1_processed['n']) * 100
        ax.scatter(df1_processed['concentration'], mortality_pct1,
                  color='red', s=100, alpha=0.6, label=f'{strain1} (Observed)', zorder=3)
        
        conc_range1 = np.logspace(np.log10(df1_processed['concentration'].min()),
                                  np.log10(df1_processed['concentration'].max()), 100)
        log_conc_range1 = np.log10(conc_range1)
        X_pred1 = sm.add_constant(log_conc_range1)
        pred_prob1 = result1.predict(X_pred1)
        pred_mortality_pct1 = pred_prob1 * 100
        
        ax.plot(conc_range1, pred_mortality_pct1,
               color='darkred', linewidth=2, label=f'{strain1} (Fitted)', zorder=2)
        
        # Dataset 2
        mortality_pct2 = (df2_processed['mortality'] / df2_processed['n']) * 100
        ax.scatter(df2_processed['concentration'], mortality_pct2,
                  color='blue', s=100, alpha=0.6, label=f'{strain2} (Observed)', zorder=3)
        
        conc_range2 = np.logspace(np.log10(df2_processed['concentration'].min()),
                                  np.log10(df2_processed['concentration'].max()), 100)
        log_conc_range2 = np.log10(conc_range2)
        X_pred2 = sm.add_constant(log_conc_range2)
        pred_prob2 = result2.predict(X_pred2)
        pred_mortality_pct2 = pred_prob2 * 100
        
        ax.plot(conc_range2, pred_mortality_pct2,
               color='darkblue', linewidth=2, label=f'{strain2} (Fitted)', zorder=2)
        
        # LD50 reference lines
        ld50_1, _, _ = compute_ldx_with_ci(result1, 50)
        ld50_2, _, _ = compute_ldx_with_ci(result2, 50)
        ax.axhline(y=50, color='gray', linestyle='--', alpha=0.5, linewidth=1, zorder=1)
        ax.axvline(x=ld50_1, color='red', linestyle=':', alpha=0.5, linewidth=1, zorder=1)
        ax.axvline(x=ld50_2, color='blue', linestyle=':', alpha=0.5, linewidth=1, zorder=1)
        
        ax.set_xlabel(f'{chemical1} Concentration', fontsize=12, fontweight='bold')
        ax.set_ylabel('Mortality (%)', fontsize=12, fontweight='bold')
        ax.set_title(f'Comparison: {strain1} vs {strain2}\nResistance Ratio: {resistance_ratio:.2f}x',
                    fontsize=14, fontweight='bold')
        ax.set_xscale('log')
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=10, loc='best')
        ax.set_ylim(-5, 105)
        
        plt.tight_layout()
        fig.savefig(tmpfile.name, format='png', dpi=150, bbox_inches='tight')
        plt.close(fig)
        pdf.image(tmpfile.name, x=10, w=190)
        os.unlink(tmpfile.name)
    
    pdf.ln(5)
    
    # Combined probit plot
    pdf.add_page()
    pdf.set_font('Arial', 'B', 14)
    pdf.cell(0, 10, 'Combined Probit Regression Lines', 0, 1)
    pdf.ln(2)
    
    with tempfile.NamedTemporaryFile(delete=False, suffix='.png') as tmpfile:
        fig = create_combined_probit_plot(df1_processed, result1, strain1,
                                         df2_processed, result2, strain2, chemical1)
        fig.savefig(tmpfile.name, format='png', dpi=150, bbox_inches='tight')
        plt.close(fig)
        pdf.image(tmpfile.name, x=10, w=190)
        os.unlink(tmpfile.name)
    
    return pdf

def fig_to_base64(fig):
    """Convert matplotlib figure to base64 string"""
    buf = BytesIO()
    fig.savefig(buf, format='png', dpi=150, bbox_inches='tight')
    buf.seek(0)
    img_str = base64.b64encode(buf.read()).decode()
    buf.close()
    return img_str

# Main App
def main():
    # Header
    st.markdown('<div class="main-header">📊 Probit Analysis Tool</div>', unsafe_allow_html=True)
    st.markdown("**Web Version 10.3** - Professional bioassay analysis in your browser")
    
    # Sidebar
    with st.sidebar:
        st.image("https://via.placeholder.com/300x100/1f77b4/ffffff?text=Probit+Analysis", 
                use_container_width=True)
        st.markdown("### 🎯 Quick Start")
        st.markdown("""
        1. Upload your data file(s)
        2. Review data summary
        3. Run analysis
        4. Download results
        
        **No installation required!**
        """)
        
        st.markdown("---")
        st.markdown("### 📋 Data Format")
        st.markdown("""
        **Tab-delimited text file:**
        ```
        strain_name
        chemical_name
        concentration	n	mortality
        0.250	100	95
        0.125	80	60
        ```
        """)
        
        st.markdown("---")
        st.markdown("### ℹ️ About")
        st.markdown("**Version 10.3 Web**")
        st.markdown("Probit regression analysis for bioassay data")
        st.markdown("USDA ARS Cattle Fever Tick Research Unit")
        st.markdown("Edinburg, TX, USA")
    
    # Main content
    tabs = st.tabs(["📁 Upload Data", "📊 Single Analysis", "⚖️ Compare Two Datasets", "📖 Help"])
    
    # Tab 1: Upload Data
    with tabs[0]:
        st.markdown("## Upload Your Data")
        
        col1, col2 = st.columns(2)
        
        with col1:
            st.markdown("### 📤 Dataset 1")
            file1 = st.file_uploader("Upload first dataset", type=['txt', 'tsv', 'csv'], key='file1')
            
            if file1:
                try:
                    # Read header lines
                    content = file1.getvalue().decode('utf-8')
                    lines = content.split('\n')
                    strain1 = lines[0].strip()
                    chemical1 = lines[1].strip()
                    
                    # IMPROVEMENT: Validate strain and chemical names
                    if not strain1:
                        st.warning("⚠️ Strain name is empty - using 'Unknown Strain'")
                        strain1 = "Unknown Strain"
                    if not chemical1:
                        st.warning("⚠️ Chemical name is empty - using 'Unknown Chemical'")
                        chemical1 = "Unknown Chemical"
                    
                    # Read data
                    df1 = pd.read_csv(StringIO(content), sep='\t', skiprows=2)
                    
                    st.success(f"✓ Loaded: {strain1} - {chemical1}")
                    st.info(f"📊 {len(df1)} observations")
                    
                    # Validate
                    errors, warnings = validate_bioassay_data(df1, strain1, chemical1)
                    
                    if errors:
                        st.error("❌ **Validation Errors:**")
                        for error in errors:
                            st.error(f"  • {error}")
                    else:
                        st.success("✓ Data validation passed")
                    
                    if warnings:
                        st.warning("⚠️ **Warnings:**")
                        for warning in warnings:
                            st.warning(f"  • {warning}")
                    
                    # Preview data
                    with st.expander("👁️ Preview Data"):
                        df1_preview = df1.copy()
                        df1_preview['mortality_%'] = (df1_preview['mortality'] / df1_preview['n'] * 100).round(2)
                        st.dataframe(df1_preview, use_container_width=True)
                    
                    # Store in session state
                    st.session_state.df1 = df1
                    st.session_state.strain1 = strain1
                    st.session_state.chemical1 = chemical1
                    st.session_state.errors1 = errors
                    
                except Exception as e:
                    st.error(f"❌ Error reading file: {str(e)}")
        
        with col2:
            st.markdown("### 📤 Dataset 2 (Optional)")
            file2 = st.file_uploader("Upload second dataset for comparison", type=['txt', 'tsv', 'csv'], key='file2')
            
            if file2:
                try:
                    # Read header lines
                    content = file2.getvalue().decode('utf-8')
                    lines = content.split('\n')
                    strain2 = lines[0].strip()
                    chemical2 = lines[1].strip()
                    
                    # IMPROVEMENT: Validate strain and chemical names
                    if not strain2:
                        st.warning("⚠️ Strain name is empty - using 'Unknown Strain'")
                        strain2 = "Unknown Strain"
                    if not chemical2:
                        st.warning("⚠️ Chemical name is empty - using 'Unknown Chemical'")
                        chemical2 = "Unknown Chemical"
                    
                    # Read data
                    df2 = pd.read_csv(StringIO(content), sep='\t', skiprows=2)
                    
                    st.success(f"✓ Loaded: {strain2} - {chemical2}")
                    st.info(f"📊 {len(df2)} observations")
                    
                    # Validate
                    errors, warnings = validate_bioassay_data(df2, strain2, chemical2)
                    
                    if errors:
                        st.error("❌ **Validation Errors:**")
                        for error in errors:
                            st.error(f"  • {error}")
                    else:
                        st.success("✓ Data validation passed")
                    
                    if warnings:
                        st.warning("⚠️ **Warnings:**")
                        for warning in warnings:
                            st.warning(f"  • {warning}")
                    
                    # Preview data
                    with st.expander("👁️ Preview Data"):
                        df2_preview = df2.copy()
                        df2_preview['mortality_%'] = (df2_preview['mortality'] / df2_preview['n'] * 100).round(2)
                        st.dataframe(df2_preview, use_container_width=True)
                    
                    # Store in session state
                    st.session_state.df2 = df2
                    st.session_state.strain2 = strain2
                    st.session_state.chemical2 = chemical2
                    st.session_state.errors2 = errors
                    
                    # Ask about reference type
                    st.markdown("---")
                    reference_type = st.selectbox(
                        "📋 What does Dataset 2 represent?",
                        [
                            "Susceptible reference strain (standard resistance test)",
                            "Previous generation/timepoint (monitoring resistance changes)",
                            "Different treatment/condition (experimental comparison)",
                            "Another strain (no assumption about susceptibility)"
                        ],
                        key="reference_type",
                        help="This determines how resistance ratios are interpreted"
                    )
                    # No need to manually assign - Streamlit does this automatically with key="reference_type"
                    
                    if reference_type.startswith("Susceptible"):
                        st.info("💡 Results will include resistance classification (RR ≥ 2.0 = resistant)")
                    else:
                        st.info("💡 Results will show relative comparison without resistance classification")
                    
                except Exception as e:
                    st.error(f"❌ Error reading file: {str(e)}")
    
    # Tab 2: Single Analysis
    with tabs[1]:
        st.markdown("## 📊 Single Dataset Analysis")
        
        if 'df1' not in st.session_state:
            st.info("👆 Please upload data in the 'Upload Data' tab first")
        elif st.session_state.get('errors1'):
            st.error("❌ Please fix data validation errors before running analysis")
        else:
            df = st.session_state.df1
            strain = st.session_state.strain1
            chemical = st.session_state.chemical1
            
            if st.button("🚀 Run Analysis", type="primary", use_container_width=True):
                with st.spinner("Analyzing data..."):
                    try:
                        # Preprocess
                        df_processed = preprocess_data(df)
                        
                        if len(df_processed) == 0:
                            st.error("❌ No valid data remaining after preprocessing. Check your data.")
                            st.stop()
                        
                        # Fit model
                        result, chi_square, df_degrees, p_value = fit_probit_model(df_processed)
                        
                    except (ValueError, RuntimeError) as e:
                        st.error("❌ Analysis failed:")
                        st.error(str(e))
                        st.error("\n**Troubleshooting Tips:**")
                        st.error("• Check that concentrations span mortality range from ~10% to ~90%")
                        st.error("• Ensure you have at least 3-4 different concentrations")
                        st.error("• Verify mortality counts don't exceed sample sizes")
                        st.error("• Consider adding more replicates if variability is high")
                        st.stop()
                    except Exception as e:
                        st.error(f"❌ Unexpected error during analysis: {str(e)}")
                        st.exception(e)
                        st.stop()
                    
                    try:
                        
                        # Data Summary
                        st.markdown("### 📋 Data Summary")
                        col1, col2, col3, col4 = st.columns(4)
                        with col1:
                            st.metric("Observations", len(df))
                        with col2:
                            st.metric("Total Tested", int(df['n'].sum()))
                        with col3:
                            st.metric("Total Deaths", int(df['mortality'].sum()))
                        with col4:
                            overall_mort = (df['mortality'].sum() / df['n'].sum() * 100)
                            st.metric("Overall Mortality", f"{overall_mort:.1f}%")
                        
                        # LD Estimates
                        st.markdown("### 💊 Lethal Dose Estimates")
                        try:
                            ld_data = []
                            for ld_level in LD_LEVELS:
                                ld, lower, upper = compute_ldx_with_ci(result, ld_level)
                                ld_data.append({
                                    'LD Level': f'LD{ld_level}',
                                    'Estimate': f'{ld:.6f}',
                                    '95% CI Lower': f'{lower:.6f}',
                                    '95% CI Upper': f'{upper:.6f}'
                                })
                            st.table(pd.DataFrame(ld_data))
                        except ValueError as e:
                            st.error("❌ Error calculating LD estimates:")
                            st.error(str(e))
                            st.error("Analysis cannot continue with invalid LD estimates.")
                            st.stop()
                        
                        # Model Fit
                        st.markdown("### 📈 Model Fit")
                        col1, col2 = st.columns(2)
                        with col1:
                            st.metric("Chi-Square", f"{chi_square:.4f}")
                            st.metric("Degrees of Freedom", df_degrees)
                        with col2:
                            st.metric("p-value", f"{p_value:.4f}")
                            if p_value < ALPHA_LEVEL:
                                st.error("❌ Model does NOT fit well (p < 0.05)")
                            else:
                                st.success("✓ Model fits well (p > 0.05)")
                        
                        # Model Parameters
                        st.markdown("### 🔢 Model Parameters")
                        params_data = []
                        for param in result.params.index:
                            params_data.append({
                                'Parameter': param,
                                'Estimate': f'{result.params[param]:.4f}',
                                'Std Error': f'{result.bse[param]:.4f}',
                                'p-value': f'{result.pvalues[param]:.4f}'
                            })
                        st.table(pd.DataFrame(params_data))
                        
                        # Replicate Variability
                        st.markdown("### 📊 Replicate Variability")
                        var_df = calculate_replicate_variability(df)
                        var_display = var_df[['concentration', 'mean_mortality_pct', 'std_mortality_pct', 
                                              'cv_pct', 'n_replicates']].copy()
                        var_display.columns = ['Concentration', 'Mean %', 'Std Dev %', 'CV%', 'Replicates']
                        var_display['Mean %'] = var_display['Mean %'].round(2)
                        var_display['Std Dev %'] = var_display['Std Dev %'].round(2)
                        var_display['CV%'] = var_display['CV%'].round(1)
                        
                        # Highlight high CV
                        def highlight_high_cv(row):
                            if row['CV%'] > 20:
                                return ['background-color: #ffcccc'] * len(row)
                            return [''] * len(row)
                        
                        st.dataframe(var_display.style.apply(highlight_high_cv, axis=1), 
                                   use_container_width=True)
                        st.caption("⚠️ Red = CV% > 20% (high variability)")
                        
                        # Plot
                        st.markdown("### 📉 Mortality Curve")
                        fig = create_mortality_plot(df_processed, result, strain, chemical)
                        st.pyplot(fig)
                        
                        # Probit plot
                        st.markdown("### 📈 Probit Regression Line")
                        fig_probit = create_probit_plot(df_processed, result, strain, chemical)
                        st.pyplot(fig_probit)
                        
                        # Generate PDF Report
                        st.markdown("### 📄 Download Report")
                        
                        # PDF customization options
                        with st.expander("⚙️ Customize PDF Report", expanded=False):
                            st.markdown("**Select what to include in the PDF:**")
                            
                            col1, col2 = st.columns(2)
                            with col1:
                                include_data_summary = st.checkbox("Data Summary", value=True, 
                                    help="Observations, total tested, mortality")
                                include_raw_data = st.checkbox("Raw Data Table", value=True,
                                    help="Table of all concentrations and mortality")
                                include_ld_estimates = st.checkbox("LD Estimates", value=True,
                                    help="LD1, LD50, LD99 with confidence intervals")
                                include_model_fit = st.checkbox("Model Fit Statistics", value=True,
                                    help="Chi-square, p-value, goodness of fit")
                            
                            with col2:
                                include_parameters = st.checkbox("Model Parameters", value=True,
                                    help="Intercept and slope with standard errors")
                                include_mortality_plot = st.checkbox("Mortality Curve Plot", value=True,
                                    help="Observed vs. fitted mortality curve")
                                include_probit_plot = st.checkbox("Probit Regression Plot", value=True,
                                    help="Probit transformation with fitted line")
                        
                        # Create customized PDF
                        pdf_options = {
                            'include_data_summary': include_data_summary,
                            'include_raw_data': include_raw_data,
                            'include_ld_estimates': include_ld_estimates,
                            'include_model_fit': include_model_fit,
                            'include_parameters': include_parameters,
                            'include_mortality_plot': include_mortality_plot,
                            'include_probit_plot': include_probit_plot
                        }
                        
                        pdf = create_pdf_single_with_plots(df, df_processed, result, chi_square, df_degrees, p_value, 
                                                          strain, chemical, options=pdf_options)
                        
                        # Save PDF to bytes
                        pdf_output = BytesIO()
                        pdf_bytes = pdf.output(dest='S').encode('latin-1')
                        pdf_output.write(pdf_bytes)
                        pdf_output.seek(0)
                        
                        # Download button
                        st.download_button(
                            label="📥 Download PDF Report (with plots)",
                            data=pdf_output,
                            file_name=f"probit_analysis_{strain}_{datetime.now().strftime('%Y%m%d_%H%M%S')}.pdf",
                            mime="application/pdf",
                            type="primary"
                        )
                        
                        # Success message
                        st.success("✓ Analysis complete!")
                        
                    except Exception as e:
                        st.error(f"❌ Error during analysis: {str(e)}")
                        st.exception(e)
    
    # Tab 3: Compare Two Datasets
    with tabs[2]:
        st.markdown("## ⚖️ Compare Two Datasets")
        
        if 'df1' not in st.session_state or 'df2' not in st.session_state:
            st.info("👆 Please upload two datasets in the 'Upload Data' tab first")
        elif st.session_state.get('errors1') or st.session_state.get('errors2'):
            st.error("❌ Please fix data validation errors before running comparison")
        else:
            df1 = st.session_state.df1
            strain1 = st.session_state.strain1
            chemical1 = st.session_state.chemical1
            
            df2 = st.session_state.df2
            strain2 = st.session_state.strain2
            chemical2 = st.session_state.chemical2
            
            # Get reference type (default to susceptible if not set)
            reference_type = st.session_state.get('reference_type', 
                                                  'Susceptible reference strain (standard resistance test)')
            is_susceptible_reference = reference_type.startswith("Susceptible")
            
            if st.button("🚀 Run Comparison", type="primary", use_container_width=True):
                with st.spinner("Comparing datasets..."):
                    try:
                        # Process both datasets
                        df1_processed = preprocess_data(df1)
                        df2_processed = preprocess_data(df2)
                        
                        if len(df1_processed) == 0 or len(df2_processed) == 0:
                            st.error("❌ No valid data remaining after preprocessing. Check your data.")
                            st.stop()
                        
                        # Fit models with error handling
                        try:
                            result1, chi_sq1, df_deg1, p_val1 = fit_probit_model(df1_processed)
                        except (ValueError, RuntimeError) as e:
                            st.error(f"❌ Model fitting failed for {strain1}:")
                            st.error(str(e))
                            st.stop()
                        
                        try:
                            result2, chi_sq2, df_deg2, p_val2 = fit_probit_model(df2_processed)
                        except (ValueError, RuntimeError) as e:
                            st.error(f"❌ Model fitting failed for {strain2}:")
                            st.error(str(e))
                            st.stop()
                        
                    except Exception as e:
                        st.error(f"❌ Unexpected error during model fitting: {str(e)}")
                        st.exception(e)
                        st.stop()
                    
                    try:
                        
                        # Side-by-side data summary
                        st.markdown("### 📋 Data Summary Comparison")
                        col1, col2 = st.columns(2)
                        
                        with col1:
                            st.markdown(f"**{strain1} - {chemical1}**")
                            st.metric("Observations", len(df1))
                            st.metric("Total Tested", int(df1['n'].sum()))
                            st.metric("Total Deaths", int(df1['mortality'].sum()))
                            overall_mort1 = (df1['mortality'].sum() / df1['n'].sum() * 100)
                            st.metric("Overall Mortality", f"{overall_mort1:.1f}%")
                        
                        with col2:
                            st.markdown(f"**{strain2} - {chemical2}**")
                            st.metric("Observations", len(df2))
                            st.metric("Total Tested", int(df2['n'].sum()))
                            st.metric("Total Deaths", int(df2['mortality'].sum()))
                            overall_mort2 = (df2['mortality'].sum() / df2['n'].sum() * 100)
                            st.metric("Overall Mortality", f"{overall_mort2:.1f}%")
                        
                        # Compute LD50 for both with error handling
                        try:
                            ld50_1, ld50_1_lower, ld50_1_upper = compute_ldx_with_ci(result1, 50)
                            ld50_2, ld50_2_lower, ld50_2_upper = compute_ldx_with_ci(result2, 50)
                        except ValueError as e:
                            st.error(f"❌ Error computing LD50 values: {str(e)}")
                            st.error("Cannot proceed with resistance ratio calculation.")
                            st.stop()
                        
                        # CRITICAL FIX: Validate LD50 values before resistance ratio calculation
                        if not (np.isfinite(ld50_1) and np.isfinite(ld50_2) and ld50_1 > 0 and ld50_2 > 0):
                            st.error("❌ Invalid LD50 values detected:")
                            st.error(f"  • {strain1} LD50: {ld50_1}")
                            st.error(f"  • {strain2} LD50: {ld50_2}")
                            st.error("This may indicate:")
                            st.error("  • Model fitting problems")
                            st.error("  • Extreme resistance levels")
                            st.error("  • Poor data quality")
                            st.error("Cannot calculate resistance ratio.")
                            st.stop()
                        
                        # CRITICAL FIX: Check for division by zero
                        if ld50_2 == 0:
                            st.error("❌ Reference strain LD50 is zero - cannot calculate resistance ratio.")
                            st.stop()
                        
                        # Resistance ratio calculation
                        resistance_ratio = ld50_1 / ld50_2
                        
                        # CRITICAL FIX: Validate resistance ratio
                        if not np.isfinite(resistance_ratio) or resistance_ratio <= 0:
                            st.error("❌ Invalid resistance ratio calculated:")
                            st.error(f"  • Resistance ratio: {resistance_ratio}")
                            st.error(f"  • LD50 ratio: {ld50_1} / {ld50_2}")
                            st.error("Check LD50 estimates and model parameters.")
                            st.stop()
                        
                        # Resistance ratio confidence interval (Fieller's method approximation)
                        try:
                            log_rr = np.log(resistance_ratio)
                            se_log_rr = np.sqrt(
                                (result1.bse[0]**2 + result1.bse[1]**2 * (norm.ppf(0.5) - result1.params[0])**2 / result1.params[1]**2) / ld50_1**2 +
                                (result2.bse[0]**2 + result2.bse[1]**2 * (norm.ppf(0.5) - result2.params[0])**2 / result2.params[1]**2) / ld50_2**2
                            )
                            
                            # CRITICAL FIX: Validate confidence interval calculation
                            if not np.isfinite(se_log_rr) or se_log_rr <= 0:
                                st.warning("⚠️ Could not calculate reliable confidence intervals for resistance ratio.")
                                st.warning("Using point estimate only.")
                                rr_lower = resistance_ratio * 0.5  # Rough approximation
                                rr_upper = resistance_ratio * 2.0  
                            else:
                                rr_lower = np.exp(log_rr - 1.96 * se_log_rr)
                                rr_upper = np.exp(log_rr + 1.96 * se_log_rr)
                                
                                # Validate bounds
                                if not (np.isfinite(rr_lower) and np.isfinite(rr_upper) and rr_lower > 0 and rr_upper > 0):
                                    st.warning("⚠️ Confidence interval calculation produced invalid bounds.")
                                    st.warning("Using point estimate only.")
                                    rr_lower = resistance_ratio * 0.5
                                    rr_upper = resistance_ratio * 2.0
                                    
                        except Exception as e:
                            st.warning(f"⚠️ Error in confidence interval calculation: {str(e)}")
                            st.warning("Using point estimate only.")
                            rr_lower = resistance_ratio * 0.5
                            rr_upper = resistance_ratio * 2.0
                        
                        # LD Estimates Comparison
                        st.markdown("### 💊 Lethal Dose Comparison")
                        
                        comparison_data = []
                        for ld_level in LD_LEVELS:
                            ld1, lower1, upper1 = compute_ldx_with_ci(result1, ld_level)
                            ld2, lower2, upper2 = compute_ldx_with_ci(result2, ld_level)
                            ratio = ld1 / ld2
                            
                            comparison_data.append({
                                'LD Level': f'LD{ld_level}',
                                f'{strain1}': f'{ld1:.6f}',
                                f'{strain1} 95% CI': f'({lower1:.6f} - {upper1:.6f})',
                                f'{strain2}': f'{ld2:.6f}',
                                f'{strain2} 95% CI': f'({lower2:.6f} - {upper2:.6f})',
                                'Ratio': f'{ratio:.2f}'
                            })
                        
                        st.table(pd.DataFrame(comparison_data))
                        
                        # Resistance Ratio
                        st.markdown("### 🔬 Resistance Ratio")
                        st.success(f"""
                        **Resistance Ratio (LD50 basis):** {resistance_ratio:.2f}x
                        
                        **95% Confidence Interval:** ({rr_lower:.2f} - {rr_upper:.2f})
                        
                        **Interpretation:** {strain1} is {resistance_ratio:.2f}x {'more resistant' if resistance_ratio > 1 else 'more susceptible'} 
                        than {strain2} to {chemical1}
                        """)
                        
                        # Statistical Tests
                        st.markdown("### 📊 Statistical Tests")
                        
                        # Test for equality of slopes (parallelism test)
                        slope1 = result1.params[1]
                        slope2 = result2.params[1]
                        se_slope1 = result1.bse[1]
                        se_slope2 = result2.bse[1]
                        
                        z_parallel = (slope1 - slope2) / np.sqrt(se_slope1**2 + se_slope2**2)
                        p_parallel = 2 * (1 - norm.cdf(abs(z_parallel)))
                        
                        col1, col2 = st.columns(2)
                        
                        with col1:
                            st.markdown("**Parallelism Test**")
                            st.caption("*Z-test for equality of slopes (Robertson et al., 2007)*")
                            st.metric("Z-statistic", f"{z_parallel:.4f}")
                            st.metric("p-value", f"{p_parallel:.4f}")
                            if p_parallel > 0.05:
                                st.success("✓ Slopes are parallel (p > 0.05)")
                                st.caption("Lines have similar slopes - good for comparing potency")
                            else:
                                st.warning("⚠️ Slopes differ significantly (p < 0.05)")
                                st.caption("Lines are not parallel - interpret resistance ratio with caution")
                        
                        # Test for equality of intercepts (given parallel slopes)
                        intercept1 = result1.params[0]
                        intercept2 = result2.params[0]
                        se_intercept1 = result1.bse[0]
                        se_intercept2 = result2.bse[0]
                        
                        z_equality = (intercept1 - intercept2) / np.sqrt(se_intercept1**2 + se_intercept2**2)
                        p_equality = 2 * (1 - norm.cdf(abs(z_equality)))
                        
                        with col2:
                            st.markdown("**Equality Test**")
                            st.caption("*Z-test for equality of intercepts*")
                            st.metric("Z-statistic", f"{z_equality:.4f}")
                            st.metric("p-value", f"{p_equality:.4f}")
                            if p_equality > 0.05:
                                st.info("No significant difference (p > 0.05)")
                                st.caption("Strains show similar susceptibility")
                            else:
                                st.success("✓ Significant difference (p < 0.05)")
                                st.caption("Strains have different susceptibility levels")
                        
                        # Model Fit Comparison
                        st.markdown("### 📈 Model Fit Comparison")
                        
                        fit_comparison = pd.DataFrame({
                            'Dataset': [strain1, strain2],
                            'Chi-Square': [f'{chi_sq1:.4f}', f'{chi_sq2:.4f}'],
                            'df': [df_deg1, df_deg2],
                            'p-value': [f'{p_val1:.4f}', f'{p_val2:.4f}'],
                            'Fit Quality': [
                                '✓ Good' if p_val1 >= 0.05 else '❌ Poor',
                                '✓ Good' if p_val2 >= 0.05 else '❌ Poor'
                            ]
                        })
                        st.table(fit_comparison)
                        
                        # Replicate Variability Comparison
                        st.markdown("### 📊 Replicate Variability Comparison")
                        
                        col1, col2 = st.columns(2)
                        
                        with col1:
                            st.markdown(f"**{strain1}**")
                            var_df1 = calculate_replicate_variability(df1)
                            var_display1 = var_df1[['concentration', 'mean_mortality_pct', 'cv_pct', 'n_replicates']].copy()
                            var_display1.columns = ['Conc.', 'Mean %', 'CV%', 'Reps']
                            var_display1['Mean %'] = var_display1['Mean %'].round(1)
                            var_display1['CV%'] = var_display1['CV%'].round(1)
                            
                            def highlight_cv(row):
                                if row['CV%'] > 20:
                                    return ['background-color: #ffcccc'] * len(row)
                                return [''] * len(row)
                            
                            st.dataframe(var_display1.style.apply(highlight_cv, axis=1), 
                                       use_container_width=True)
                        
                        with col2:
                            st.markdown(f"**{strain2}**")
                            var_df2 = calculate_replicate_variability(df2)
                            var_display2 = var_df2[['concentration', 'mean_mortality_pct', 'cv_pct', 'n_replicates']].copy()
                            var_display2.columns = ['Conc.', 'Mean %', 'CV%', 'Reps']
                            var_display2['Mean %'] = var_display2['Mean %'].round(1)
                            var_display2['CV%'] = var_display2['CV%'].round(1)
                            
                            st.dataframe(var_display2.style.apply(highlight_cv, axis=1), 
                                       use_container_width=True)
                        
                        # Combined Plot
                        st.markdown("### 📉 Combined Mortality Curves")
                        
                        fig, ax = plt.subplots(figsize=(12, 7))
                        
                        # Dataset 1
                        mortality_pct1 = (df1_processed['mortality'] / df1_processed['n']) * 100
                        ax.scatter(df1_processed['concentration'], mortality_pct1,
                                  color='red', s=100, alpha=0.6, label=f'{strain1} (Observed)', zorder=3)
                        
                        conc_range1 = np.logspace(np.log10(df1_processed['concentration'].min()),
                                                  np.log10(df1_processed['concentration'].max()), 100)
                        log_conc_range1 = np.log10(conc_range1)
                        X_pred1 = sm.add_constant(log_conc_range1)
                        pred_prob1 = result1.predict(X_pred1)
                        pred_mortality_pct1 = pred_prob1 * 100
                        
                        ax.plot(conc_range1, pred_mortality_pct1,
                               color='darkred', linewidth=2, label=f'{strain1} (Fitted)', zorder=2)
                        
                        # Dataset 2
                        mortality_pct2 = (df2_processed['mortality'] / df2_processed['n']) * 100
                        ax.scatter(df2_processed['concentration'], mortality_pct2,
                                  color='blue', s=100, alpha=0.6, label=f'{strain2} (Observed)', zorder=3)
                        
                        conc_range2 = np.logspace(np.log10(df2_processed['concentration'].min()),
                                                  np.log10(df2_processed['concentration'].max()), 100)
                        log_conc_range2 = np.log10(conc_range2)
                        X_pred2 = sm.add_constant(log_conc_range2)
                        pred_prob2 = result2.predict(X_pred2)
                        pred_mortality_pct2 = pred_prob2 * 100
                        
                        ax.plot(conc_range2, pred_mortality_pct2,
                               color='darkblue', linewidth=2, label=f'{strain2} (Fitted)', zorder=2)
                        
                        # LD50 reference lines
                        ax.axhline(y=50, color='gray', linestyle='--', alpha=0.5, linewidth=1, zorder=1)
                        ax.axvline(x=ld50_1, color='red', linestyle=':', alpha=0.5, linewidth=1, zorder=1)
                        ax.axvline(x=ld50_2, color='blue', linestyle=':', alpha=0.5, linewidth=1, zorder=1)
                        
                        ax.set_xlabel(f'{chemical1} Concentration', fontsize=12, fontweight='bold')
                        ax.set_ylabel('Mortality (%)', fontsize=12, fontweight='bold')
                        ax.set_title(f'Comparison: {strain1} vs {strain2}\nResistance Ratio: {resistance_ratio:.2f}x',
                                    fontsize=14, fontweight='bold')
                        ax.set_xscale('log')
                        ax.grid(True, alpha=0.3)
                        ax.legend(fontsize=10, loc='best')
                        ax.set_ylim(-5, 105)
                        
                        plt.tight_layout()
                        st.pyplot(fig)
                        
                        # Combined Probit Plot
                        st.markdown("### 📈 Combined Probit Regression Lines")
                        fig_probit = create_combined_probit_plot(df1_processed, result1, strain1,
                                                                 df2_processed, result2, strain2, chemical1)
                        st.pyplot(fig_probit)
                        
                        # Summary interpretation
                        st.markdown("### 📝 Summary Interpretation")
                        
                        # Determine resistance status based on reference type
                        ci_includes_one = (rr_lower <= 1.0 <= rr_upper)
                        biologically_significant = (resistance_ratio >= 2.0 or resistance_ratio <= 0.5)
                        
                        # Check for extremely wide CI (indicates unreliable estimate)
                        ci_width = rr_upper - rr_lower
                        ci_is_unreliable = (ci_width > 1000) or (rr_upper > 1000000) or (rr_lower < 0.001)
                        
                        # Check for extreme RR values
                        extreme_resistance = resistance_ratio >= 10.0
                        extreme_susceptibility = resistance_ratio <= 0.1
                        
                        if is_susceptible_reference:
                            # Standard resistance classification with improved logic
                            if extreme_resistance:
                                # Very high RR overrides CI uncertainty
                                resistance_status = f"{resistance_ratio:.2f}-fold resistance"
                                resistance_classification = "highly resistant"
                                ci_note = f'- **Note:** CI is very wide ({rr_lower:.2f} - {rr_upper:.2f}), indicating high uncertainty in the exact value, but resistance level is clearly high' if ci_is_unreliable else ''
                            elif extreme_susceptibility:
                                # Very low RR overrides CI uncertainty  
                                resistance_status = f"{resistance_ratio:.2f}-fold increased susceptibility"
                                resistance_classification = "highly susceptible"
                                ci_note = f'- **Note:** CI is very wide ({rr_lower:.2f} - {rr_upper:.2f}), indicating high uncertainty in the exact value, but susceptibility is clearly high' if ci_is_unreliable else ''
                            elif ci_includes_one and not ci_is_unreliable:
                                # Normal CI that includes 1.0 and is reliable
                                resistance_status = f"no significant difference from {strain2}"
                                resistance_classification = "similar susceptibility"
                                ci_note = f'- **Note:** CI includes 1.0, indicating no statistically significant difference'
                            elif ci_includes_one and ci_is_unreliable:
                                # Unreliable CI that includes 1.0 - use RR magnitude
                                if biologically_significant:
                                    if resistance_ratio > 1:
                                        resistance_status = f"{resistance_ratio:.2f}-fold resistance (CI unreliable)"
                                        resistance_classification = "resistant"
                                    else:
                                        resistance_status = f"{resistance_ratio:.2f}-fold increased susceptibility (CI unreliable)"
                                        resistance_classification = "more susceptible"
                                    ci_note = f'- **Note:** CI is extremely wide ({rr_lower:.2f} - {rr_upper:.2f}) and unreliable. Interpretation based on RR point estimate.'
                                else:
                                    resistance_status = f"uncertain - CI too wide for reliable interpretation"
                                    resistance_classification = "uncertain"
                                    ci_note = f'- **Note:** CI is extremely wide ({rr_lower:.2f} - {rr_upper:.2f}). Consider additional replicates for reliable estimate.'
                            elif not biologically_significant:
                                # Statistically significant but small effect
                                if resistance_ratio > 1:
                                    resistance_status = f"statistically different but low-level resistance ({resistance_ratio:.2f}-fold)"
                                    resistance_classification = "low-level resistance (RR < 2)"
                                else:
                                    resistance_status = f"statistically different but similar susceptibility ({resistance_ratio:.2f}-fold)"
                                    resistance_classification = "similar to susceptible"
                                ci_note = f'- **Note:** RR < 2, indicating low-level or no practical resistance'
                            else:
                                # Clear biological significance with reliable CI
                                if resistance_ratio > 1:
                                    resistance_status = f"{resistance_ratio:.2f}-fold resistance"
                                    resistance_classification = "resistant"
                                else:
                                    resistance_status = f"{resistance_ratio:.2f}-fold increased susceptibility"
                                    resistance_classification = "more susceptible"
                                ci_note = ''
                            
                            interpretation = f"""
                            **Comparison Type:** Resistance test against susceptible reference
                            
                            **Key Findings:**
                            
                            1. **Resistance Level:** {strain1} shows {resistance_status} 
                               compared to {strain2} (RR = {resistance_ratio:.2f}, 95% CI: {rr_lower:.2f} - {rr_upper:.2f})
                               {ci_note}
                            
                            2. **Statistical Significance:** The difference between strains is 
                               {'statistically significant (p < 0.05)' if p_equality < 0.05 else 'not statistically significant (p > 0.05)'}
                            
                            3. **Dose-Response Relationship:** The regression lines are 
                               {'parallel (p > 0.05)' if p_parallel > 0.05 else 'not parallel (p < 0.05)'}, 
                               {'indicating similar modes of action' if p_parallel > 0.05 else 'suggesting different mechanisms may be involved (Robertson et al., 2007)'}
                            
                            4. **Model Quality:** 
                               - {strain1}: {'Good fit (p > 0.05)' if p_val1 >= 0.05 else f'Poor fit (p = {p_val1:.4f}) - check variability'}
                               - {strain2}: {'Good fit (p > 0.05)' if p_val2 >= 0.05 else f'Poor fit (p = {p_val2:.4f}) - check variability'}
                            
                            **For Publication:**
                            
                            {f'{strain1} showed no significant difference in susceptibility to {chemical1} compared to {strain2} (LD50: {ld50_1:.3f} vs {ld50_2:.3f}, RR = {resistance_ratio:.2f}, 95% CI: {rr_lower:.2f}-{rr_upper:.2f}).' if (ci_includes_one and not ci_is_unreliable and not extreme_resistance and not extreme_susceptibility) else
                             f'{strain1} exhibited {resistance_ratio:.2f}-fold resistance to {chemical1} compared to {strain2} (LD50: {ld50_1:.3f} vs {ld50_2:.3f}, RR = {resistance_ratio:.2f}, 95% CI: {rr_lower:.2f}-{rr_upper:.2f}){", though confidence intervals are very wide indicating high uncertainty" if ci_is_unreliable else ""}.' if resistance_ratio > 1 else
                             f'{strain1} exhibited {resistance_ratio:.2f}-fold increased susceptibility to {chemical1} compared to {strain2} (LD50: {ld50_1:.3f} vs {ld50_2:.3f}, RR = {resistance_ratio:.2f}, 95% CI: {rr_lower:.2f}-{rr_upper:.2f}){", though confidence intervals are very wide indicating high uncertainty" if ci_is_unreliable else ""}.'}
                            {'The dose-response curves were parallel (p > 0.05), indicating similar modes of action.' if p_parallel > 0.05 else 
                            'The dose-response curves showed significantly different slopes (p < 0.05), suggesting potential mechanistic differences.'}"
                            """
                        else:
                            # Relative comparison without resistance classification
                            comparison_direction = "higher" if resistance_ratio > 1 else "lower"
                            percent_change = abs((resistance_ratio - 1.0) * 100)
                            
                            interpretation = f"""
                            **Comparison Type:** {reference_type}
                            
                            **Key Findings:**
                            
                            1. **Relative Difference:** {strain1} has {resistance_ratio:.2f}-fold {comparison_direction} LD50 than {strain2}
                               (LD50: {ld50_1:.3f} vs {ld50_2:.3f}, RR = {resistance_ratio:.2f}, 95% CI: {rr_lower:.2f} - {rr_upper:.2f})
                               - This represents a {percent_change:.1f}% {'increase' if resistance_ratio > 1 else 'decrease'} in LD50
                               {f'- **Note:** CI includes 1.0, indicating no statistically significant difference' if ci_includes_one else ''}
                               {f'- **Note:** CI excludes 1.0, indicating a statistically significant difference' if not ci_includes_one else ''}
                            
                            2. **Statistical Significance:** The difference between datasets is 
                               {'statistically significant (p < 0.05)' if p_equality < 0.05 else 'not statistically significant (p > 0.05)'}
                            
                            3. **Dose-Response Relationship:** The regression lines are 
                               {'parallel (p > 0.05)' if p_parallel > 0.05 else 'not parallel (p < 0.05)'}, 
                               {'indicating similar response patterns' if p_parallel > 0.05 else 'suggesting different response patterns (Robertson et al., 2007)'}
                            
                            4. **Model Quality:** 
                               - {strain1}: {'Good fit (p > 0.05)' if p_val1 >= 0.05 else f'Poor fit (p = {p_val1:.4f}) - check variability'}
                               - {strain2}: {'Good fit (p > 0.05)' if p_val2 >= 0.05 else f'Poor fit (p = {p_val2:.4f}) - check variability'}
                            
                            **For Publication:**
                            
                            {f'{strain1} showed no significant difference in susceptibility to {chemical1} compared to {strain2} (LD50: {ld50_1:.3f} vs {ld50_2:.3f}, RR = {resistance_ratio:.2f}, 95% CI: {rr_lower:.2f}-{rr_upper:.2f}).' if ci_includes_one else
                             f'{strain1} had a {resistance_ratio:.2f}-fold {comparison_direction} LD50 for {chemical1} compared to {strain2} (LD50: {ld50_1:.3f} vs {ld50_2:.3f}, RR = {resistance_ratio:.2f}, 95% CI: {rr_lower:.2f}-{rr_upper:.2f}), representing a {percent_change:.1f}% {"increase" if resistance_ratio > 1 else "decrease"} in tolerance.'}
                            {'The dose-response curves were parallel (p > 0.05), indicating similar response patterns.' if p_parallel > 0.05 else 
                            'The dose-response curves showed significantly different slopes (p < 0.05), suggesting different response patterns.'}"
                            """
                        
                        st.info(interpretation)
                        
                        # Generate PDF Report
                        st.markdown("### 📄 Download Report")
                        
                        # PDF customization options for comparison
                        with st.expander("⚙️ Customize PDF Report", expanded=False):
                            st.markdown("**Select what to include in the PDF:**")
                            
                            col1, col2 = st.columns(2)
                            with col1:
                                comp_include_raw_data = st.checkbox("Raw Data Tables", value=True,
                                    help="Tables for both datasets", key="comp_raw_data")
                                comp_include_rr = st.checkbox("Resistance Ratio", value=True,
                                    help="RR with confidence interval", key="comp_rr")
                                comp_include_ld_comparison = st.checkbox("LD Comparison Table", value=True,
                                    help="Side-by-side LD estimates", key="comp_ld")
                                comp_include_statistical_tests = st.checkbox("Statistical Tests", value=True,
                                    help="Parallelism and equality tests", key="comp_stats")
                            
                            with col2:
                                comp_include_model_fit = st.checkbox("Model Fit Comparison", value=True,
                                    help="Chi-square for both datasets", key="comp_fit")
                                comp_include_mortality_plot = st.checkbox("Combined Mortality Plot", value=True,
                                    help="Both curves on one plot", key="comp_mort_plot")
                                comp_include_probit_plot = st.checkbox("Combined Probit Plot", value=True,
                                    help="Both regression lines", key="comp_prob_plot")
                        
                        # Create customized comparison PDF
                        comp_pdf_options = {
                            'include_raw_data': comp_include_raw_data,
                            'include_rr': comp_include_rr,
                            'include_ld_comparison': comp_include_ld_comparison,
                            'include_statistical_tests': comp_include_statistical_tests,
                            'include_model_fit': comp_include_model_fit,
                            'include_mortality_plot': comp_include_mortality_plot,
                            'include_probit_plot': comp_include_probit_plot
                        }
                        
                        pdf = create_pdf_comparison_with_plots(
                            df1, df1_processed, result1, chi_sq1, df_deg1, p_val1, strain1, chemical1,
                            df2, df2_processed, result2, chi_sq2, df_deg2, p_val2, strain2, chemical2,
                            resistance_ratio, rr_lower, rr_upper, p_parallel, p_equality,
                            options=comp_pdf_options
                        )
                        
                        # Save PDF to bytes
                        pdf_output = BytesIO()
                        pdf_bytes = pdf.output(dest='S').encode('latin-1')
                        pdf_output.write(pdf_bytes)
                        pdf_output.seek(0)
                        
                        # Download button
                        st.download_button(
                            label="📥 Download Comparison PDF Report (with plots)",
                            data=pdf_output,
                            file_name=f"probit_comparison_{strain1}_vs_{strain2}_{datetime.now().strftime('%Y%m%d_%H%M%S')}.pdf",
                            mime="application/pdf",
                            type="primary"
                        )
                        
                        st.success("✓ Comparison analysis complete!")
                        
                    except Exception as e:
                        st.error(f"❌ Error during comparison: {str(e)}")
                        st.exception(e)
    
    # Tab 4: Help
    with tabs[3]:
        st.markdown("## 📖 Help & Documentation")
        
        st.markdown("### 🎯 How to Use")
        st.markdown("""
        1. **Upload Data**: Click the 'Upload Data' tab and upload your tab-delimited text file(s)
        2. **Review**: Check that your data loaded correctly and passed validation
        3. **Analyze**: Go to 'Single Analysis' or 'Compare Two Datasets' and click 'Run Analysis'
        4. **Download**: Download your results as PDF
        
        ### 📋 Data Format Requirements
        
        Your file must be a tab-delimited text file with this format:
        ```
        Line 1: Strain name
        Line 2: Chemical name
        Line 3: concentration	n	mortality
        Line 4+: Data rows
        ```
        
        **Example:**
        ```
        Yucatan_Strain
        Coumaphos
        concentration	n	mortality
        0.250	100	95
        0.125	80	60
        0.063	75	20
        ```
        
        ### ✅ What Gets Analyzed
        
        - **LD1, LD50, LD99**: Lethal dose estimates with 95% confidence intervals
        - **Model Fit**: Chi-square goodness-of-fit test
        - **Parameters**: Probit regression coefficients
        - **Variability**: Replicate variability analysis (CV%)
        - **Plots**: Observed vs. fitted mortality curves
        
        ### ⚠️ Common Issues & Solutions
        
        **"Model did not converge"**
        - Insufficient concentration spacing (use 4-6 concentrations spanning 10-90% mortality)
        - Poor concentration selection (too narrow range)
        - High variability between replicates (check CV% > 20%)
        - Try logarithmic concentration spacing (e.g., 0.1, 0.3, 1.0, 3.0, 10.0)
        
        **"Slope is too close to zero"**
        - Dose-response curve is flat (all mortalities similar)
        - Concentration range too narrow
        - Chemical may be ineffective at tested doses
        - Try wider concentration range or different chemical concentrations
        
        **"Invalid LD50 values detected"**
        - Model fitting problems due to extreme data
        - Complete separation (all 0% or 100% mortality)
        - Numerical instability from poor data quality
        - Ensure mortality ranges from ~10% to ~90% across concentrations

        **"Mortality exceeds sample size"**
        - Check that mortality counts are ≤ n for all rows
        - Example: n=96, mortality=98 is impossible!

        **"Model does NOT fit well"**
        - High variability between replicates (p < 0.05)
        - Check CV% in variability table
        - Results are still valid, but confidence intervals may be underestimated
        - Consider additional replicates to reduce variability

        **"High variability (CV% > 20%)"**
        - Some concentrations show inconsistent responses
        - May indicate experimental issues or biological variability
        - Consider additional replicates or review experimental protocol
        
        ### 🌐 Web Version Features
        
        **Advantages:**
        - ✅ No installation required
        - ✅ Works on any device with a browser
        - ✅ No IT approval needed
        - ✅ Always up-to-date
        - ✅ Secure - data never leaves your browser
        
        **Limitations:**
        - PDF generation under development
        - Analysis runs in your browser (slower for large datasets)
        - No local file saving (use download buttons)
        
        ### 🔒 Privacy & Security
        
        - All analysis runs in your browser
        - Your data is NEVER uploaded to any server
        - Results are computed locally
        - Completely private and secure
        
        ### 📊 Statistical Methods
        
        #### **Probit Regression**
        The tool fits a probit regression model using generalized linear models (GLM) with a probit link function:
        
        ```
        Probit(mortality) = β₀ + β₁ × log₁₀(concentration)
        ```
        
        **Reference:** Finney, D.J. (1971). *Probit Analysis*, 3rd ed. Cambridge University Press.
        
        #### **LD Estimation**
        Lethal doses (LD1, LD50, LD99) are estimated from the fitted probit model with 95% confidence intervals 
        calculated using the delta method—a standard approach for propagating uncertainty from regression parameters 
        to derived quantities.
        
        **Formula:** 
        ```
        LDₓ = 10^[(Probit(x) - β₀) / β₁]
        ```
        
        where Probit(x) is the probit value for proportion x (e.g., 0.50 for LD50).
        
        **Delta Method:** The confidence intervals account for uncertainty in both β₀ and β₁ using a Taylor series 
        approximation to propagate their standard errors to the LD estimate.
        
        **References:** Finney, D.J. (1971). *Probit Analysis*, 3rd ed. Cambridge University Press; 
        Ver Hoef, J.M. (2012). Who invented the delta method? *The American Statistician* 66:124-127.
        
        #### **Resistance Ratios**
        Resistance ratios compare LD50 values between two strains using Fieller's method for confidence intervals:
        
        ```
        RR = LD50(Strain 1) / LD50(Strain 2)
        ```
        
        - **RR > 1:** Strain 1 is more resistant
        - **RR < 1:** Strain 1 is more susceptible
        
        **Reference:** Fieller, E.C. (1940). The biological standardization of insulin. *J. R. Stat. Soc. Suppl.* 7:1-64.
        
        #### **Parallelism Test**
        Tests whether two dose-response curves have equal slopes using a Z-test:
        
        ```
        Z = (slope₁ - slope₂) / √(SE₁² + SE₂²)
        p-value = 2 × Φ(-|Z|)
        ```
        
        - **p > 0.05:** Slopes are parallel (similar mode of action)
        - **p < 0.05:** Slopes differ (potentially different mechanisms)
        
        **Interpretation:** Parallel slopes suggest both strains respond to the chemical via the same biological 
        mechanism, just at different dose levels. Non-parallel slopes may indicate different modes of action or 
        resistance mechanisms.
        
        **Reference:** Robertson, J.L., Russell, R.M., Preisler, H.K., and Savin, N.E. (2007). 
        *Bioassays with Arthropods*, 2nd ed. CRC Press.
        
        #### **Equality Test**
        Tests whether two strains have significantly different susceptibility (different intercepts) using a Z-test:
        
        ```
        Z = (intercept₁ - intercept₂) / √(SE₁² + SE₂²)
        ```
        
        - **p > 0.05:** No significant difference in susceptibility
        - **p < 0.05:** Significant difference detected
        
        #### **Goodness-of-Fit Test**
        Pearson's chi-square test evaluates how well the model fits the observed data:
        
        ```
        χ² = Σ[(Observed - Expected)² / Expected]
        ```
        
        - **p > 0.05:** Good fit (model adequately describes the data)
        - **p < 0.05:** Poor fit (high variability or model inadequacy)
        
        **Note:** Poor fit doesn't invalidate results but suggests caution with confidence intervals.
        
        ### 📝 How to Cite Statistical Methods
        
        **Example Methods Section:**
        
        "Probit regression analysis was performed using a generalized linear model with probit link function 
        (Finney, 1971). Lethal dose estimates (LD1, LD50, LD99) and their 95% confidence intervals were calculated 
        using the delta method (Finney, 1971; Ver Hoef, 2012). Resistance ratios were computed as the ratio of LD50 
        values with confidence intervals calculated using Fieller's method (Fieller, 1940). Parallelism of dose-response 
        curves was assessed using a Z-test for equality of slopes (Robertson et al., 2007). Goodness-of-fit was evaluated 
        using Pearson's chi-square test. Statistical significance was set at α = 0.05."
        
        **Tool Citation:**
        ```
        Tidwell, J. (2024). Probit Analysis Tool for Bioassay Research. 
        USDA ARS Cattle Fever Tick Research Unit, Edinburg, TX, USA.
        Available at: [URL]
        ```
        
        ### 📚 Key References
        
        1. **Finney, D.J. (1971).** *Probit Analysis*, 3rd ed. Cambridge University Press.
           - Classic reference for probit methodology and delta method
        
        2. **Robertson, J.L., Russell, R.M., Preisler, H.K., and Savin, N.E. (2007).** 
           *Bioassays with Arthropods*, 2nd ed. CRC Press.
           - Comprehensive guide to bioassay analysis with arthropods
           - Source for parallelism test methodology
        
        3. **Fieller, E.C. (1940).** The biological standardization of insulin. 
           *Journal of the Royal Statistical Society, Supplement* 7:1-64.
           - Original method for ratio confidence intervals
        
        4. **Ver Hoef, J.M. (2012).** Who invented the delta method? 
           *The American Statistician* 66:124-127.
           - Historical review and technical details of delta method
        
        5. **Abbott, W.S. (1925).** A method of computing the effectiveness of an insecticide. 
           *Journal of Economic Entomology* 18:265-267.
           - Abbott's correction for control mortality
        
        ### 💡 Tips
        
        - Use example data to test the tool first
        - Validate data before running analysis
        - Check for high CV% in variability table
        - Read model fit interpretation carefully
        - Save results by downloading when feature is available
        """)

if __name__ == "__main__":
    main()
