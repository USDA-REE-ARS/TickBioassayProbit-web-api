"""
Probit Analysis Tool - Web Application
Streamlit-based web interface for bioassay probit regression analysis
Version: 11.7 Web

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
import hashlib
import re
import json
import urllib.request
import urllib.error
from itertools import combinations
from statsmodels.stats.multitest import multipletests

# Set page configuration
st.set_page_config(
    page_title="Probit Analysis Tool",
    page_icon="📊",
    layout="wide",
    initial_sidebar_state="expanded"
)

# Constants
ALPHA_LEVEL = 0.05
DEFAULT_LC_LEVELS = [1, 50, 99]

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


def parse_lc_levels(value):
    """Parse user-specified LC levels and return sorted unique percentages (0 < x < 100)."""
    if value is None:
        return DEFAULT_LC_LEVELS.copy()

    if isinstance(value, (list, tuple, np.ndarray)):
        raw_values = value
    else:
        raw_values = [part.strip() for part in str(value).replace(';', ',').split(',') if part.strip()]

    levels = []
    for item in raw_values:
        try:
            level = float(item)
        except (TypeError, ValueError):
            raise ValueError(f"'{item}' is not a valid LC percentage.")
        if not (0 < level < 100):
            raise ValueError(f"LC levels must be greater than 0 and less than 100. Invalid value: {level:g}")
        levels.append(level)

    if not levels:
        raise ValueError("Enter at least one LC level (for example: 1, 50, 99).")

    return sorted(set(levels))


def format_lc_level(level):
    """Format an LC level cleanly for display (e.g., 50 instead of 50.0)."""
    level = float(level)
    return str(int(level)) if level.is_integer() else f"{level:g}"


def format_slope(value):
    """Format slope estimates for display without changing calculation precision.

    Ordinary slopes use three decimal places. Very small non-zero values and
    very large values use scientific notation so meaningful magnitude is not
    hidden by fixed-width rounding.
    """
    value = float(value)
    if not np.isfinite(value):
        return "N/A"
    magnitude = abs(value)
    if value != 0 and (magnitude < 0.001 or magnitude >= 1000):
        return f"{value:.3e}"
    return f"{value:.3f}"


def format_p_value(value):
    """Format p-values compactly for display without changing numeric precision.

    Values below 0.001 are shown in scientific notation so very small p-values
    are not rendered as long strings of zeros. A computed value that underflows
    to exactly zero is displayed as <1e-300 rather than as 0.000... .
    """
    try:
        value = float(value)
    except (TypeError, ValueError):
        return "N/A"
    if not np.isfinite(value) or value < 0 or value > 1:
        return "N/A"
    if value == 0:
        return "<1e-300"
    if value < 0.001:
        return f"{value:.3e}"
    return f"{value:.4g}"


def _normalize_column_name(name):
    """Normalize a column heading for tolerant matching."""
    return ''.join(ch for ch in str(name).strip().lower() if ch.isalnum())


COLUMN_ALIASES = {
    'concentration': {
        'concentration', 'conc', 'dose', 'dosage', 'testconcentration',
        'chemicalconcentration', 'treatmentconcentration', 'concentrationppm',
        'ppm', 'percentconcentration', 'concentrationpercent', 'ugml', 'mgml'
    },
    'n': {
        'n', 'ntested', 'tested', 'total', 'totaltested', 'numbertested',
        'samplesize', 'samplecount', 'numberexposed', 'nexposed', 'exposed',
        'totalexposed', 'numberassayed', 'nassayed'
    },
    'mortality': {
        'mortality', 'dead', 'deaths', 'death', 'ndead', 'numberdead',
        'mortalitycount', 'died', 'killed', 'numberkilled', 'nkilled',
        'deadcount', 'numbermortality'
    }
}

METADATA_ALIASES = {
    'strain': {
        'strain', 'strainname', 'population', 'populationname', 'colony',
        'colonyname', 'isolate'
    },
    'chemical': {
        'chemical', 'chemicalname', 'acaricide', 'acaricidename',
        'insecticide', 'compound', 'compoundname', 'activeingredient', 'ai'
    },
    'units': {
        'unit', 'units', 'concentrationunit', 'concentrationunits', 'doseunit',
        'doseunits'
    }
}


def detect_bioassay_columns(df):
    """Return best-effort mapping from canonical bioassay fields to uploaded columns."""
    normalized = {col: _normalize_column_name(col) for col in df.columns}
    mapping = {}
    for canonical, aliases in COLUMN_ALIASES.items():
        # Exact canonical name wins, then normalized aliases.
        for col, norm_name in normalized.items():
            if norm_name == canonical or norm_name in aliases:
                mapping[canonical] = col
                break
    return mapping


def apply_bioassay_column_mapping(df, mapping):
    """Create canonical concentration/n/mortality columns from a user/auto mapping."""
    out = df.copy()
    for canonical in ('concentration', 'n', 'mortality'):
        source = mapping.get(canonical)
        if source is not None and source in out.columns:
            out[canonical] = pd.to_numeric(out[source], errors='coerce')
    return out


def _read_delimited_table(text):
    """Read a delimited text table using delimiter inference with a conservative fallback."""
    try:
        df = pd.read_csv(StringIO(text), sep=None, engine='python')
    except Exception:
        df = pd.read_csv(StringIO(text), sep=r'\s+|,|;|\t', engine='python')
    df.columns = [str(col).strip() for col in df.columns]
    return df


def _line_looks_like_bioassay_header(line):
    """Return True when a line appears to contain the bioassay table header."""
    tokens = [t.strip() for t in re.split(r'\t|,|;|\s{2,}', line.strip()) if t.strip()]
    normalized = {_normalize_column_name(t) for t in tokens}
    hits = 0
    for canonical, aliases in COLUMN_ALIASES.items():
        if canonical in normalized or normalized.intersection(aliases):
            hits += 1
    return hits >= 2


def detect_metadata_columns(df):
    """Extract constant strain/chemical/unit metadata from ordinary table columns when present."""
    normalized = {col: _normalize_column_name(col) for col in df.columns}
    metadata = {'strain': '', 'chemical': '', 'units': ''}
    for field, aliases in METADATA_ALIASES.items():
        source = next((col for col, norm in normalized.items() if norm in aliases), None)
        if source is None:
            continue
        values = df[source].dropna().astype(str).str.strip()
        values = values[values != '']
        unique = list(dict.fromkeys(values.tolist()))
        if len(unique) == 1:
            metadata[field] = unique[0]
    return metadata


def detect_metadata_conflicts(df):
    """Identify recognized metadata columns containing multiple non-empty values.

    A single probit dataset should represent one strain/population, one chemical,
    and one concentration scale. Multiple values in recognized metadata columns
    therefore indicate that observations may have been unintentionally pooled.
    """
    normalized = {col: _normalize_column_name(col) for col in df.columns}
    conflicts = []
    for field, aliases in METADATA_ALIASES.items():
        source = next((col for col, norm in normalized.items() if norm in aliases), None)
        if source is None:
            continue
        values = df[source].dropna().astype(str).str.strip()
        values = values[values != '']
        unique = list(dict.fromkeys(values.tolist()))
        if len(unique) > 1:
            preview = ', '.join(unique[:4])
            if len(unique) > 4:
                preview += ', ...'
            conflicts.append(
                f"Recognized {field} metadata column '{source}' contains multiple values ({preview}). "
                "Filter or split the file so each analysis contains one population, chemical, and concentration scale."
            )
    return conflicts


def read_bioassay_file(content):
    """Read either a legacy two-line-metadata file or a conventional delimited table.

    Legacy format:
        line 1 = strain/population
        line 2 = chemical/acaricide
        line 3 = table header

    Conventional format:
        line 1 = table header
        data may optionally contain constant strain/population, chemical/acaricide,
        and concentration-unit columns.
    """
    lines = [line for line in content.splitlines() if line.strip()]
    if len(lines) < 2:
        raise ValueError("The uploaded file does not contain enough data to identify a table.")

    # Preserve backward compatibility when the third non-empty line is clearly
    # the table header and the first line is not.
    legacy = len(lines) >= 3 and not _line_looks_like_bioassay_header(lines[0]) \
             and _line_looks_like_bioassay_header(lines[2])

    if legacy:
        strain_hint = lines[0].strip()
        chemical_hint = lines[1].strip()
        table_text = '\n'.join(lines[2:])
        df = _read_delimited_table(table_text)
        column_meta = detect_metadata_columns(df)
        units_hint = column_meta.get('units', '')
        source_format = 'Legacy metadata + table'
    else:
        df = _read_delimited_table('\n'.join(lines))
        column_meta = detect_metadata_columns(df)
        strain_hint = column_meta.get('strain', '')
        chemical_hint = column_meta.get('chemical', '')
        units_hint = column_meta.get('units', '')
        source_format = 'Header-first table'

    if df.empty:
        raise ValueError("No bioassay observations were found in the uploaded table.")

    return strain_hint, chemical_hint, units_hint, df, source_format

def validate_bioassay_data(df, strain='', chemical=''):
    """Validate canonicalized bioassay data for structural and biological constraints."""
    errors = []
    warnings = []
    required_cols = ['concentration', 'n', 'mortality']

    missing = [col for col in required_cols if col not in df.columns]
    if missing:
        errors.append(f"Missing required mapped field(s): {', '.join(missing)}")
        return errors, warnings

    for col in required_cols:
        if df[col].isna().any():
            errors.append(f"Column '{col}' contains missing or non-numeric value(s)")
    if errors:
        return errors, warnings

    if (df['concentration'] < 0).any():
        errors.append("Concentration values cannot be negative")
    if (df['n'] <= 0).any():
        errors.append("Number tested (n) must be greater than zero in every row")
    if (df['mortality'] < 0).any():
        errors.append("Mortality counts cannot be negative")

    # Raw n and mortality are counts and therefore should be whole numbers.
    if not np.allclose(df['n'], np.round(df['n'])):
        errors.append("Number tested (n) must contain whole-number counts")
    if not np.allclose(df['mortality'], np.round(df['mortality'])):
        errors.append("Mortality must contain whole-number counts")

    invalid = df[df['mortality'] > df['n']]
    if len(invalid) > 0:
        errors.append(f"Mortality exceeds sample size in {len(invalid)} row(s)")

    positive_concentrations = df.loc[df['concentration'] > 0, 'concentration'].nunique()
    if positive_concentrations < 3:
        errors.append(
            f"Need at least 3 distinct positive treatment concentrations; found {positive_concentrations}"
        )
    elif positive_concentrations < 4:
        warnings.append("Only 3 distinct treatment concentrations were supplied; LC estimates may be imprecise")

    if (df['concentration'] == 0).any():
        warnings.append(
            "Concentration = 0 row(s) detected as untreated controls; these rows are excluded "
            "from the dose-response fit and may be used for Abbott correction"
        )

    # Actual coefficient of variation among replicate mortality percentages.
    treatment = df[df['concentration'] > 0].copy()
    if len(treatment) > 1:
        treatment['mort_pct'] = treatment['mortality'] / treatment['n'] * 100
        by_conc = treatment.groupby('concentration')['mort_pct'].agg(['mean', 'std', 'count'])
        by_conc['cv'] = (by_conc['std'] / by_conc['mean'].replace(0, np.nan)) * 100
        high_var = by_conc[(by_conc['count'] > 1) & (by_conc['cv'] > 20)]
        if len(high_var) > 0:
            warnings.append(f"High replicate variability (CV > 20%) at {len(high_var)} concentration(s)")

    if not str(strain).strip():
        warnings.append("Strain/population name is blank")
    if not str(chemical).strip():
        warnings.append("Chemical/acaricide name is blank")

    return errors, warnings

def detect_control_group(df):
    """Detect an untreated control group (concentration == 0) in uploaded
    bioassay data and compute its pooled mortality rate, for Abbott's
    (1925) correction of background/control mortality. Returns a dict the
    UI and preprocess_data() both use, so detection happens once per
    upload rather than being silently redone (or skipped) at analysis time.
    """
    control_rows = df[(df['concentration'] == 0) & (df['n'] > 0)]
    if len(control_rows) == 0:
        return {
            'has_control': False,
            'n_control_rows': 0,
            'total_control_n': 0,
            'total_control_mortality': 0,
            'control_mortality_pct': None,
        }
    total_n = int(control_rows['n'].sum())
    total_mortality = int(control_rows['mortality'].sum())
    control_pct = (total_mortality / total_n * 100) if total_n > 0 else None
    return {
        'has_control': True,
        'n_control_rows': int(len(control_rows)),
        'total_control_n': total_n,
        'total_control_mortality': total_mortality,
        'control_mortality_pct': control_pct,
    }


def apply_abbott_correction(mortality, n, control_pct):
    """Correct observed treatment mortality counts for background/control
    mortality using Abbott's (1925) formula:

        corrected% = (observed% - control%) / (100 - control%) * 100

    mortality and n are pandas Series of equal length; control_pct is the
    pooled control-group mortality percentage (0-100). Returns corrected
    mortality counts (float, clipped to [0, n]) on the original count scale
    so downstream code can keep treating them like any other mortality
    count.
    """
    if control_pct >= 100:
        raise ValueError(
            "Control mortality is 100%; Abbott's correction is undefined "
            "(there is no surviving control to normalize against)."
        )
    observed_pct = mortality / n * 100
    corrected_pct = (observed_pct - control_pct) / (100 - control_pct) * 100
    corrected_pct = corrected_pct.clip(lower=0, upper=100)
    return corrected_pct / 100 * n


def preprocess_data(df, control_pct=None):
    """Prepare positive-concentration rows for grouped-binomial probit regression.

    Raw 0% and 100% treatment responses are retained for model fitting.  A
    continuity adjustment is used only later when empirical probits are plotted.
    If an untreated control is available, Abbott correction is applied before
    fitting; corrected mortality counts may therefore be fractional.
    """
    df = df[(df['n'] > 0) & (df['concentration'] > 0)].copy()
    df['mortality_raw'] = df['mortality'].copy()

    if control_pct is not None and control_pct > 0:
        df['mortality'] = apply_abbott_correction(df['mortality'], df['n'], control_pct)

    df['alive'] = df['n'] - df['mortality']
    df['log_concentration'] = np.log10(df['concentration'])
    return df


def fit_probit_model(df):
    """Fit a grouped-binomial GLM with probit link and Pearson goodness-of-fit."""
    X = sm.add_constant(df['log_concentration'])
    y = np.column_stack([df['mortality'].to_numpy(), df['alive'].to_numpy()])

    try:
        model = sm.GLM(y, X, family=sm.families.Binomial(link=sm.families.links.Probit()))
        result = model.fit()
        if not result.converged:
            raise RuntimeError(
                "Model did not converge. Check concentration spacing, response range, "
                "replication, and variability."
            )
        if not np.isfinite(result.params).all():
            raise RuntimeError("Model produced non-finite regression parameters.")
        if not np.isfinite(result.bse).all() or (result.bse <= 0).any():
            raise RuntimeError("Model produced invalid parameter standard errors.")
    except Exception as e:
        if isinstance(e, RuntimeError):
            raise
        raise RuntimeError(
            f"Model fitting failed: {e}. Check data format, treatment range, and mortality counts."
        )

    predicted_prob = result.predict(X)
    if not np.isfinite(predicted_prob).all() or (predicted_prob < 0).any() or (predicted_prob > 1).any():
        raise RuntimeError("Model produced invalid probability predictions.")

    # statsmodels computes the Pearson statistic using the binomial variance,
    # which is the appropriate grouped-binomial GLM goodness-of-fit statistic.
    chi_square = float(result.pearson_chi2)
    df_degrees = int(result.df_resid)
    if df_degrees <= 0:
        raise RuntimeError("Insufficient residual degrees of freedom for goodness-of-fit assessment.")
    p_value = float(chi2.sf(chi_square, df_degrees))
    return result, chi_square, df_degrees, p_value


def assess_lc_estimate(df_processed, level, estimate):
    """Assess whether an LC estimate is interpolation or extrapolation.

    The primary classification is based on whether the estimated concentration
    lies within the tested positive-concentration range. Response coverage is
    reported separately because replicate variation can make the observed
    mortality range non-monotonic even when the estimate lies within the tested
    concentration range.
    """
    observed_pct = (df_processed['mortality'] / df_processed['n']) * 100
    response_low = float(observed_pct.min())
    response_high = float(observed_pct.max())
    conc_low = float(df_processed['concentration'].min())
    conc_high = float(df_processed['concentration'].max())
    response_covered = response_low <= float(level) <= response_high
    concentration_covered = conc_low <= float(estimate) <= conc_high
    status = 'Interpolated' if concentration_covered else 'Extrapolated'
    return {
        'status': status,
        'response_covered': response_covered,
        'response_low': response_low,
        'response_high': response_high,
        'concentration_low': conc_low,
        'concentration_high': conc_high,
    }


def compute_lcx_with_ci(result, x, alpha=ALPHA_LEVEL):
    """Compute LCx with confidence intervals with robust error checking"""
    intercept, slope = result.params
    
    # CRITICAL FIX: Validate model parameters
    if not np.isfinite(intercept) or not np.isfinite(slope):
        raise ValueError("Model parameters are not finite. Check data quality and model convergence.")
    
    # LC estimates are meaningful here only for an increasing concentration-mortality relationship.
    if slope <= 0:
        raise ValueError(
            f"LC estimates are not reported because the fitted slope is non-positive ({format_slope(slope)}). "
            "Review concentration coding, assay response, and data quality."
        )
    
    target_quantile = norm.ppf(x / 100.0)
    
    # Point estimate with validation
    log_conc_lc = (target_quantile - intercept) / slope
    ld = 10 ** log_conc_lc
    
    # CRITICAL FIX: Validate LC estimate
    if not np.isfinite(ld) or ld <= 0:
        raise ValueError(f"Invalid LC{format_lc_level(x)} estimate: {ld}. " +
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
    
    var_log_lc = (grad_intercept ** 2) * var_intercept + \
                  (grad_slope ** 2) * var_slope + \
                  2 * grad_intercept * grad_slope * cov_intercept_slope
    
    # CRITICAL FIX: Check for negative variance
    if var_log_lc < 0:
        raise ValueError("Negative variance calculated in LC estimation. Check model covariance matrix.")
    
    se_log_lc = np.sqrt(var_log_lc)
    z_critical = norm.ppf(1 - alpha / 2)
    
    log_lc_lower = log_conc_lc - z_critical * se_log_lc
    log_lc_upper = log_conc_lc + z_critical * se_log_lc
    
    lc_lower = 10 ** log_lc_lower
    lc_upper = 10 ** log_lc_upper
    
    # CRITICAL FIX: Validate confidence interval bounds
    if not (np.isfinite(lc_lower) and np.isfinite(lc_upper) and lc_lower > 0 and lc_upper > 0):
        raise ValueError(f"Invalid confidence interval bounds: [{lc_lower}, {lc_upper}]")
    
    return ld, lc_lower, lc_upper

def calculate_r_squared(result, df_processed):
    """Calculate R-squared for probit regression model"""
    # Get predicted probabilities
    X = sm.add_constant(df_processed['log_concentration'])
    predicted_prob = result.predict(X)
    predicted_prob = np.clip(np.asarray(predicted_prob, dtype=float), 1e-6, 1 - 1e-6)
    
    # Convert to probit scale for R² calculation
    observed_mortality_pct = (df_processed['mortality'] / df_processed['n']) * 100
    # Continuity correction for 0%/100% boundary values so ppf() doesn't
    # return +/- infinity. (df_processed['mortality'] already reflects
    # Abbott's correction for control mortality, if one was applied.)
    observed_mortality_pct_adj = observed_mortality_pct.copy()
    observed_mortality_pct_adj[observed_mortality_pct_adj <= 0] = 0.1
    observed_mortality_pct_adj[observed_mortality_pct_adj >= 100] = 99.9
    observed_probits = norm.ppf(observed_mortality_pct_adj / 100)
    
    # Predicted probits
    predicted_probits = norm.ppf(predicted_prob)
    
    # Calculate R²
    ss_res = np.sum((observed_probits - predicted_probits) ** 2)
    ss_tot = np.sum((observed_probits - np.mean(observed_probits)) ** 2)
    
    r_squared = 1 - (ss_res / ss_tot) if ss_tot > 0 else 0
    
    # Ensure R² is between 0 and 1
    r_squared = max(0, min(1, r_squared))
    
    return r_squared

def interpret_slope_biology(slope, chemical="chemical"):
    """Provide a conservative interpretation of dose-response steepness."""
    magnitude = abs(float(slope))
    if slope <= 0:
        return (
            f"The fitted slope is non-positive (slope = {format_slope(slope)}), which is atypical for a standard "
            "increasing concentration-mortality relationship. Review concentration coding and assay data."
        )
    if magnitude >= 10:
        descriptor = "very steep"
    elif magnitude >= 5:
        descriptor = "steep"
    elif magnitude >= 2:
        descriptor = "moderate"
    else:
        descriptor = "shallow"
    return (
        f"The {chemical} concentration-response is {descriptor} (probit slope = {format_slope(slope)}). "
        "Slope describes how rapidly mortality changes across log10 concentration and can reflect "
        "response heterogeneity and experimental precision; it does not by itself identify a molecular "
        "target, mode of action, or resistance mechanism."
    )


def interpret_r_squared_biology(r_squared):
    """Describe the auxiliary probit-scale R-squared without treating it as GLM goodness-of-fit."""
    return (
        f"Descriptive probit-scale R² = {r_squared:.3f}. This summarizes how closely empirical probits "
        "track the fitted probit line and is provided as a descriptive aid only. Formal model-fit "
        "assessment in this application uses the Pearson chi-square statistic and its residual degrees of freedom."
    )

def compute_resistance_ratio_ci(result1, result2, lc_level=50, alpha=ALPHA_LEVEL):
    """
    Compute resistance ratio (LCx_1 / LCx_2) with a delta-method 95% CI,
    consistent with the log10-scale approach used in compute_lcx_with_ci().

    Both LC estimates are computed on the log10(concentration) scale, so the
    ratio's confidence interval is derived by combining the two independent
    log10-scale variances rather than mixing log10 and natural-log scales.
    """
    intercept1, slope1 = result1.params
    intercept2, slope2 = result2.params

    if not (np.isfinite(intercept1) and np.isfinite(slope1) and
            np.isfinite(intercept2) and np.isfinite(slope2)):
        raise ValueError("Model parameters are not finite. Check data quality and model convergence.")

    if slope1 <= 0 or slope2 <= 0:
        raise ValueError(
            "LC50 ratios are not reported when either fitted concentration-response slope is non-positive."
        )

    target_quantile = norm.ppf(lc_level / 100.0)

    # Point estimates on log10(concentration) scale
    log_conc_lc1 = (target_quantile - intercept1) / slope1
    log_conc_lc2 = (target_quantile - intercept2) / slope2

    # Delta-method variance for each LC estimate (log10 scale) - same derivation as compute_lcx_with_ci
    cov1 = result1.cov_params()
    var_intercept1 = cov1.iloc[0, 0]
    var_slope1 = cov1.iloc[1, 1]
    cov_is1 = cov1.iloc[0, 1]

    cov2 = result2.cov_params()
    var_intercept2 = cov2.iloc[0, 0]
    var_slope2 = cov2.iloc[1, 1]
    cov_is2 = cov2.iloc[0, 1]

    if not all(np.isfinite([var_intercept1, var_slope1, cov_is1, var_intercept2, var_slope2, cov_is2])):
        raise ValueError("Invalid covariance matrix values. Check model fitting.")

    grad_intercept1 = -1 / slope1
    grad_slope1 = -(target_quantile - intercept1) / (slope1 ** 2)
    var_log_lc1 = (grad_intercept1 ** 2) * var_intercept1 + \
                  (grad_slope1 ** 2) * var_slope1 + \
                  2 * grad_intercept1 * grad_slope1 * cov_is1

    grad_intercept2 = -1 / slope2
    grad_slope2 = -(target_quantile - intercept2) / (slope2 ** 2)
    var_log_lc2 = (grad_intercept2 ** 2) * var_intercept2 + \
                  (grad_slope2 ** 2) * var_slope2 + \
                  2 * grad_intercept2 * grad_slope2 * cov_is2

    if var_log_lc1 < 0 or var_log_lc2 < 0:
        raise ValueError("Negative variance calculated in LC estimation. Check model covariance matrix.")

    # log10(ratio) = log10(LCx_1) - log10(LCx_2); variances add for independent samples
    log10_rr = log_conc_lc1 - log_conc_lc2
    se_log10_rr = np.sqrt(var_log_lc1 + var_log_lc2)

    if not np.isfinite(se_log10_rr) or se_log10_rr < 0:
        raise ValueError("Could not calculate a valid standard error for the resistance ratio.")

    z_critical = norm.ppf(1 - alpha / 2)

    rr = 10 ** log10_rr
    rr_lower = 10 ** (log10_rr - z_critical * se_log10_rr)
    rr_upper = 10 ** (log10_rr + z_critical * se_log10_rr)

    if not (np.isfinite(rr) and np.isfinite(rr_lower) and np.isfinite(rr_upper) and rr > 0):
        raise ValueError("Resistance ratio confidence interval calculation produced invalid bounds.")

    return rr, rr_lower, rr_upper


def compute_lc50_fold_difference(result1, result2, name1='Dataset 1', name2='Dataset 2',
                                 lc_level=50, alpha=ALPHA_LEVEL):
    """Return an unsigned LCx fold-difference with the smaller LCx in the denominator.

    This is intended for non-directional pairwise comparisons. The returned
    fold-difference is always >= 1. When the directional result1/result2 ratio
    is < 1, the ratio and confidence interval are inverted so the larger LCx is
    divided by the smaller LCx:

        fold = 1 / ratio
        CI = (1 / upper, 1 / lower)

    Directional susceptible-reference resistance ratios should continue to use
    compute_resistance_ratio_ci() without inversion.
    """
    ratio, lower, upper = compute_resistance_ratio_ci(
        result1, result2, lc_level=lc_level, alpha=alpha
    )

    if ratio >= 1:
        return {
            'fold_difference': float(ratio),
            'lower': float(lower),
            'upper': float(upper),
            'higher_name': str(name1),
            'lower_name': str(name2),
            'inverted': False,
        }

    return {
        'fold_difference': float(1.0 / ratio),
        'lower': float(1.0 / upper),
        'upper': float(1.0 / lower),
        'higher_name': str(name2),
        'lower_name': str(name1),
        'inverted': True,
    }


def compare_dose_response_models(df1, df2):
    """Compare two grouped-binomial probit curves with nested combined GLMs.

    The full model contains dataset, log10(concentration), and their interaction.
    A likelihood-ratio test of full vs. parallel models tests the interaction
    (slope difference). If a common-slope representation is reasonable, a second
    likelihood-ratio test of parallel vs. common models tests the dataset shift.
    """
    a = df1[['log_concentration', 'mortality', 'alive']].copy()
    b = df2[['log_concentration', 'mortality', 'alive']].copy()
    a['dataset'] = 0.0
    b['dataset'] = 1.0
    combined = pd.concat([a, b], ignore_index=True)
    combined['interaction'] = combined['log_concentration'] * combined['dataset']
    y = np.column_stack([combined['mortality'].to_numpy(), combined['alive'].to_numpy()])

    x_full = sm.add_constant(combined[['log_concentration', 'dataset', 'interaction']], has_constant='add')
    x_parallel = sm.add_constant(combined[['log_concentration', 'dataset']], has_constant='add')
    x_common = sm.add_constant(combined[['log_concentration']], has_constant='add')

    family = sm.families.Binomial(link=sm.families.links.Probit())
    full = sm.GLM(y, x_full, family=family).fit()
    parallel = sm.GLM(y, x_parallel, family=family).fit()
    common = sm.GLM(y, x_common, family=family).fit()
    if not (full.converged and parallel.converged and common.converged):
        raise RuntimeError("Combined comparison model did not converge.")

    lr_slope = max(0.0, 2.0 * (full.llf - parallel.llf))
    p_slope = float(chi2.sf(lr_slope, 1))
    lr_shift = max(0.0, 2.0 * (parallel.llf - common.llf))
    p_shift = float(chi2.sf(lr_shift, 1))

    return {
        'full_result': full,
        'parallel_result': parallel,
        'common_result': common,
        'lr_slope': lr_slope,
        'p_slope': p_slope,
        'lr_shift': lr_shift,
        'p_shift': p_shift,
    }


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

def create_probit_plot(df, result, strain, chemical, lc_levels=None):
    """Create probit transformation plot"""
    lc_levels = lc_levels or DEFAULT_LC_LEVELS
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # Calculate empirical probits
    mortality_pct = (df['mortality'] / df['n']) * 100
    # Continuity correction for 0%/100% boundary values so ppf() doesn't
    # return +/- infinity. (df['mortality'] already reflects Abbott's
    # correction for control mortality, if one was applied.)
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
    fitted_probits = result.params.iloc[0] + result.params.iloc[1] * log_conc_range
    
    ax.plot(log_conc_range, fitted_probits,
           color='blue', linewidth=2, label=f'{strain} (Fitted)', zorder=2)
    
    # Reference lines for user-selected LC levels
    for lc_level in lc_levels:
        probit_value = norm.ppf(lc_level / 100.0)
        ax.axhline(y=probit_value, color='gray', linestyle='--', alpha=0.3, linewidth=1)
        ax.text(df['log_concentration'].max(), probit_value, f' LC{format_lc_level(lc_level)}',
               verticalalignment='center', fontsize=9, alpha=0.7)
    
    ax.set_xlabel(f'Log10({chemical} Concentration)', fontsize=12, fontweight='bold')
    ax.set_ylabel('Probit (Mortality)', fontsize=12, fontweight='bold')
    ax.set_title(f'Probit Regression Line\n{strain} - {chemical}',
                fontsize=14, fontweight='bold')
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=10)
    
    plt.tight_layout()
    return fig

def pdf_to_bytes(pdf):
    """Render an FPDF document to bytes in a way that works with fpdf2 and PyFPDF.

    fpdf2 returns bytes/bytearray from ``output()`` whereas legacy PyFPDF may
    write to an output destination instead of returning the document payload.
    Writing to a temporary file and reading it back avoids package-version
    differences and prevents zero-byte/blank Streamlit downloads.
    """
    temp_path = None
    try:
        with tempfile.NamedTemporaryFile(delete=False, suffix='.pdf') as tmpfile:
            temp_path = tmpfile.name

        # Both fpdf2 and legacy PyFPDF support writing directly to a filename.
        pdf.output(temp_path)

        with open(temp_path, 'rb') as handle:
            data = handle.read()

        if not data:
            raise RuntimeError('PDF generation produced an empty file.')
        if not data.startswith(b'%PDF'):
            raise RuntimeError('PDF generation produced an invalid file header.')
        return data
    finally:
        if temp_path and os.path.exists(temp_path):
            try:
                os.unlink(temp_path)
            except OSError:
                pass


def _write_control_info(pdf, control_info):
    if not control_info or not control_info.get('has_control'):
        pdf.cell(0, 6, 'Untreated control: not supplied', 0, 1)
        pdf.cell(0, 6, "Abbott correction: not applied", 0, 1)
        return
    pct = control_info.get('control_mortality_pct')
    deaths = control_info.get('total_control_mortality', 0)
    total = control_info.get('total_control_n', 0)
    pdf.cell(0, 6, f'Untreated control mortality: {pct:.2f}% ({deaths}/{total})', 0, 1)
    pdf.cell(0, 6, f"Abbott correction: {'applied' if pct and pct > 0 else 'not required (0% control mortality)'}", 0, 1)


def create_pdf_single_with_plots(df, df_processed, result, chi_square, df_degrees, p_value,
                                 strain, chemical, options=None, lc_levels=None,
                                 units='', control_info=None):
    """Create a customizable PDF report for one dataset."""
    lc_levels = lc_levels or DEFAULT_LC_LEVELS
    options = options or {
        'include_data_summary': True, 'include_raw_data': True,
        'include_lc_estimates': True, 'include_model_fit': True,
        'include_parameters': True, 'include_mortality_plot': True,
        'include_probit_plot': True
    }
    pdf = FPDF()
    pdf.set_auto_page_break(auto=True, margin=15)
    pdf.add_page()
    pdf.set_font('Arial', 'B', 16)
    pdf.cell(0, 10, 'Probit Analysis Report', 0, 1, 'C')
    pdf.set_font('Arial', 'B', 12)
    pdf.cell(0, 8, f'Strain/Population: {strain or "Not specified"}', 0, 1)
    pdf.cell(0, 8, f'Chemical/Acaricide: {chemical or "Not specified"}', 0, 1)
    if units:
        pdf.cell(0, 8, f'Concentration units: {units}', 0, 1)
    pdf.cell(0, 8, f'Date: {datetime.now().strftime("%Y-%m-%d %H:%M")}', 0, 1)
    pdf.ln(3)

    if options.get('include_data_summary', True):
        pdf.set_font('Arial', 'B', 14)
        pdf.cell(0, 9, 'Data Summary', 0, 1)
        pdf.set_font('Arial', '', 10)
        pdf.cell(0, 6, f'Rows supplied: {len(df)}', 0, 1)
        pdf.cell(0, 6, f'Treatment rows fitted: {len(df_processed)}', 0, 1)
        pdf.cell(0, 6, f'Total individuals supplied: {int(df["n"].sum())}', 0, 1)
        pdf.cell(0, 6, f'Total raw mortality supplied: {int(df["mortality"].sum())}', 0, 1)
        _write_control_info(pdf, control_info)
        pdf.ln(3)

    if options.get('include_raw_data', True):
        pdf.set_font('Arial', 'B', 14)
        pdf.cell(0, 9, 'Raw Data', 0, 1)
        pdf.set_font('Arial', 'B', 9)
        for width, label in [(45, 'Concentration'), (35, 'N Tested'), (35, 'Mortality'), (45, 'Mortality %')]:
            pdf.cell(width, 6, label, 1)
        pdf.ln()
        pdf.set_font('Arial', '', 9)
        for _, row in df.iterrows():
            mort_pct = row['mortality'] / row['n'] * 100
            pdf.cell(45, 5, f"{row['concentration']:.6g}", 1)
            pdf.cell(35, 5, f"{int(row['n'])}", 1)
            pdf.cell(35, 5, f"{int(row['mortality'])}", 1)
            pdf.cell(45, 5, f"{mort_pct:.2f}%", 1)
            pdf.ln()
        pdf.ln(3)

    if options.get('include_lc_estimates', True):
        pdf.set_font('Arial', 'B', 14)
        pdf.cell(0, 9, 'Lethal Concentration Estimates (95% CI)', 0, 1)
        pdf.set_font('Arial', '', 10)
        for level in lc_levels:
            est, lower, upper = compute_lcx_with_ci(result, level)
            coverage = assess_lc_estimate(df_processed, level, est)
            response_note = 'response covered' if coverage['response_covered'] else 'response target outside observed range'
            pdf.cell(0, 6, f"LC{format_lc_level(level)}: {est:.6g} ({lower:.6g} - {upper:.6g}) [{coverage['status']}; {response_note}]", 0, 1)
            if coverage['status'] == 'Extrapolated' or not coverage['response_covered']:
                pdf.set_font('Arial', 'I', 8)
                pdf.multi_cell(
                    0, 4,
                    f"  Tested concentration range {coverage['concentration_low']:.6g}-{coverage['concentration_high']:.6g}; "
                    f"observed mortality range {coverage['response_low']:.1f}-{coverage['response_high']:.1f}%."
                )
                pdf.set_font('Arial', '', 10)
        pdf.ln(3)

    if options.get('include_model_fit', True):
        pdf.set_font('Arial', 'B', 14)
        pdf.cell(0, 9, 'Model Fit Statistics', 0, 1)
        pdf.set_font('Arial', '', 10)
        pdf.cell(0, 6, f'Pearson chi-square: {chi_square:.4f}', 0, 1)
        pdf.cell(0, 6, f'Residual degrees of freedom: {df_degrees}', 0, 1)
        pdf.cell(0, 6, f'p-value: {format_p_value(p_value)}', 0, 1)
        dispersion = chi_square / df_degrees if df_degrees > 0 else float('nan')
        pdf.cell(0, 6, f'Pearson dispersion (chi-square/df): {dispersion:.3f}', 0, 1)
        pdf.multi_cell(0, 5, 'Interpretation: p < 0.05 indicates evidence of lack of fit/extra-binomial variation; p >= 0.05 does not demonstrate that the model is true. Dispersion >1 indicates residual variation above the nominal binomial expectation.')
        pdf.ln(3)

    if options.get('include_parameters', True):
        pdf.set_font('Arial', 'B', 14)
        pdf.cell(0, 9, 'Model Parameters', 0, 1)
        pdf.set_font('Arial', '', 10)
        pdf.cell(0, 6, f'Intercept: {result.params.iloc[0]:.4f} (SE {result.bse.iloc[0]:.4f})', 0, 1)
        pdf.cell(0, 6, f'Slope: {format_slope(result.params.iloc[1])} (SE {format_slope(result.bse.iloc[1])})', 0, 1)
        try:
            r2 = calculate_r_squared(result, df_processed)
            pdf.cell(0, 6, f'Descriptive probit-scale R-squared: {r2:.4f}', 0, 1)
        except Exception:
            pass
        pdf.ln(3)

    if options.get('include_mortality_plot', True):
        pdf.add_page()
        pdf.set_font('Arial', 'B', 14)
        pdf.cell(0, 10, 'Mortality Curve', 0, 1)
        with tempfile.NamedTemporaryFile(delete=False, suffix='.png') as tmpfile:
            fig = create_mortality_plot(df_processed, result, strain, chemical)
            fig.savefig(tmpfile.name, format='png', dpi=150, bbox_inches='tight')
            plt.close(fig)
            pdf.image(tmpfile.name, x=10, w=190)
            os.unlink(tmpfile.name)

    if options.get('include_probit_plot', True):
        pdf.add_page()
        pdf.set_font('Arial', 'B', 14)
        pdf.cell(0, 10, 'Probit Regression Line', 0, 1)
        with tempfile.NamedTemporaryFile(delete=False, suffix='.png') as tmpfile:
            fig = create_probit_plot(df_processed, result, strain, chemical, lc_levels=lc_levels)
            fig.savefig(tmpfile.name, format='png', dpi=150, bbox_inches='tight')
            plt.close(fig)
            pdf.image(tmpfile.name, x=10, w=190)
            os.unlink(tmpfile.name)
    return pdf


def render_column_mapping(df_raw, dataset_key):
    """Render column mapping controls and return canonical data plus a validity flag."""
    auto_mapping = detect_bioassay_columns(df_raw)
    columns = list(df_raw.columns)
    placeholder = "— Select column —"
    options = [placeholder] + columns
    labels = {
        'concentration': 'Concentration column',
        'n': 'Number tested (n) column',
        'mortality': 'Mortality / deaths column'
    }
    missing_auto = [field for field in ('concentration', 'n', 'mortality') if field not in auto_mapping]
    with st.expander("🔧 Column Mapping", expanded=bool(missing_auto)):
        st.caption("Common headings are detected automatically. Confirm or change each assignment.")
        chosen = {}
        for field in ('concentration', 'n', 'mortality'):
            detected = auto_mapping.get(field)
            default_index = options.index(detected) if detected in columns else 0
            source = st.selectbox(labels[field], options, index=default_index, key=f"{dataset_key}_map_{field}")
            if source != placeholder:
                chosen[field] = source
        selected_sources = list(chosen.values())
        complete = len(chosen) == 3
        unique = len(selected_sources) == len(set(selected_sources))
        mapping_valid = complete and unique
        if not complete:
            st.warning("Map all three required fields before analysis.")
        if not unique:
            st.error("Each required field must map to a different uploaded column.")
        detected_text = ", ".join(f"{field} → {source}" for field, source in auto_mapping.items()) or "none"
        st.caption(f"Auto-detected mapping: {detected_text}")
    return apply_bioassay_column_mapping(df_raw, chosen), chosen, mapping_valid


def fig_to_base64(fig):
    """Convert matplotlib figure to base64 string"""
    buf = BytesIO()
    fig.savefig(buf, format='png', dpi=150, bbox_inches='tight')
    buf.seek(0)
    img_str = base64.b64encode(buf.read()).decode()
    buf.close()
    return img_str

# -----------------------------------------------------------------------------
# Version 11 multi-dataset application helpers
# -----------------------------------------------------------------------------


def _decode_uploaded_text(uploaded_file):
    """Decode an uploaded text file with a conservative fallback."""
    raw = uploaded_file.getvalue()
    try:
        return raw.decode('utf-8-sig'), None
    except UnicodeDecodeError:
        return raw.decode('latin-1'), "File was decoded as Latin-1 because UTF-8 decoding failed."


def clear_dynamic_dataset_state(dataset_key):
    """Clear widget/session keys associated with one dynamic multi-upload dataset."""
    prefixes = [
        f'{dataset_key}_file_signature',
        f'{dataset_key}_strain_input',
        f'{dataset_key}_chemical_input',
        f'{dataset_key}_map_concentration',
        f'{dataset_key}_map_n',
        f'{dataset_key}_map_mortality',
    ]
    for key in prefixes:
        st.session_state.pop(key, None)


def render_uploaded_dataset_dynamic(uploaded_file, dataset_key):
    """Render one dataset in the v11 multi-file uploader and return a record dict."""
    content, decode_warning = _decode_uploaded_text(uploaded_file)
    strain_hint, chemical_hint, units_hint, df_raw, source_format = read_bioassay_file(content)
    metadata_conflicts = detect_metadata_conflicts(df_raw)
    signature = hashlib.sha256(uploaded_file.getvalue()).hexdigest()[:16]
    sig_key = f'{dataset_key}_file_signature'

    if st.session_state.get(sig_key) != signature:
        st.session_state[sig_key] = signature
        st.session_state[f'{dataset_key}_strain_input'] = strain_hint
        st.session_state[f'{dataset_key}_chemical_input'] = chemical_hint
        for field in ('concentration', 'n', 'mortality'):
            st.session_state.pop(f'{dataset_key}_map_{field}', None)

    with st.expander(f"📄 {uploaded_file.name}", expanded=True):
        st.caption(f"Detected file structure: {source_format}")
        if decode_warning:
            st.warning(decode_warning)

        col_a, col_b = st.columns(2)
        with col_a:
            strain = st.text_input(
                "Strain / population",
                key=f'{dataset_key}_strain_input',
                help="Auto-populated when found in the file; edit or enter manually as needed."
            ).strip()
        with col_b:
            chemical = st.text_input(
                "Chemical / acaricide",
                key=f'{dataset_key}_chemical_input',
                help="Auto-populated when found in the file; edit or enter manually as needed."
            ).strip()

        if units_hint:
            st.caption(f"Detected concentration units in file: {units_hint}")
        else:
            st.caption(
                "Concentration units were not detected in this file. Enter them once in the shared "
                "Concentration Units section below the uploaded datasets."
            )

        df, mapping, mapping_valid = render_column_mapping(df_raw, dataset_key)
        if mapping_valid:
            errors, warnings = validate_bioassay_data(df, strain, chemical)
        else:
            errors, warnings = (["Column mapping is incomplete or duplicated"], [])
        errors.extend(metadata_conflicts)
        if decode_warning:
            warnings.append(decode_warning)

        if errors:
            st.error("❌ Validation errors")
            for error in errors:
                st.error(f"• {error}")
        else:
            st.success("✓ Data validation passed")

        if warnings:
            st.warning("⚠️ Warnings")
            for warning in warnings:
                st.warning(f"• {warning}")

        control_info = {'has_control': False}
        canonical_ready = (
            mapping_valid
            and all(col in df.columns for col in ('concentration', 'n', 'mortality'))
            and not df[['concentration', 'n', 'mortality']].isna().any().any()
        )
        if canonical_ready:
            control_info = detect_control_group(df)
            if control_info['has_control']:
                ctrl_pct = control_info['control_mortality_pct']
                st.info(
                    f"🧪 Untreated control: {control_info['n_control_rows']} row(s), pooled mortality "
                    f"{ctrl_pct:.1f}% ({control_info['total_control_mortality']}/{control_info['total_control_n']}). "
                    "Abbott correction is applied when control mortality is >0%."
                )
                if ctrl_pct > 20:
                    st.warning(
                        f"Control mortality is {ctrl_pct:.1f}%. This exceeds a commonly used 20% "
                        "acceptability threshold; consult the protocol governing this assay."
                    )
            else:
                st.caption("No concentration = 0 untreated-control row detected; raw treatment mortality will be modeled.")

        if mapping_valid and not errors:
            with st.expander("👁️ Preview canonical data", expanded=False):
                preview = df.copy()
                preview['mortality_%'] = (preview['mortality'] / preview['n'] * 100).round(2)
                st.dataframe(preview, use_container_width=True)

    return {
        'key': dataset_key,
        'signature': signature,
        'file_name': uploaded_file.name,
        'df': df,
        'strain': strain,
        'chemical': chemical,
        'units': units_hint,
        'units_detected': units_hint,
        'units_source': 'Embedded in file' if units_hint else 'Shared entry required',
        'errors': errors,
        'warnings': warnings,
        'control_info': control_info,
        'source_format': source_format,
        'mapping': mapping,
    }


def dataset_display_name(record):
    """Return a stable human-readable name for a dataset record."""
    strain = (record.get('strain') or '').strip()
    return strain if strain else record.get('file_name', record.get('key', 'Dataset'))


def dataset_sort_key(record):
    """Case-insensitive alphabetical sort key for dataset records.

    Strain/population is the primary key; file name and internal key provide
    deterministic tie-breakers when display names are duplicated.
    """
    display = dataset_display_name(record)
    file_name = str(record.get('file_name', ''))
    internal_key = str(record.get('key', ''))
    return (display.casefold(), file_name.casefold(), internal_key.casefold())


def sort_dataset_records(records):
    """Return dataset records in stable alphabetical display order."""
    return sorted(records, key=dataset_sort_key)


def sorted_analysis_items(analyses):
    """Return (key, analysis) pairs in alphabetical dataset order."""
    return sorted(
        analyses.items(),
        key=lambda item: dataset_sort_key(item[1]['record'])
    )


def fit_dataset_record(record, lc_levels):
    """Fit one validated dataset and return rich analysis objects and serializable summaries."""
    if record.get('errors'):
        raise ValueError("Dataset has validation errors and cannot be analyzed.")

    control_info = record.get('control_info') or {'has_control': False}
    control_pct = control_info.get('control_mortality_pct') if control_info.get('has_control') else None
    df_processed = preprocess_data(record['df'], control_pct=control_pct)
    if len(df_processed) == 0:
        raise ValueError("No positive-concentration observations remain after preprocessing.")

    result, chi_square, df_degrees, p_value = fit_probit_model(df_processed)
    intercept = float(result.params.iloc[0])
    slope = float(result.params.iloc[1])
    dispersion = float(chi_square / df_degrees) if df_degrees > 0 else np.nan
    try:
        r_squared = float(calculate_r_squared(result, df_processed))
    except Exception:
        r_squared = np.nan

    lc_rows = []
    lc_error = None
    if slope > 0:
        for level in lc_levels:
            estimate, lower, upper = compute_lcx_with_ci(result, level)
            coverage = assess_lc_estimate(df_processed, level, estimate)
            lc_rows.append({
                'level': float(level),
                'label': f"LC{format_lc_level(level)}",
                'estimate': float(estimate),
                'lower': float(lower),
                'upper': float(upper),
                'status': coverage['status'],
                'response_covered': bool(coverage['response_covered']),
                'response_low': float(coverage['response_low']),
                'response_high': float(coverage['response_high']),
                'concentration_low': float(coverage['concentration_low']),
                'concentration_high': float(coverage['concentration_high']),
            })
    else:
        lc_error = "LC estimates are suppressed because the fitted concentration-response slope is non-positive."

    return {
        'key': record['key'],
        'record': record,
        'df_processed': df_processed,
        'result': result,
        'chi_square': float(chi_square),
        'df_degrees': int(df_degrees),
        'p_value': float(p_value),
        'dispersion': dispersion,
        'intercept': intercept,
        'slope': slope,
        'r_squared': r_squared,
        'lc_rows': lc_rows,
        'lc_error': lc_error,
    }


def serializable_analysis_summary(analysis):
    """Convert an analysis result to a compact structure suitable for session state/AI context."""
    record = analysis['record']
    return {
        'dataset_key': analysis['key'],
        'strain_population': record.get('strain', ''),
        'chemical_acaricide': record.get('chemical', ''),
        'concentration_units': record.get('units', ''),
        'n_rows_supplied': int(len(record['df'])),
        'n_treatment_rows_fitted': int(len(analysis['df_processed'])),
        'slope': analysis['slope'],
        'intercept': analysis['intercept'],
        'pearson_chi_square': analysis['chi_square'],
        'residual_df': analysis['df_degrees'],
        'pearson_p_value': analysis['p_value'],
        'pearson_dispersion': analysis['dispersion'],
        'descriptive_probit_r_squared': None if not np.isfinite(analysis['r_squared']) else analysis['r_squared'],
        'lc_estimates': analysis['lc_rows'],
        'lc_warning': analysis['lc_error'],
        'control': analysis['record'].get('control_info', {}),
        'validation_warnings': analysis['record'].get('warnings', []),
    }


def run_dataset_analyses(records, lc_levels):
    """Fit all supplied records, returning successes and failures keyed by dataset key."""
    analyses = {}
    failures = {}
    for record in sort_dataset_records(records):
        try:
            analyses[record['key']] = fit_dataset_record(record, lc_levels)
        except Exception as exc:
            failures[record['key']] = str(exc)
    return analyses, failures


def compare_multi_dose_response_models(analyses):
    """Global multi-population nested grouped-binomial probit comparison.

    Full model:     logC + dataset + logC:dataset
    Parallel model: logC + dataset
    Common model:   logC

    Full vs parallel tests the global dataset-by-concentration interaction.
    Parallel vs common tests the global dataset shift under a common-slope model.
    """
    if len(analyses) < 2:
        raise ValueError("At least two successfully fitted datasets are required for a global comparison.")

    parts = []
    ordered_keys = [key for key, _ in sorted_analysis_items(analyses)]
    for key in ordered_keys:
        data = analyses[key]['df_processed'][['log_concentration', 'mortality', 'alive']].copy()
        data['dataset'] = key
        parts.append(data)
    combined = pd.concat(parts, ignore_index=True)

    # Stable category order keeps the first selected dataset as the baseline.
    combined['dataset'] = pd.Categorical(combined['dataset'], categories=ordered_keys, ordered=True)
    dummies = pd.get_dummies(combined['dataset'], drop_first=True, dtype=float)
    logc = combined[['log_concentration']].astype(float)

    x_common = sm.add_constant(logc, has_constant='add')
    x_parallel = sm.add_constant(pd.concat([logc, dummies], axis=1), has_constant='add')

    interaction = pd.DataFrame(index=combined.index)
    for column in dummies.columns:
        interaction[f'{column}_x_logc'] = dummies[column] * combined['log_concentration']
    x_full = sm.add_constant(pd.concat([logc, dummies, interaction], axis=1), has_constant='add')

    y = np.column_stack([combined['mortality'].to_numpy(), combined['alive'].to_numpy()])
    family = sm.families.Binomial(link=sm.families.links.Probit())

    full = sm.GLM(y, x_full, family=family).fit()
    parallel = sm.GLM(y, x_parallel, family=family).fit()
    common = sm.GLM(y, x_common, family=family).fit()
    if not (full.converged and parallel.converged and common.converged):
        raise RuntimeError("Global multi-dataset comparison model did not converge.")

    interaction_df = len(ordered_keys) - 1
    shift_df = len(ordered_keys) - 1
    lr_interaction = max(0.0, 2.0 * (full.llf - parallel.llf))
    lr_shift = max(0.0, 2.0 * (parallel.llf - common.llf))

    return {
        'full_result': full,
        'parallel_result': parallel,
        'common_result': common,
        'lr_interaction': float(lr_interaction),
        'interaction_df': int(interaction_df),
        'p_interaction': float(chi2.sf(lr_interaction, interaction_df)),
        'lr_shift': float(lr_shift),
        'shift_df': int(shift_df),
        'p_shift': float(chi2.sf(lr_shift, shift_df)),
    }


def _adjust_pvalue_list(values, method):
    """Adjust a list of p-values while preserving missing values."""
    adjusted = [np.nan] * len(values)
    valid_idx = [i for i, value in enumerate(values) if value is not None and np.isfinite(value)]
    if not valid_idx:
        return adjusted

    raw = [float(values[i]) for i in valid_idx]
    if method == 'None':
        corrected = raw
    else:
        method_name = 'holm' if method == 'Holm' else 'bonferroni'
        corrected = multipletests(raw, alpha=ALPHA_LEVEL, method=method_name)[1]

    for idx, value in zip(valid_idx, corrected):
        adjusted[idx] = float(value)
    return adjusted


def build_multi_comparisons(analyses, selected_keys, mode, reference_key=None,
                            susceptible_reference=False, correction='Holm'):
    """Build reference-vs-all or all-pairwise LC50 comparisons.

    Reference-vs-all comparisons preserve the selected reference as the
    denominator. If that reference is designated susceptible, the directional
    result is labeled a Resistance Ratio.

    All-pairwise comparisons are reported as LC50 fold-differences: the dataset
    with the smaller LC50 is always the denominator, so the reported magnitude
    is always >= 1 and directly expresses the fold separation between datasets.
    """
    selected_keys = sorted(
        selected_keys,
        key=lambda key: dataset_sort_key(analyses[key]['record'])
    )

    if mode == 'Reference vs all':
        if reference_key is None:
            raise ValueError("Select a reference dataset.")
        pairs = [(key, reference_key) for key in selected_keys if key != reference_key]
    else:
        pairs = list(combinations(selected_keys, 2))

    rows = []
    for first_key, second_key in pairs:
        a = analyses[first_key]
        b = analyses[second_key]
        a_name = dataset_display_name(a['record'])
        b_name = dataset_display_name(b['record'])

        ratio = lower = upper = np.nan
        ratio_error = ''
        numerator_key = first_key
        denominator_key = second_key
        numerator_name = a_name
        denominator_name = b_name
        inverted_for_fold = False

        try:
            if mode == 'All pairwise':
                fold = compute_lc50_fold_difference(
                    a['result'], b['result'], name1=a_name, name2=b_name, lc_level=50
                )
                ratio = fold['fold_difference']
                lower = fold['lower']
                upper = fold['upper']
                inverted_for_fold = bool(fold['inverted'])
                if inverted_for_fold:
                    numerator_key, denominator_key = second_key, first_key
                numerator_name = fold['higher_name']
                denominator_name = fold['lower_name']
            else:
                ratio, lower, upper = compute_resistance_ratio_ci(
                    a['result'], b['result'], lc_level=50
                )
        except Exception as exc:
            ratio_error = str(exc)

        try:
            pair_model = compare_dose_response_models(a['df_processed'], b['df_processed'])
            p_slope = float(pair_model['p_slope'])
            p_shift = float(pair_model['p_shift'])
            lr_slope = float(pair_model['lr_slope'])
            lr_shift = float(pair_model['lr_shift'])
        except Exception as exc:
            p_slope = p_shift = lr_slope = lr_shift = np.nan
            ratio_error = (ratio_error + '; ' if ratio_error else '') + f"Curve comparison failed: {exc}"

        is_rr = bool(
            mode == 'Reference vs all'
            and susceptible_reference
            and denominator_key == reference_key
        )

        if mode == 'All pairwise':
            ratio_type = 'LC50 Fold-Difference'
            alphabetical_names = sorted([a_name, b_name], key=str.casefold)
            comparison = f'{alphabetical_names[0]} vs {alphabetical_names[1]}'
            direction = f'{numerator_name} has the higher LC50'
        else:
            ratio_type = 'Resistance Ratio' if is_rr else 'LC50 Ratio'
            comparison = f'{a_name} / {b_name}'
            direction = (
                f'{a_name} has the higher LC50' if np.isfinite(ratio) and ratio > 1
                else f'{b_name} has the higher LC50' if np.isfinite(ratio) and ratio < 1
                else 'LC50 point estimates are equal' if np.isfinite(ratio)
                else ''
            )

        rows.append({
            'numerator_key': numerator_key,
            'denominator_key': denominator_key,
            'comparison': comparison,
            'numerator': numerator_name if mode == 'All pairwise' else a_name,
            'denominator': denominator_name if mode == 'All pairwise' else b_name,
            'ratio_type': ratio_type,
            'ratio': float(ratio) if np.isfinite(ratio) else np.nan,
            'ratio_lower': float(lower) if np.isfinite(lower) else np.nan,
            'ratio_upper': float(upper) if np.isfinite(upper) else np.nan,
            'higher_lc50_dataset': numerator_name if mode == 'All pairwise' else (
                a_name if np.isfinite(ratio) and ratio > 1 else
                b_name if np.isfinite(ratio) and ratio < 1 else ''
            ),
            'lower_lc50_dataset': denominator_name if mode == 'All pairwise' else (
                b_name if np.isfinite(ratio) and ratio > 1 else
                a_name if np.isfinite(ratio) and ratio < 1 else ''
            ),
            'direction': direction,
            'inverted_for_fold_difference': inverted_for_fold,
            'lr_slope': lr_slope,
            'p_slope_raw': p_slope,
            'lr_shift': lr_shift,
            'p_shift_raw': p_shift,
            'note': ratio_error,
            '_sort_a': min(a_name.casefold(), b_name.casefold()),
            '_sort_b': max(a_name.casefold(), b_name.casefold()),
        })

    rows.sort(key=lambda row: (row['_sort_a'], row['_sort_b']))
    slope_adjusted = _adjust_pvalue_list([row['p_slope_raw'] for row in rows], correction)
    shift_adjusted = _adjust_pvalue_list([row['p_shift_raw'] for row in rows], correction)
    for row, p_slope_adj, p_shift_adj in zip(rows, slope_adjusted, shift_adjusted):
        row['p_slope_adjusted'] = p_slope_adj
        row['p_shift_adjusted'] = p_shift_adj
        row['correction'] = correction
        row.pop('_sort_a', None)
        row.pop('_sort_b', None)
    return rows

def create_multi_mortality_plot(analyses):
    """Create observed/fitted mortality curves for multiple datasets."""
    fig, ax = plt.subplots(figsize=(12, 7))
    for _, analysis in sorted_analysis_items(analyses):
        data = analysis['df_processed']
        name = dataset_display_name(analysis['record'])
        conc = np.logspace(np.log10(data['concentration'].min()), np.log10(data['concentration'].max()), 150)
        pred = analysis['result'].predict(sm.add_constant(np.log10(conc))) * 100
        line, = ax.plot(conc, pred, linewidth=2, label=f'{name} fitted')
        ax.scatter(
            data['concentration'],
            data['mortality'] / data['n'] * 100,
            s=55,
            alpha=0.7,
            color=line.get_color(),
            label=f'{name} observed'
        )
    ax.set_xscale('log')
    ax.set_ylim(-5, 105)
    ax.set_xlabel('Concentration')
    ax.set_ylabel('Mortality (%)')
    ax.set_title('Multi-Dataset Mortality Curves')
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8, ncol=2)
    plt.tight_layout()
    return fig


def create_multi_probit_plot(analyses):
    """Create observed/fitted probit lines for multiple datasets."""
    fig, ax = plt.subplots(figsize=(12, 7))
    for _, analysis in sorted_analysis_items(analyses):
        data = analysis['df_processed']
        name = dataset_display_name(analysis['record'])
        mortality_pct = (data['mortality'] / data['n']) * 100
        mortality_adj = mortality_pct.clip(lower=0.1, upper=99.9)
        empirical = norm.ppf(mortality_adj / 100)
        log_range = np.linspace(data['log_concentration'].min(), data['log_concentration'].max(), 150)
        fitted = analysis['intercept'] + analysis['slope'] * log_range
        line, = ax.plot(log_range, fitted, linewidth=2, label=f'{name} fitted')
        ax.scatter(
            data['log_concentration'], empirical,
            s=55, alpha=0.7, color=line.get_color(), label=f'{name} observed'
        )
    ax.set_xlabel('Log10(Concentration)')
    ax.set_ylabel('Probit (Mortality)')
    ax.set_title('Multi-Dataset Probit Regression')
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8, ncol=2)
    plt.tight_layout()
    return fig


def _find_lc_row(analysis, level=50):
    for row in analysis.get('lc_rows', []):
        if np.isclose(row['level'], float(level)):
            return row
    return None


def create_lc50_forest_plot(analyses):
    """Create a forest plot of LC50 estimates and 95% confidence intervals."""
    rows = []
    for _, analysis in sorted_analysis_items(analyses):
        row = _find_lc_row(analysis, 50)
        if row:
            rows.append((dataset_display_name(analysis['record']), row))
    if not rows:
        raise ValueError("No valid LC50 estimates are available for the forest plot.")

    fig_height = max(4.5, 0.55 * len(rows) + 2)
    fig, ax = plt.subplots(figsize=(10, fig_height))
    y = np.arange(len(rows))
    estimates = np.array([row['estimate'] for _, row in rows])
    lower = np.array([row['lower'] for _, row in rows])
    upper = np.array([row['upper'] for _, row in rows])
    xerr = np.vstack([estimates - lower, upper - estimates])
    ax.errorbar(estimates, y, xerr=xerr, fmt='o', capsize=4)
    ax.set_yticks(y)
    ax.set_yticklabels([name for name, _ in rows])
    ax.set_xscale('log')
    ax.set_xlabel('LC50 (log scale)')
    ax.set_title('LC50 Estimates with 95% Confidence Intervals')
    ax.grid(True, axis='x', alpha=0.3)
    ax.invert_yaxis()
    plt.tight_layout()
    return fig


def create_ratio_forest_plot(comparisons, title=None):
    """Create a forest plot for finite resistance ratios or LC50 fold-differences."""
    rows = [row for row in comparisons if np.isfinite(row.get('ratio', np.nan))]
    if not rows:
        raise ValueError("No valid ratios/fold-differences are available for the forest plot.")

    has_fold = any(row.get('ratio_type') == 'LC50 Fold-Difference' for row in rows)
    has_rr = any(row.get('ratio_type') == 'Resistance Ratio' for row in rows)
    if title is None:
        if has_fold and not has_rr:
            title = 'LC50 Fold-Differences'
        elif has_rr and not has_fold:
            title = 'Resistance Ratios'
        else:
            title = 'LC50 Comparisons'

    fig_height = max(4.5, 0.55 * len(rows) + 2)
    fig, ax = plt.subplots(figsize=(10, fig_height))
    y = np.arange(len(rows))
    estimates = np.array([row['ratio'] for row in rows])
    lower = np.array([row['ratio_lower'] for row in rows])
    upper = np.array([row['ratio_upper'] for row in rows])
    xerr = np.vstack([estimates - lower, upper - estimates])
    ax.errorbar(estimates, y, xerr=xerr, fmt='o', capsize=4)
    ax.axvline(1.0, linestyle='--', linewidth=1)
    ax.set_yticks(y)
    ax.set_yticklabels([row['comparison'] for row in rows])
    ax.set_xscale('log')
    if has_fold and not has_rr:
        ax.set_xlabel('LC50 fold-difference (larger LC50 / smaller LC50; log scale)')
    elif has_rr and not has_fold:
        ax.set_xlabel('Resistance ratio (test / susceptible reference; log scale)')
    else:
        ax.set_xlabel('LC50 comparison magnitude (log scale)')
    ax.set_title(title)
    ax.grid(True, axis='x', alpha=0.3)
    ax.invert_yaxis()
    plt.tight_layout()
    return fig

def create_pdf_multi_report(analyses, global_result, comparisons, comparison_mode,
                            correction, lc_levels, chemical='', units='', reference_name=''):
    """Create a consolidated landscape PDF for multi-dataset analysis."""
    pdf = FPDF(orientation='L')
    pdf.set_auto_page_break(auto=True, margin=12)
    pdf.add_page()
    pdf.set_font('Arial', 'B', 16)
    pdf.cell(0, 10, 'Multi-Dataset Probit Analysis Report', 0, 1, 'C')
    pdf.set_font('Arial', '', 10)
    pdf.cell(0, 6, f'Datasets analyzed: {len(analyses)}', 0, 1)
    if chemical:
        pdf.cell(0, 6, f'Chemical/Acaricide: {chemical}', 0, 1)
    if units:
        pdf.cell(0, 6, f'Concentration units: {units}', 0, 1)
    pdf.cell(0, 6, f'Comparison mode: {comparison_mode}', 0, 1)
    if reference_name:
        pdf.cell(0, 6, f'Reference dataset: {reference_name}', 0, 1)
    pdf.cell(0, 6, f'Multiple-testing adjustment: {correction}', 0, 1)
    pdf.cell(0, 6, f'Date: {datetime.now().strftime("%Y-%m-%d %H:%M")}', 0, 1)
    pdf.ln(3)

    pdf.set_font('Arial', 'B', 13)
    pdf.cell(0, 8, 'Global Curve Comparison', 0, 1)
    pdf.set_font('Arial', '', 10)
    pdf.cell(0, 6, f'Dataset x log10(concentration) interaction: LR = {global_result["lr_interaction"]:.4f}, df = {global_result["interaction_df"]}, p = {format_p_value(global_result["p_interaction"])}', 0, 1)
    pdf.cell(0, 6, f'Dataset shift under common slope: LR = {global_result["lr_shift"]:.4f}, df = {global_result["shift_df"]}, p = {format_p_value(global_result["p_shift"])}', 0, 1)
    if global_result['p_interaction'] < ALPHA_LEVEL:
        pdf.multi_cell(0, 5, 'Interpretation: significant slope heterogeneity was detected; common-slope shift tests and single LC50 comparison metrics should be interpreted cautiously.')
    else:
        pdf.multi_cell(0, 5, 'Interpretation: no significant global slope interaction was detected; a common-slope representation is not rejected by this test.')
    pdf.ln(3)

    pdf.set_font('Arial', 'B', 13)
    pdf.cell(0, 8, 'Per-Dataset Model Summary', 0, 1)
    headers = [('Dataset', 52), ('Slope', 24), ('LC50', 32), ('95% CI', 55), ('Dispersion', 28), ('Pearson p', 30), ('LC50 coverage', 35)]
    pdf.set_font('Arial', 'B', 8)
    for label, width in headers:
        pdf.cell(width, 6, label, 1)
    pdf.ln()
    pdf.set_font('Arial', '', 8)
    for _, analysis in sorted_analysis_items(analyses):
        name = dataset_display_name(analysis['record'])[:28]
        lc50 = _find_lc_row(analysis, 50)
        if lc50:
            lc50_text = f'{lc50["estimate"]:.5g}'
            ci_text = f'{lc50["lower"]:.5g}-{lc50["upper"]:.5g}'
            coverage = lc50['status']
        else:
            lc50_text = ci_text = 'Not reported'
            coverage = 'Non-positive slope'
        values = [
            (name, 52),
            (format_slope(analysis["slope"]), 24),
            (lc50_text, 32),
            (ci_text, 55),
            (f'{analysis["dispersion"]:.3f}', 28),
            (format_p_value(analysis["p_value"]), 30),
            (coverage, 35),
        ]
        for value, width in values:
            pdf.cell(width, 6, value, 1)
        pdf.ln()
    pdf.ln(4)

    if comparisons:
        pdf.set_font('Arial', 'B', 13)
        pdf.cell(0, 8, 'Dataset Comparisons', 0, 1)
        headers = [('Comparison', 52), ('Type', 38), ('Fold/Ratio', 22), ('95% CI', 40), ('Higher LC50', 42), ('Slope p adj', 26), ('Shift p adj', 26)]
        pdf.set_font('Arial', 'B', 7)
        for label, width in headers:
            pdf.cell(width, 6, label, 1)
        pdf.ln()
        pdf.set_font('Arial', '', 7)
        for row in comparisons:
            ratio_text = f'{row["ratio"]:.3g}' if np.isfinite(row['ratio']) else 'N/A'
            ci_text = f'{row["ratio_lower"]:.3g}-{row["ratio_upper"]:.3g}' if np.isfinite(row['ratio_lower']) else 'N/A'
            ps_text = format_p_value(row['p_slope_adjusted'])
            ph_text = format_p_value(row['p_shift_adjusted'])
            values = [
                (row['comparison'][:27], 52),
                (row['ratio_type'][:21], 38),
                (ratio_text, 22),
                (ci_text, 40),
                ((row.get('higher_lc50_dataset') or '')[:22], 42),
                (ps_text, 26),
                (ph_text, 26),
            ]
            for value, width in values:
                pdf.cell(width, 6, value, 1)
            pdf.ln()

    # Publication-quality plots on separate pages.
    plot_builders = [
        ('Multi-Dataset Mortality Curves', lambda: create_multi_mortality_plot(analyses)),
        ('Multi-Dataset Probit Regression', lambda: create_multi_probit_plot(analyses)),
    ]
    try:
        plot_builders.append(('LC50 Forest Plot', lambda: create_lc50_forest_plot(analyses)))
    except Exception:
        pass
    if comparisons and any(np.isfinite(row.get('ratio', np.nan)) for row in comparisons):
        plot_builders.append(('LC50 Comparison Forest Plot', lambda: create_ratio_forest_plot(comparisons)))

    for title, builder in plot_builders:
        try:
            fig = builder()
        except Exception:
            continue
        pdf.add_page()
        pdf.set_font('Arial', 'B', 13)
        pdf.cell(0, 8, title, 0, 1)
        with tempfile.NamedTemporaryFile(delete=False, suffix='.png') as tmpfile:
            fig.savefig(tmpfile.name, format='png', dpi=150, bbox_inches='tight')
            plt.close(fig)
            pdf.image(tmpfile.name, x=10, y=25, w=275)
            os.unlink(tmpfile.name)

    return pdf


def _ai_configuration():
    """Return optional AI endpoint configuration supplied by the deployment environment."""
    return {
        'endpoint': os.environ.get('PROBIT_AI_ENDPOINT', '').strip(),
        'api_key': os.environ.get('PROBIT_AI_API_KEY', '').strip(),
        'model': os.environ.get('PROBIT_AI_MODEL', '').strip(),
    }


def _json_safe(value):
    """Recursively replace non-finite numeric values with None for strict JSON payloads."""
    if isinstance(value, dict):
        return {str(k): _json_safe(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_safe(v) for v in value]
    if isinstance(value, (np.floating, float)):
        numeric = float(value)
        return numeric if np.isfinite(numeric) else None
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.bool_,)):
        return bool(value)
    return value


def build_ai_context():
    """Build a structured, raw-data-free context for the optional AI assistant."""
    individual = st.session_state.get('v11_individual_summaries', {})
    multi = st.session_state.get('v11_multi_summary')
    individual_values = sorted(
        individual.values(),
        key=lambda item: (
            str(item.get('strain_population') or item.get('dataset_key', '')).casefold(),
            str(item.get('dataset_key', '')).casefold(),
        )
    )
    return _json_safe({
        'tool_version': '11.7',
        'individual_analyses': individual_values,
        'multi_dataset_analysis': multi,
        'scope_note': (
            'This context contains derived statistics and assay metadata only. Raw uploaded observations are not included.'
        )
    })


def call_ai_assistant(question, prior_messages, context):
    """Call an administrator-configured OpenAI-compatible chat-completions endpoint.

    No third-party Python client is required. The exact endpoint, API key, and
    model are supplied by deployment environment variables. Only structured
    analysis summaries are transmitted, not uploaded raw assay tables.
    """
    config = _ai_configuration()
    if not config['endpoint'] or not config['model']:
        raise RuntimeError("AI assistant endpoint/model is not configured by the application administrator.")

    system_prompt = (
        "You are the explanatory assistant for the Probit Analysis Tool. "
        "Use only the structured analysis context supplied by the application. "
        "Do not recalculate or replace numerical results. Do not invent experimental details. "
        "Clearly distinguish statistical results from biological hypotheses. "
        "Do not claim that slope, parallelism, or resistance ratio alone identifies a molecular mechanism. "
        "If a requested conclusion is unsupported by the supplied context, say so. "
        "When discussing resistance, use the term resistance ratio only when the supplied comparison labels it that way. "
        "For LC50 fold-differences, describe the magnitude and identify which dataset has the higher LC50; do not call the fold-difference a resistance ratio."
    )
    messages = [{'role': 'system', 'content': system_prompt}]
    messages.append({'role': 'system', 'content': 'Structured analysis context:\n' + json.dumps(context, default=str)})
    for message in prior_messages[-10:]:
        if message.get('role') in ('user', 'assistant'):
            messages.append({'role': message['role'], 'content': str(message.get('content', ''))})
    messages.append({'role': 'user', 'content': question})

    payload = {
        'model': config['model'],
        'messages': messages,
        'temperature': 0.2,
    }
    headers = {'Content-Type': 'application/json'}
    if config['api_key']:
        headers['Authorization'] = f"Bearer {config['api_key']}"
    request = urllib.request.Request(
        config['endpoint'],
        data=json.dumps(payload).encode('utf-8'),
        headers=headers,
        method='POST'
    )
    try:
        with urllib.request.urlopen(request, timeout=60) as response:
            body = json.loads(response.read().decode('utf-8'))
    except urllib.error.HTTPError as exc:
        detail = exc.read().decode('utf-8', errors='replace')
        raise RuntimeError(f"AI endpoint returned HTTP {exc.code}: {detail[:500]}") from exc
    except Exception as exc:
        raise RuntimeError(f"AI endpoint request failed: {exc}") from exc

    try:
        return str(body['choices'][0]['message']['content']).strip()
    except Exception as exc:
        raise RuntimeError("AI endpoint response did not match the expected chat-completions schema.") from exc


def _normalize_units_value(value):
    """Normalize concentration-unit labels for comparison without changing display text."""
    text = str(value or '').strip().lower()
    text = text.replace('μ', 'u').replace('µ', 'u')
    text = text.replace('percent', '%')
    text = re.sub(r'\s+', '', text)
    return text


def _compatibility_status(records):
    """Return chemical/unit compatibility information for selected datasets."""
    chemicals = [record.get('chemical', '').strip() for record in records]
    units = [record.get('units', '').strip() for record in records]
    nonblank_chem = {_normalize_column_name(x) for x in chemicals if x}
    nonblank_units = {_normalize_units_value(x) for x in units if x}
    return {
        'chemical_conflict': len(nonblank_chem) > 1,
        'units_conflict': len(nonblank_units) > 1,
        'chemical_missing': any(not x for x in chemicals),
        'units_missing': any(not x for x in units),
        'chemical_display': next((x for x in chemicals if x), ''),
        'units_display': next((x for x in units if x), ''),
    }


def _format_analysis_table(analyses):
    rows = []
    for _, analysis in sorted_analysis_items(analyses):
        lc50 = _find_lc_row(analysis, 50)
        rows.append({
            'Dataset': dataset_display_name(analysis['record']),
            'Slope': format_slope(analysis['slope']),
            'LC50': lc50['estimate'] if lc50 else np.nan,
            'LC50 Lower': lc50['lower'] if lc50 else np.nan,
            'LC50 Upper': lc50['upper'] if lc50 else np.nan,
            'LC50 Coverage': lc50['status'] if lc50 else 'Not reported',
            'Pearson χ²': analysis['chi_square'],
            'df': analysis['df_degrees'],
            'Pearson p': format_p_value(analysis['p_value']),
            'Dispersion φ': analysis['dispersion'],
        })
    return pd.DataFrame(rows)


def _format_comparison_table(comparisons):
    rows = []
    for row in comparisons:
        rows.append({
            'Comparison': row['comparison'],
            'Type': row['ratio_type'],
            'LC50 Fold / Ratio': row['ratio'],
            '95% CI Lower': row['ratio_lower'],
            '95% CI Upper': row['ratio_upper'],
            'Higher LC50 Dataset': row.get('higher_lc50_dataset', ''),
            'Lower LC50 Dataset': row.get('lower_lc50_dataset', ''),
            'Slope p (raw)': format_p_value(row['p_slope_raw']),
            'Slope p (adjusted)': format_p_value(row['p_slope_adjusted']),
            'Shift p (raw)': format_p_value(row['p_shift_raw']),
            'Shift p (adjusted)': format_p_value(row['p_shift_adjusted']),
            'Note': row['note'],
        })
    return pd.DataFrame(rows)


# -----------------------------------------------------------------------------
# Main Version 11 application
# -----------------------------------------------------------------------------

def main():
    st.markdown('<div class="main-header">📊 Probit Analysis Tool</div>', unsafe_allow_html=True)
    st.markdown("**Web Version 11.7** — single- and multi-dataset bioassay probit analysis")

    with st.sidebar:
        st.markdown("### 🎯 Version 11")
        st.markdown(
            "Upload multiple bioassay datasets, analyze individual populations, compare populations to a reference "
            "or pairwise, and optionally enable an administrator-configured AI results assistant."
        )
        st.markdown("---")
        st.markdown("### 📋 Supported input")
        st.markdown(
            "- Legacy files: strain on line 1, chemical on line 2, then table\n"
            "- Header-first CSV/TSV/text files\n"
            "- Required analytical fields: concentration, n tested, mortality\n"
            "- Optional metadata columns: strain/population, chemical/acaricide, units\n"
            "- If units are missing, enter them once for the analysis after upload"
        )
        st.markdown("---")
        st.markdown("### ℹ️ About")
        st.markdown("**Version 11.7 Web**")
        st.markdown("USDA ARS Cattle Fever Tick Research Unit")
        st.markdown("Edinburg, TX, USA")

    tabs = st.tabs([
        "📁 Upload Data",
        "📊 Individual Analyses",
        "⚖️ Multi-Dataset Comparison",
        "💬 AI Assistant",
        "📖 Help",
    ])

    # ------------------------------------------------------------------
    # Upload Data
    # ------------------------------------------------------------------
    with tabs[0]:
        st.markdown("## Upload Bioassay Datasets")
        st.markdown("### ⚙️ Lethal Concentration Estimate Settings")
        lc_levels_input = st.text_input(
            "LC levels to estimate (%)",
            value=st.session_state.get('lc_levels_input', '1, 50, 99'),
            key='lc_levels_input',
            help="Enter percentages separated by commas; values must be >0 and <100."
        )
        try:
            lc_levels = parse_lc_levels(lc_levels_input)
            st.session_state['v11_lc_levels'] = lc_levels
            st.caption("Selected estimates: " + ", ".join(f"LC{format_lc_level(x)}" for x in lc_levels))
        except ValueError as exc:
            st.error(f"Invalid LC settings: {exc}")
            lc_levels = DEFAULT_LC_LEVELS.copy()
            st.session_state['v11_lc_levels'] = lc_levels

        uploaded_files = st.file_uploader(
            "Upload one or more datasets",
            type=['txt', 'tsv', 'csv'],
            accept_multiple_files=True,
            key='v11_multi_uploader',
            help="Each file should represent one population/strain and one chemical concentration scale."
        )

        records = []
        active_keys = []
        seen_signatures = set()
        if uploaded_files:
            st.caption(f"{len(uploaded_files)} file(s) selected")
            for uploaded_file in uploaded_files:
                signature = hashlib.sha256(uploaded_file.getvalue()).hexdigest()[:16]
                if signature in seen_signatures:
                    st.warning(f"Duplicate file content skipped: {uploaded_file.name}")
                    continue
                seen_signatures.add(signature)
                dataset_key = f"v11ds_{signature}"
                active_keys.append(dataset_key)
                try:
                    record = render_uploaded_dataset_dynamic(uploaded_file, dataset_key)
                    records.append(record)
                except Exception as exc:
                    st.error(f"Could not read {uploaded_file.name}: {exc}")

        previous_keys = set(st.session_state.get('v11_active_dataset_keys', []))
        current_keys = set(active_keys)
        for removed_key in previous_keys - current_keys:
            clear_dynamic_dataset_state(removed_key)
        st.session_state['v11_active_dataset_keys'] = active_keys
        records = sort_dataset_records(records)

        # --------------------------------------------------------------
        # Shared concentration-unit configuration
        # --------------------------------------------------------------
        units_confirmed = False
        units_configuration_valid = False
        effective_units = ''
        if records:
            st.markdown("### ⚖️ Concentration Units")
            detected_values = [
                str(record.get('units_detected') or '').strip()
                for record in records
                if str(record.get('units_detected') or '').strip()
            ]
            detected_by_normalized = {}
            for value in detected_values:
                detected_by_normalized.setdefault(_normalize_units_value(value), value)
            detected_unique = list(detected_by_normalized.values())
            units_missing_from_files = any(
                not str(record.get('units_detected') or '').strip() for record in records
            )
            embedded_units_conflict = len(detected_unique) > 1
            detected_common = detected_unique[0] if len(detected_unique) == 1 else ''

            if embedded_units_conflict:
                st.error(
                    "Uploaded files contain conflicting embedded concentration units: "
                    + ", ".join(detected_unique)
                    + ". Correct or separate these datasets before comparison; the shared-unit field will not override embedded units."
                )
                effective_units = ''
            elif units_missing_from_files:
                # Auto-fill from a detected common unit when available, but preserve a manually entered value.
                previous_auto = st.session_state.get('v11_shared_units_auto_value', '')
                current_shared = st.session_state.get('v11_shared_units_input', '')
                if detected_common != previous_auto:
                    if not current_shared or current_shared == previous_auto:
                        st.session_state['v11_shared_units_input'] = detected_common
                    st.session_state['v11_shared_units_auto_value'] = detected_common

                if detected_common:
                    st.info(
                        f"Detected `{detected_common}` in one or more uploaded files. "
                        "The shared entry below will be applied to files that do not contain unit metadata."
                    )
                else:
                    st.info(
                        "No concentration-unit metadata were detected. Enter the units once below; "
                        "the value will be applied to all uploaded datasets."
                    )

                effective_units = st.text_input(
                    "Concentration units for this analysis",
                    key='v11_shared_units_input',
                    help="Examples: ppm, %, µg/mL. This value is applied only to datasets that do not provide units in the uploaded file."
                ).strip()

                if detected_common and effective_units and (
                    _normalize_units_value(effective_units) != _normalize_units_value(detected_common)
                ):
                    st.error(
                        f"The shared unit entry `{effective_units}` does not match the embedded unit `{detected_common}`. "
                        "Use the same concentration scale for all datasets or separate the analyses."
                    )
                elif effective_units:
                    units_configuration_valid = True
                else:
                    st.warning("Enter concentration units before running a multi-dataset comparison.")
            else:
                effective_units = detected_common
                units_configuration_valid = bool(effective_units)
                if effective_units:
                    st.success(f"Concentration units detected in all uploaded datasets: {effective_units}")
                else:
                    st.warning("Concentration units could not be determined from the uploaded datasets.")

            # Apply the one resolved unit value only where metadata were absent. Embedded values remain authoritative.
            if not embedded_units_conflict:
                for record in records:
                    detected = str(record.get('units_detected') or '').strip()
                    if detected:
                        record['units'] = detected
                        record['units_source'] = 'Embedded in file'
                    else:
                        record['units'] = effective_units
                        record['units_source'] = 'Shared analysis entry' if effective_units else 'Not specified'

            # Require an explicit user confirmation for multiple uploaded datasets.
            if len(records) >= 2:
                confirmation_material = '|'.join(
                    [record['signature'] for record in records]
                    + [str(effective_units), ','.join(sorted(detected_unique, key=str.casefold))]
                )
                confirmation_token = hashlib.sha256(confirmation_material.encode('utf-8')).hexdigest()[:12]
                units_confirmed = st.checkbox(
                    (
                        f"I confirm all uploaded datasets use the same concentration units ({effective_units})."
                        if effective_units else
                        "I confirm all uploaded datasets use the same concentration units."
                    ),
                    value=False,
                    key=f'v11_confirm_common_units_{confirmation_token}',
                    disabled=not units_configuration_valid or embedded_units_conflict,
                    help="Required for multi-dataset LC50 ratio/fold-difference comparisons."
                )
            else:
                units_confirmed = bool(units_configuration_valid)

            st.session_state['v11_units_confirmed'] = bool(units_confirmed)
            st.session_state['v11_units_configuration_valid'] = bool(units_configuration_valid and not embedded_units_conflict)
            st.session_state['v11_effective_units'] = effective_units
        else:
            st.session_state['v11_units_confirmed'] = False
            st.session_state['v11_units_configuration_valid'] = False
            st.session_state['v11_effective_units'] = ''

        st.session_state['v11_datasets'] = records

        if records:
            summary_rows = []
            for record in records:
                summary_rows.append({
                    'Dataset': dataset_display_name(record),
                    'File': record['file_name'],
                    'Chemical': record['chemical'] or 'Not specified',
                    'Units': record['units'] or 'Not specified',
                    'Unit Source': record.get('units_source', ''),
                    'Rows': len(record['df']),
                    'Status': 'Ready' if not record['errors'] else 'Needs attention',
                })
            st.markdown("### Dataset Summary")
            st.dataframe(pd.DataFrame(summary_rows), use_container_width=True, hide_index=True)
            ready_count = sum(not record['errors'] for record in records)
            st.info(f"{ready_count} of {len(records)} dataset(s) are ready for analysis.")
        else:
            st.info("Upload one or more files to begin.")

    # Common record state for later tabs.
    records = sort_dataset_records(st.session_state.get('v11_datasets', []))
    valid_records = [record for record in records if not record.get('errors')]
    lc_levels = st.session_state.get('v11_lc_levels', DEFAULT_LC_LEVELS)

    # ------------------------------------------------------------------
    # Individual Analyses
    # ------------------------------------------------------------------
    with tabs[1]:
        st.markdown("## Individual Dataset Analysis")
        if not valid_records:
            st.info("Upload and validate at least one dataset first.")
        else:
            label_to_key = {}
            for i, record in enumerate(valid_records, start=1):
                label = dataset_display_name(record)
                if label in label_to_key:
                    label = f"{label} ({record['file_name']})"
                label_to_key[label] = record['key']
            selected_label = st.selectbox("Dataset", list(label_to_key.keys()), key='v11_individual_dataset')
            selected_key = label_to_key[selected_label]
            selected_record = next(record for record in valid_records if record['key'] == selected_key)

            if st.button("🚀 Run Individual Analysis", type='primary', key='v11_run_individual'):
                try:
                    analysis = fit_dataset_record(selected_record, lc_levels)
                except Exception as exc:
                    st.error(f"Analysis failed: {exc}")
                else:
                    st.session_state.setdefault('v11_individual_summaries', {})[selected_key] = serializable_analysis_summary(analysis)
                    name = dataset_display_name(selected_record)
                    st.markdown(f"### {name}")
                    c1, c2, c3, c4 = st.columns(4)
                    c1.metric("Slope", format_slope(analysis['slope']))
                    c2.metric("Pearson χ²", f"{analysis['chi_square']:.4f}")
                    c3.metric("Pearson p", format_p_value(analysis['p_value']))
                    c4.metric("Dispersion φ", f"{analysis['dispersion']:.3f}")

                    if analysis['slope'] <= 0:
                        st.error(
                            "LC estimates are not reported because the fitted concentration-response slope is non-positive. "
                            "Review concentration coding and assay response before interpreting this dataset."
                        )
                    else:
                        lc_table = pd.DataFrame([{
                            'LC Level': row['label'],
                            'Estimate': row['estimate'],
                            '95% CI Lower': row['lower'],
                            '95% CI Upper': row['upper'],
                            'Concentration Coverage': row['status'],
                            'Target Response Observed': 'Yes' if row['response_covered'] else 'No',
                        } for row in analysis['lc_rows']])
                        st.markdown("### Lethal Concentration Estimates")
                        st.dataframe(lc_table, use_container_width=True, hide_index=True)
                        for row in analysis['lc_rows']:
                            if row['status'] == 'Extrapolated':
                                st.warning(
                                    f"{row['label']} is outside the tested concentration range "
                                    f"({row['concentration_low']:.6g}-{row['concentration_high']:.6g})."
                                )

                    st.markdown("### Model Interpretation")
                    st.info(interpret_slope_biology(analysis['slope'], selected_record.get('chemical') or 'chemical'))
                    if np.isfinite(analysis['r_squared']):
                        st.info(interpret_r_squared_biology(analysis['r_squared']))
                    if analysis['p_value'] < ALPHA_LEVEL:
                        st.warning("Pearson goodness-of-fit indicates evidence of lack of fit and/or extra-binomial variation.")
                    else:
                        st.info("Pearson goodness-of-fit does not detect significant lack of fit.")

                    st.markdown("### Replicate Variability")
                    var_df = calculate_replicate_variability(selected_record['df'][selected_record['df']['concentration'] > 0])
                    st.dataframe(var_df, use_container_width=True, hide_index=True)

                    st.markdown("### Plots")
                    fig1 = create_mortality_plot(analysis['df_processed'], analysis['result'], name, selected_record.get('chemical', ''))
                    st.pyplot(fig1)
                    plt.close(fig1)
                    fig2 = create_probit_plot(analysis['df_processed'], analysis['result'], name, selected_record.get('chemical', ''), lc_levels=lc_levels)
                    st.pyplot(fig2)
                    plt.close(fig2)

                    pdf_options = {
                        'include_data_summary': True,
                        'include_raw_data': True,
                        'include_lc_estimates': True,
                        'include_model_fit': True,
                        'include_parameters': True,
                        'include_mortality_plot': True,
                        'include_probit_plot': True,
                    }
                    pdf = create_pdf_single_with_plots(
                        selected_record['df'], analysis['df_processed'], analysis['result'],
                        analysis['chi_square'], analysis['df_degrees'], analysis['p_value'],
                        name, selected_record.get('chemical', ''), options=pdf_options,
                        lc_levels=lc_levels, units=selected_record.get('units', ''),
                        control_info=selected_record.get('control_info')
                    )
                    st.download_button(
                        "📥 Download Individual PDF Report",
                        data=pdf_to_bytes(pdf),
                        file_name=f"probit_{re.sub(r'[^A-Za-z0-9_.-]+', '_', name)}_{datetime.now().strftime('%Y%m%d_%H%M%S')}.pdf",
                        mime='application/pdf',
                        key='v11_individual_pdf'
                    )

    # ------------------------------------------------------------------
    # Multi-Dataset Comparison
    # ------------------------------------------------------------------
    with tabs[2]:
        st.markdown("## Multi-Dataset Comparison")
        if len(valid_records) < 2:
            st.info("Upload and validate at least two datasets to use multi-dataset comparison.")
        else:
            key_to_label = {}
            label_to_key = {}
            for record in valid_records:
                base = dataset_display_name(record)
                label = base
                counter = 2
                while label in label_to_key:
                    label = f"{base} ({counter})"
                    counter += 1
                key_to_label[record['key']] = label
                label_to_key[label] = record['key']

            selected_labels = st.multiselect(
                "Datasets to compare",
                options=list(label_to_key.keys()),
                default=list(label_to_key.keys()),
                key='v11_compare_selection'
            )
            selected_keys = [label_to_key[label] for label in selected_labels]
            selected_records = sort_dataset_records([record for record in valid_records if record['key'] in selected_keys])

            comparison_mode = st.radio(
                "Comparison mode",
                ['Reference vs all', 'All pairwise'],
                horizontal=True,
                key='v11_comparison_mode'
            )
            correction = st.selectbox(
                "Multiple-comparison adjustment",
                ['Holm', 'Bonferroni', 'None'],
                index=0,
                key='v11_p_adjustment',
                help="Holm controls family-wise error while being less conservative than standard Bonferroni."
            )

            reference_key = None
            susceptible_reference = False
            if comparison_mode == 'Reference vs all' and selected_labels:
                reference_label = st.selectbox("Reference dataset", selected_labels, key='v11_reference_dataset')
                reference_key = label_to_key[reference_label]
                susceptible_reference = st.checkbox(
                    "This is a susceptible reference population",
                    value=False,
                    key='v11_susceptible_reference',
                    help="When selected, reference-based LC50 ratios are labeled Resistance Ratios."
                )

            compatibility = _compatibility_status(selected_records)
            compatible = True
            if compatibility['chemical_conflict']:
                st.error("Selected datasets contain different chemical/acaricide names. Analyze each chemical separately.")
                compatible = False
            if compatibility['units_conflict']:
                st.error("Selected datasets contain different concentration units. Convert them to a common scale before comparison.")
                compatible = False

            if selected_records and compatibility['chemical_missing']:
                confirmed = st.checkbox(
                    "I confirm all selected datasets use the same chemical/active ingredient",
                    value=False,
                    key='v11_confirm_chemical'
                )
                compatible = compatible and confirmed

            if selected_records and compatibility['units_missing']:
                st.error(
                    "Concentration units are unresolved for one or more selected datasets. "
                    "Enter the shared units in the Upload Data tab before comparison."
                )
                compatible = False
            elif selected_records:
                if not st.session_state.get('v11_units_confirmed', False):
                    st.warning(
                        "Confirm in the Upload Data tab that all uploaded datasets use the same concentration units "
                        "before running a multi-dataset comparison."
                    )
                    compatible = False
                else:
                    st.success(f"Common concentration units confirmed: {compatibility['units_display']}")

            if len(selected_records) < 2:
                compatible = False
                st.warning("Select at least two datasets.")

            if st.button(
                "🚀 Run Multi-Dataset Analysis",
                type='primary',
                disabled=not compatible,
                key='v11_run_multi'
            ):
                analyses, failures = run_dataset_analyses(selected_records, lc_levels)
                if failures:
                    for key, message in failures.items():
                        label = key_to_label.get(key, key)
                        st.error(f"{label}: {message}")
                if len(analyses) < 2:
                    st.error("Fewer than two datasets could be fitted successfully; multi-dataset analysis cannot continue.")
                else:
                    selected_success_keys = [
                        key for key, _ in sorted_analysis_items(
                            {key: analyses[key] for key in selected_keys if key in analyses}
                        )
                    ]
                    if comparison_mode == 'Reference vs all' and reference_key not in analyses:
                        st.error("The selected reference dataset did not fit successfully.")
                    else:
                        try:
                            global_result = compare_multi_dose_response_models({key: analyses[key] for key in selected_success_keys})
                            comparisons = build_multi_comparisons(
                                analyses, selected_success_keys, comparison_mode,
                                reference_key=reference_key,
                                susceptible_reference=susceptible_reference,
                                correction=correction
                            )
                        except Exception as exc:
                            st.error(f"Multi-dataset comparison failed: {exc}")
                        else:
                            st.markdown("### Global Curve Comparison")
                            g1, g2 = st.columns(2)
                            with g1:
                                st.metric("Slope interaction LR", f"{global_result['lr_interaction']:.4f}")
                                st.metric("Slope interaction p", format_p_value(global_result['p_interaction']))
                                st.caption(f"df = {global_result['interaction_df']}")
                            with g2:
                                st.metric("Dataset shift LR", f"{global_result['lr_shift']:.4f}")
                                st.metric("Dataset shift p", format_p_value(global_result['p_shift']))
                                st.caption(f"df = {global_result['shift_df']}")

                            if global_result['p_interaction'] < ALPHA_LEVEL:
                                st.warning(
                                    "A significant global dataset × log10(concentration) interaction was detected. "
                                    "The dose-response slopes are heterogeneous, so common-slope shifts and single LC50 comparison metrics "
                                    "do not fully summarize the differences among populations."
                                )
                            else:
                                st.info(
                                    "No significant global slope interaction was detected. A common-slope representation is "
                                    "not rejected by this test."
                                )

                            st.markdown("### Per-Dataset Summary")
                            analysis_table = _format_analysis_table({key: analyses[key] for key in selected_success_keys})
                            st.dataframe(analysis_table, use_container_width=True, hide_index=True)

                            st.markdown("### Dataset Comparisons")
                            comparison_table = _format_comparison_table(comparisons)
                            st.dataframe(comparison_table, use_container_width=True, hide_index=True)
                            st.caption(
                                f"Pairwise slope and common-slope shift p-values are adjusted using {correction}. "
                                "Resistance-ratio and LC50 fold-difference confidence intervals remain individual 95% confidence intervals."
                            )

                            st.markdown("### Multi-Dataset Plots")
                            fig = create_multi_mortality_plot({key: analyses[key] for key in selected_success_keys})
                            st.pyplot(fig)
                            plt.close(fig)
                            fig = create_multi_probit_plot({key: analyses[key] for key in selected_success_keys})
                            st.pyplot(fig)
                            plt.close(fig)
                            try:
                                fig = create_lc50_forest_plot({key: analyses[key] for key in selected_success_keys})
                                st.pyplot(fig)
                                plt.close(fig)
                            except Exception as exc:
                                st.info(f"LC50 forest plot unavailable: {exc}")
                            try:
                                if comparison_mode == 'All pairwise':
                                    ratio_title = 'LC50 Fold-Differences'
                                elif any(row['ratio_type'] == 'Resistance Ratio' for row in comparisons):
                                    ratio_title = 'Resistance Ratios'
                                else:
                                    ratio_title = 'LC50 Ratios'
                                fig = create_ratio_forest_plot(comparisons, title=ratio_title)
                                st.pyplot(fig)
                                plt.close(fig)
                            except Exception as exc:
                                st.info(f"Ratio forest plot unavailable: {exc}")

                            reference_name = ''
                            if reference_key and reference_key in analyses:
                                reference_name = dataset_display_name(analyses[reference_key]['record'])
                            pdf = create_pdf_multi_report(
                                {key: analyses[key] for key in selected_success_keys},
                                global_result,
                                comparisons,
                                comparison_mode,
                                correction,
                                lc_levels,
                                chemical=compatibility['chemical_display'],
                                units=compatibility['units_display'],
                                reference_name=reference_name,
                            )
                            st.download_button(
                                "📥 Download Multi-Dataset PDF Report",
                                data=pdf_to_bytes(pdf),
                                file_name=f"probit_multi_dataset_{datetime.now().strftime('%Y%m%d_%H%M%S')}.pdf",
                                mime='application/pdf',
                                key='v11_multi_pdf'
                            )

                            individual_summaries = {
                                key: serializable_analysis_summary(analyses[key]) for key in selected_success_keys
                            }
                            st.session_state.setdefault('v11_individual_summaries', {}).update(individual_summaries)
                            st.session_state['v11_multi_summary'] = {
                                'comparison_mode': comparison_mode,
                                'reference_dataset': reference_name,
                                'susceptible_reference': bool(susceptible_reference),
                                'multiple_testing_adjustment': correction,
                                'chemical_acaricide': compatibility['chemical_display'],
                                'concentration_units': compatibility['units_display'],
                                'global_slope_interaction': {
                                    'lr': global_result['lr_interaction'],
                                    'df': global_result['interaction_df'],
                                    'p_value': global_result['p_interaction'],
                                },
                                'global_dataset_shift': {
                                    'lr': global_result['lr_shift'],
                                    'df': global_result['shift_df'],
                                    'p_value': global_result['p_shift'],
                                    'interpret_with_caution_if_slope_interaction_significant': True,
                                },
                                'comparisons': comparisons,
                            }
                            st.success("✓ Multi-dataset analysis complete")

    # ------------------------------------------------------------------
    # Optional AI Assistant
    # ------------------------------------------------------------------
    with tabs[3]:
        st.markdown("## AI Results Assistant")
        st.info(
            "The assistant is optional and disabled unless the deployment administrator configures an approved endpoint. "
            "When enabled, the application sends structured analysis summaries and metadata—not the uploaded raw assay table—to the endpoint."
        )
        context = build_ai_context()
        has_results = bool(context['individual_analyses'] or context['multi_dataset_analysis'])
        if not has_results:
            st.warning("Run at least one individual or multi-dataset analysis before asking questions about results.")

        config = _ai_configuration()
        configured = bool(config['endpoint'] and config['model'])
        if not configured:
            st.markdown("### Assistant not configured")
            st.markdown(
                "An administrator can enable the chat by setting these server-side environment variables:\n\n"
                "- `PROBIT_AI_ENDPOINT` — full approved OpenAI-compatible chat-completions endpoint\n"
                "- `PROBIT_AI_MODEL` — approved model/deployment name\n"
                "- `PROBIT_AI_API_KEY` — optional bearer token/API key\n\n"
                "No AI credential is stored in the source code or repository."
            )
        elif not has_results:
            st.caption("The endpoint is configured; analysis results are needed before chat can begin.")
        else:
            st.caption(f"Configured model/deployment: {config['model']}")
            st.caption("Questions should focus on interpretation of the statistics supplied by the application.")
            history = st.session_state.setdefault('v11_ai_messages', [])
            for message in history:
                with st.chat_message(message['role']):
                    st.markdown(message['content'])

            prompt = st.chat_input("Ask about the analysis results")
            if prompt:
                with st.chat_message('user'):
                    st.markdown(prompt)
                prior = list(history)
                history.append({'role': 'user', 'content': prompt})
                with st.chat_message('assistant'):
                    with st.spinner("Interpreting the analysis..."):
                        try:
                            answer = call_ai_assistant(prompt, prior, context)
                        except Exception as exc:
                            answer = f"The AI assistant could not complete the request: {exc}"
                    st.markdown(answer)
                history.append({'role': 'assistant', 'content': answer})

    # ------------------------------------------------------------------
    # Help
    # ------------------------------------------------------------------
    with tabs[4]:
        st.markdown("## Help & Documentation")
        st.markdown(
            """
### Version 11 workflow

1. **Upload Data** — Upload one or many CSV/TSV/text bioassay files. Metadata are auto-detected when possible and remain editable.
2. **Individual Analyses** — Select any validated dataset and run the same grouped-binomial probit analysis used in Version 10.7. Dataset lists, tables, plots, and reports are organized alphabetically by strain/population name.
3. **Multi-Dataset Comparison** — Compare all selected populations to one reference or run all pairwise comparisons.
4. **AI Assistant** — Optional. Available only when an administrator configures an approved AI endpoint.

### Multi-dataset statistics

Each dataset is fitted independently with:

```
Probit(mortality) = β₀ + β₁ × log10(concentration)
```

The global multi-population model is then fitted using dataset indicators and dataset × log10(concentration) interactions.

- **Global slope-interaction test:** likelihood-ratio test of the interaction terms.
- **Global dataset-shift test:** likelihood-ratio test of dataset effects under a common-slope model.
- **Reference vs all:** each non-reference population is compared with the selected reference.
- **All pairwise:** every selected dataset pair is compared using an **LC50 fold-difference**, with the smaller LC50 always placed in the denominator so the reported fold-difference is ≥1.
- **Multiplicity:** Holm adjustment is recommended; Bonferroni or no adjustment are also available.
- **Reference ratios:** LC50 ratios are labeled **Resistance Ratios** only when the denominator is explicitly designated as a susceptible reference population. Otherwise, reference-based comparisons preserve the selected reference as the denominator and are labeled LC50 Ratios.
- **Pairwise fold-differences:** all-pairwise comparisons report larger LC50 / smaller LC50. If the original directional ratio is <1, both the ratio and CI are inverted: fold = 1/ratio and CI = (1/upper, 1/lower). The point estimate is therefore always ≥1; its CI may still include values below 1 when the direction of the difference is statistically uncertain.
- **Confidence intervals:** resistance-ratio and LC50 fold-difference confidence intervals remain individual delta-method 95% CIs on the log10 ratio scale; p-value adjustment does not convert them into simultaneous confidence intervals.

### Interpretation safeguards

- LC estimates are suppressed when a fitted slope is non-positive.
- LC estimates outside the tested concentration range are labeled extrapolated.
- Different chemicals or incompatible concentration units cannot be combined in one ratio analysis.
- If concentration-unit metadata are missing, units are entered once in the Upload Data tab and applied to datasets lacking embedded units.
- Multi-dataset analysis requires the user to confirm that all uploaded datasets use the same concentration units.
- Missing chemical metadata still require explicit user confirmation before comparison.
- Significant slope interaction means a single LC50 ratio or fold-difference does not fully summarize curve differences.
- Pearson dispersion φ = Pearson χ² / residual df is reported as an indicator of extra-binomial variation.

### Optional AI assistant

The statistical engine remains authoritative. The optional assistant receives a structured summary of results rather than the uploaded raw table and is instructed not to recalculate or replace numerical results. The feature remains disabled unless an approved endpoint is configured outside the repository through environment variables.

### Security and deployment

- Uploaded assay files are processed by the Streamlit server hosting this application.
- The application does not intentionally persist uploaded raw assay files to permanent storage.
- AI credentials are read from server-side environment variables and are not embedded in source code.
- If AI is enabled, deployment administrators are responsible for ensuring the endpoint, data handling, retention, and authorization comply with applicable organizational requirements.

### References

- Finney, D.J. (1971). *Probit Analysis*, 3rd ed. Cambridge University Press.
- Robertson, J.L., Russell, R.M., Preisler, H.K., and Savin, N.E. (2007). *Bioassays with Arthropods*, 2nd ed. CRC Press.
- Ver Hoef, J.M. (2012). Who invented the delta method? *The American Statistician* 66:124-127.
- Abbott, W.S. (1925). A method of computing the effectiveness of an insecticide. *Journal of Economic Entomology* 18:265-267.
            """
        )


if __name__ == "__main__":
    main()
