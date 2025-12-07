"""

@author: Carolin Ector
Date: August 2025

"""

"""
Growth Analysis Script with Error Propagation
Clock-Neuroblastoma Figure S4
"""

import os
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.signal import savgol_filter
from openpyxl import Workbook
from openpyxl.utils.dataframe import dataframe_to_rows

# --- CONFIGURATION ---
data_folder = 'growth_data'
celllines = ["CHP212", "SKNAS", "SY5Y", "NGP", "GIMEN", "SKNSH", "CLBGA", "IMR5", "SKNBE", "SKNBE2"]
output_dir = 'growth_analysis_results'
os.makedirs(output_dir, exist_ok=True)

# Smoothing parameters
RAW_SMOOTH_WINDOW = 5
RAW_SMOOTH_POLY = 2
NORM_SMOOTH_WINDOW = 21
NORM_SMOOTH_POLY = 3

def logistic(t, L, k, t0):
    """Logistic growth function"""
    return L / (1 + np.exp(-k * (t - t0)))

def calculate_doubling_time(k):
    """Calculate doubling time from growth rate k"""
    return np.log(2) / k if k > 0 else np.inf

def calculate_doubling_time_error(k, k_err):
    """Calculate doubling time error from growth rate k and its error"""
    if k > 0:
        dt_err = np.log(2) * k_err / (k**2)
        return dt_err
    else:
        return np.inf

def fit_logistic_curve_with_errors(t, y, y_err=None):
    """Fit logistic curve with error propagation and return parameters with uncertainties"""
    try:
        L_init = np.max(y) * 1.1
        k_init = 0.05
        t0_init = np.median(t)
        
        bounds = ([0, 0, 0], [np.inf, 1, np.max(t)])
        
        if y_err is not None and np.any(y_err > 0):
            weights = 1.0 / np.maximum(y_err, np.min(y_err[y_err > 0]) * 0.1)
            popt, pcov = curve_fit(logistic, t, y, p0=[L_init, k_init, t0_init], 
                                 bounds=bounds, maxfev=10000, sigma=y_err, absolute_sigma=True)
        else:
            popt, pcov = curve_fit(logistic, t, y, p0=[L_init, k_init, t0_init], 
                                 bounds=bounds, maxfev=10000)
        
        param_errors = np.sqrt(np.diag(pcov))
        
        fitted = logistic(t, *popt)
        ss_res = np.sum((y - fitted) ** 2)
        ss_tot = np.sum((y - np.mean(y)) ** 2)
        r2 = 1 - (ss_res / ss_tot) if ss_tot > 0 else 0
        
        return popt, param_errors, r2, fitted, pcov
    except Exception as e:
        print(f"    Curve fitting failed: {e}")
        return None, None, np.nan, None, None

def bootstrap_fit_analysis(data_list, n_bootstrap=1000):
    """Perform bootstrap analysis on individual curves to estimate parameter uncertainties"""
    if len(data_list) < 3:
        return None, None
    
    bootstrap_params = []
    
    max_time = 92
    common_time = np.linspace(0, max_time, int(max_time * 2) + 1)
    
    aligned_data = []
    for d in data_list:
        t_orig = np.array(d['time'])
        conf_orig = np.array(d['confluence'])
        
        valid_mask = t_orig <= max_time
        t_valid = t_orig[valid_mask]
        conf_valid = conf_orig[valid_mask]
        
        if len(t_valid) < 5:
            continue
        
        conf_interp = np.interp(common_time, t_valid, conf_valid)
        
        if len(conf_interp) >= NORM_SMOOTH_WINDOW:
            conf_interp = savgol_filter(conf_interp, NORM_SMOOTH_WINDOW, NORM_SMOOTH_POLY)
        
        if conf_interp[0] != 0:
            conf_interp = conf_interp / conf_interp[0]
        
        aligned_data.append(conf_interp)
    
    if len(aligned_data) < 3:
        return None, None
    
    aligned_array = np.array(aligned_data)
    
    for _ in range(n_bootstrap):
        indices = np.random.choice(len(aligned_data), size=len(aligned_data), replace=True)
        bootstrap_sample = aligned_array[indices]
        
        bootstrap_mean = np.mean(bootstrap_sample, axis=0)
        
        if len(bootstrap_mean) >= NORM_SMOOTH_WINDOW:
            bootstrap_mean = savgol_filter(bootstrap_mean, NORM_SMOOTH_WINDOW, NORM_SMOOTH_POLY)
            bootstrap_mean = savgol_filter(bootstrap_mean, NORM_SMOOTH_WINDOW, NORM_SMOOTH_POLY)
        
        params, _, r2, _, _ = fit_logistic_curve_with_errors(common_time, bootstrap_mean)
        if params is not None and r2 > 0.5:
            bootstrap_params.append(params)
    
    if len(bootstrap_params) < 10:
        return None, None
    
    bootstrap_params = np.array(bootstrap_params)
    
    mean_params = np.mean(bootstrap_params, axis=0)
    std_params = np.std(bootstrap_params, axis=0)
    
    return mean_params, std_params

def find_normalization_point(t, threshold_hours=22, jump_threshold=4, max_search_time=26):
    """Find x0 normalization point"""
    for i, time_point in enumerate(t):
        if abs(time_point - threshold_hours) <= 0.5:
            return i
    
    for i in range(1, len(t)):
        if t[i] > max_search_time:
            break
        if t[i] - t[i-1] >= jump_threshold:
            return i
    
    return np.argmin(np.abs(t - threshold_hours))

def process_data_files(cellline):
    """Process Excel file for a given cell line"""
    print(f"  Processing cell line: {cellline}")
    
    file_pattern = os.path.join(data_folder, f"Growth_{cellline}.xlsx")
    files = glob.glob(file_pattern)
    
    if not files:
        print(f"    No file found for {cellline}")
        return []
    
    file_path = files[0]
    print(f"    Processing file: {file_path}")
    
    try:
        df = pd.read_excel(file_path)
    except Exception as e:
        print(f"    Error reading {file_path}: {e}")
        return []
    
    # Assume first column is time
    time_col = df.columns[0]
    t = np.array(df[time_col].values)
    
    data_list = []
    
    # Process each data column
    data_columns = df.columns[1:]  # All columns except first (time)
    
    for i, col in enumerate(data_columns):
        vals = pd.to_numeric(df[col], errors='coerce')
        
        valid_mask = ~np.isnan(vals)
        t_clean = np.array(t[valid_mask])
        vals_clean = np.array(vals[valid_mask])
        
        if len(vals_clean) < 10:
            print(f"    Skipping column {col} - insufficient data points")
            continue
        
        data_list.append({
            'time': t_clean,
            'confluence': vals_clean,
            'image_id': f"Image_{i+1}"
        })
    
    print(f"    Loaded {len(data_list)} valid image datasets")
    return data_list

def calculate_averages_with_x0_stats(data_list, smooth_normalized=False):
    """Calculate mean and std from list of time series, including x0 and raw max confluence statistics"""
    if not data_list:
        return None, None, None, None, None, None, None, None
    
    print(f"    Averaging {len(data_list)} series")
    
    # Extract x0 statistics
    x0_values = [d.get('x0_confluence', 0) for d in data_list if 'x0_confluence' in d]
    x0_times = [d.get('x0_time', 0) for d in data_list if 'x0_time' in d]
    
    x0_confluence_mean = np.mean(x0_values) if x0_values else 0
    x0_confluence_std = np.std(x0_values, ddof=1) if len(x0_values) > 1 else 0
    x0_time_mean = np.mean(x0_times) if x0_times else 0
    x0_time_std = np.std(x0_times, ddof=1) if len(x0_times) > 1 else 0
    
    # Extract raw maximum confluence at 92 hours from raw data
    raw_max_confluence_values = [d.get('raw_max_confluence', 0) for d in data_list if 'raw_max_confluence' in d]
    
    raw_max_confluence_mean = np.mean(raw_max_confluence_values) if raw_max_confluence_values else 0
    raw_max_confluence_std = np.std(raw_max_confluence_values, ddof=1) if len(raw_max_confluence_values) > 1 else 0
    
    if smooth_normalized:
        max_time = 92
        common_time = np.linspace(0, max_time, int(max_time * 2) + 1)
        
        aligned_data = []
        
        for d in data_list:
            t_orig = np.array(d['time'])
            conf_orig = np.array(d['confluence'])
            
            valid_mask = t_orig <= max_time
            t_valid = t_orig[valid_mask]
            conf_valid = conf_orig[valid_mask]
            
            if len(t_valid) < 5:
                continue
            
            conf_interp = np.interp(common_time, t_valid, conf_valid)
            
            if len(conf_interp) >= NORM_SMOOTH_WINDOW:
                conf_interp = savgol_filter(conf_interp, NORM_SMOOTH_WINDOW, NORM_SMOOTH_POLY)
            
            aligned_data.append(conf_interp)
        
        if not aligned_data:
            return None, None, None, x0_confluence_mean, x0_confluence_std, x0_time_mean, raw_max_confluence_mean, raw_max_confluence_std
        
        aligned_array = np.array(aligned_data)
        mean_conf = np.mean(aligned_array, axis=0)
        std_conf = np.std(aligned_array, axis=0, ddof=1)
        
        if mean_conf[0] != 0:
            mean_conf = mean_conf / mean_conf[0]
            std_conf = std_conf / mean_conf[0] if mean_conf[0] != 0 else std_conf
        
        if len(mean_conf) >= NORM_SMOOTH_WINDOW:
            mean_conf = savgol_filter(mean_conf, NORM_SMOOTH_WINDOW, NORM_SMOOTH_POLY)
            mean_conf = savgol_filter(mean_conf, NORM_SMOOTH_WINDOW, NORM_SMOOTH_POLY)
            std_conf = savgol_filter(std_conf, NORM_SMOOTH_WINDOW, NORM_SMOOTH_POLY)
        
        return common_time, mean_conf, std_conf, x0_confluence_mean, x0_confluence_std, x0_time_mean, raw_max_confluence_mean, raw_max_confluence_std
    
    else:
        min_length = min(len(d['time']) for d in data_list)
        
        aligned_data = []
        common_time = None
        
        for d in data_list:
            t_truncated = np.array(d['time'][:min_length])
            conf_truncated = np.array(d['confluence'][:min_length])
            aligned_data.append(conf_truncated)
            if common_time is None:
                common_time = t_truncated
        
        if not aligned_data:
            return None, None, None, x0_confluence_mean, x0_confluence_std, x0_time_mean, raw_max_confluence_mean, raw_max_confluence_std
        
        aligned_array = np.array(aligned_data)
        mean_conf = np.mean(aligned_array, axis=0)
        std_conf = np.std(aligned_array, axis=0, ddof=1)
        
        return common_time, mean_conf, std_conf, x0_confluence_mean, x0_confluence_std, x0_time_mean, raw_max_confluence_mean, raw_max_confluence_std

def normalize_and_truncate_data(data_list):
    """Normalize data to x0 and limit to 92 hours post-x0, while preserving raw max confluence"""
    normalized_data = []
    
    for data_entry in data_list:
        t = np.array(data_entry['time'])
        conf = np.array(data_entry['confluence'])
        
        if len(t) < 5 or len(conf) < 5:
            continue
        
        norm_idx = find_normalization_point(t)
        
        if norm_idx >= len(t) or norm_idx >= len(conf):
            continue
        
        x0_confluence = conf[norm_idx]
        x0_time = t[norm_idx]
        
        # Calculate max confluence from raw data at 92 hours post-x0
        t_post_x0 = t[norm_idx:] - t[norm_idx]  # Time relative to x0
        conf_post_x0 = conf[norm_idx:]  # Raw confluence post-x0
        
        # Find confluence at 92 hours post-x0 from raw data
        max_time = 92
        if len(t_post_x0) > 0 and np.max(t_post_x0) >= 85:  # If we have data near 92 hours
            idx_92 = np.argmin(np.abs(t_post_x0 - max_time))
            raw_max_confluence = conf_post_x0[idx_92]
        elif len(conf_post_x0) > 0:  # Take the last available point
            raw_max_confluence = conf_post_x0[-1]
        else:
            raw_max_confluence = x0_confluence
        
        # Data normalization for growth curve fitting
        t_norm = t[norm_idx:]
        conf_norm = conf[norm_idx:]
        
        time_mask = (t_norm - t_norm[0]) <= max_time
        t_norm = t_norm[time_mask]
        conf_norm = conf_norm[time_mask]
        
        if len(t_norm) == 0 or len(conf_norm) == 0:
            continue
        
        t_norm = t_norm - t_norm[0]
        
        if x0_confluence > 0:
            normalized_conf = conf_norm / x0_confluence
        else:
            normalized_conf = conf_norm
            x0_confluence = 1.0
        
        normalized_data.append({
            'time': t_norm,
            'confluence': normalized_conf,
            'x0_confluence': x0_confluence,
            'x0_time': x0_time,
            'raw_max_confluence': raw_max_confluence,  # Raw confluence at 92h
            'image_id': data_entry['image_id']
        })
    
    return normalized_data

def create_fitted_plots_and_save_params(cellline, normalized_data):
    """Create fitted plots and save parameters to Excel with error propagation and raw max confluence"""
    fig, ax = plt.subplots(1, 1, figsize=(10, 8))
    fig.suptitle(f'{cellline} - Fitted Growth Curve with Error Bars', fontsize=16)
    
    fit_results = []
    
    if normalized_data:
        t, mean_conf, std_conf, x0_mean, x0_std, x0_time_mean, raw_max_mean, raw_max_std = calculate_averages_with_x0_stats(normalized_data, smooth_normalized=True)
        
        if t is not None:
            # Calculate x0_time_std for this group
            x0_times = [d.get('x0_time', 0) for d in normalized_data if 'x0_time' in d]
            x0_time_std = np.std(x0_times, ddof=1) if len(x0_times) > 1 else 0
            
            ax.plot(t, mean_conf, color='blue', 
                   label=f'{cellline} (x0={x0_mean:.1f}±{x0_std:.1f}%)', linewidth=2)
            ax.fill_between(t, mean_conf - std_conf, mean_conf + std_conf, 
                           color='blue', alpha=0.3)
            
            params, param_errors, r2, fitted, pcov = fit_logistic_curve_with_errors(t, mean_conf, std_conf)
            bootstrap_params, bootstrap_errors = bootstrap_fit_analysis(normalized_data)
            
            if params is not None:
                L, k, t0 = params
                L_err, k_err, t0_err = param_errors if param_errors is not None else (0, 0, 0)
                
                if bootstrap_errors is not None:
                    k_err = max(k_err, bootstrap_errors[1])
                
                doubling_time = calculate_doubling_time(k)
                doubling_time_err = calculate_doubling_time_error(k, k_err)
                
                ax.plot(t, fitted, '--', color='red', 
                       alpha=0.8, linewidth=1.5, label='Logistic Fit')
                
                fit_results.append({
                    'cellline': cellline,
                    'L': L,
                    'L_err': L_err,
                    'k': k,
                    'k_err': k_err,
                    't0': t0,
                    't0_err': t0_err,
                    'doubling_time': doubling_time,
                    'doubling_time_err': doubling_time_err,
                    'r2': r2,
                    'x0_confluence_mean': x0_mean,
                    'x0_confluence_std': x0_std,
                    'x0_time_mean': x0_time_mean,
                    'x0_time_std': x0_time_std,
                    'raw_max_confluence_mean': raw_max_mean,
                    'raw_max_confluence_std': raw_max_std,
                    'n_images': len(normalized_data),
                    'bootstrap_available': bootstrap_errors is not None
                })
                
                param_text = f'k={k:.3f}±{k_err:.3f}, DT={doubling_time:.1f}±{doubling_time_err:.1f}h, R²={r2:.2f}'
                ax.text(0.02, 0.98, param_text, transform=ax.transAxes, fontsize=10,
                       verticalalignment='top', bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8))
    
    ax.set_xlabel('Time (hours)')
    ax.set_ylabel('Normalized Confluence')
    ax.set_title(f'{cellline} Growth Analysis')
    ax.grid(True, alpha=0.3)
    ax.legend()
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f'{cellline}_fitted_growth_curve.svg'), 
                format='svg', bbox_inches='tight')
    plt.close()
    
    return fit_results

def save_fit_parameters_to_excel(all_fit_results):
    """Save all fitting parameters with error estimates, x0 confluence, and raw max confluence data to Excel files"""
    if not all_fit_results:
        print("No fitting results to save")
        return
    
    df = pd.DataFrame(all_fit_results)
    
    # Enhanced column order including raw max confluence data
    column_order = [
        'cellline',
        'x0_confluence_mean', 'x0_confluence_std', 'x0_time_mean', 'x0_time_std',
        'raw_max_confluence_mean', 'raw_max_confluence_std',
        'k', 'k_err', 'doubling_time', 'doubling_time_err',
        'L', 'L_err', 't0', 't0_err', 'r2',
        'n_images', 'bootstrap_available'
    ]
    df = df[column_order]
    
    # Calculate additional metrics using raw confluence values
    df['x0_confluence_cv'] = (df['x0_confluence_std'] / df['x0_confluence_mean'] * 100).round(2)
    df['raw_max_confluence_cv'] = (df['raw_max_confluence_std'] / df['raw_max_confluence_mean'] * 100).round(2)
    df['raw_growth_fold_change'] = (df['raw_max_confluence_mean'] / df['x0_confluence_mean']).round(2)
    
    excel_path = os.path.join(output_dir, 'growth_analysis_parameters.xlsx')
    
    with pd.ExcelWriter(excel_path, engine='openpyxl') as writer:
        # All parameters sheet
        df.to_excel(writer, sheet_name='All_Parameters', index=False)
        
        # Growth summary
        growth_summary = df[['cellline', 
                           'x0_confluence_mean', 'x0_confluence_std', 'x0_confluence_cv',
                           'raw_max_confluence_mean', 'raw_max_confluence_std', 'raw_max_confluence_cv',
                           'raw_growth_fold_change',
                           'k', 'k_err', 'doubling_time', 'doubling_time_err', 'r2', 'n_images']].copy()
        growth_summary.to_excel(writer, sheet_name='Growth_Summary', index=False)
        
        # Confluence analysis
        confluence_summary = df[['cellline', 
                               'x0_confluence_mean', 'x0_confluence_std', 'x0_confluence_cv',
                               'raw_max_confluence_mean', 'raw_max_confluence_std', 'raw_max_confluence_cv',
                               'raw_growth_fold_change', 'x0_time_mean', 'x0_time_std',
                               'n_images']].copy()
        confluence_summary.to_excel(writer, sheet_name='Confluence_Analysis', index=False)
        
        # Growth kinetics
        kinetics_summary = df[['cellline',
                             'k', 'k_err', 'doubling_time', 'doubling_time_err',
                             'L', 'L_err', 't0', 't0_err', 'r2']].copy()
        kinetics_summary.to_excel(writer, sheet_name='Growth_Kinetics', index=False)
    
    print(f"Saved growth analysis parameters to {excel_path}")
    
    # Focused analysis file
    focused_path = os.path.join(output_dir, 'growth_analysis_focused.xlsx')
    focused_analysis = df[['cellline', 
                          'x0_confluence_mean', 'x0_confluence_std', 'x0_confluence_cv',
                          'raw_max_confluence_mean', 'raw_max_confluence_std', 'raw_max_confluence_cv',
                          'raw_growth_fold_change',
                          'k', 'k_err', 'doubling_time', 'doubling_time_err', 'r2', 'n_images']].copy()
    
    focused_analysis['k_relative_error_percent'] = (focused_analysis['k_err'] / focused_analysis['k'] * 100).round(2)
    focused_analysis['dt_relative_error_percent'] = (focused_analysis['doubling_time_err'] / focused_analysis['doubling_time'] * 100).round(2)
    
    focused_analysis.to_excel(focused_path, index=False)
    print(f"Saved focused growth analysis to {focused_path}")

def main():
    """Main execution function"""
    print("Starting Growth Analysis...")
    print(f"Looking for data files in: {data_folder}")
    
    if not os.path.exists(data_folder):
        print(f"ERROR: Data folder '{data_folder}' does not exist!")
        print("Please create the folder and place your Excel files there.")
        return
    
    all_fit_results = []
    
    for cellline in celllines:
        print(f"\nProcessing cell line: {cellline}")
        
        all_data = process_data_files(cellline)
        
        if not all_data:
            print(f"  No data found for {cellline}")
            continue
        
        print("  Normalizing data...")
        normalized_data = normalize_and_truncate_data(all_data)
        
        print("  Creating fitted plots and extracting parameters...")
        fit_results = create_fitted_plots_and_save_params(cellline, normalized_data)
        all_fit_results.extend(fit_results)
    
    if all_fit_results:
        print("\nSaving fitting parameters to Excel...")
        save_fit_parameters_to_excel(all_fit_results)
        
        print(f"\nAnalysis complete! Results saved in '{output_dir}' directory.")
        print(f"Generated files:")
        print(f"  - Growth curve plots: {len([r for r in all_fit_results])} SVG files")
        print(f"  - Complete parameters: growth_analysis_parameters.xlsx")
        print(f"  - Focused analysis: growth_analysis_focused.xlsx")
        print(f"\nKey features:")
        print(f"  * X0 confluence mean and standard deviation")
        print(f"  * Raw maximum confluence at 92h endpoint")
        print(f"  * Growth fold change calculation (raw max/x0 ratio)")
        print(f"  * Coefficient of variation (CV) for confluences")
        print(f"  * Error propagation for growth parameters")
        print(f"  * Bootstrap analysis for parameter uncertainty")
    else:
        print("\nNo data was successfully processed. Please check:")
        print(f"  1. Data folder path: {data_folder}")
        print(f"  2. File naming convention: Growth_(celllinename).xlsx")
        print(f"  3. Excel file format: time in first column, data in subsequent columns")

if __name__ == "__main__":
    main()