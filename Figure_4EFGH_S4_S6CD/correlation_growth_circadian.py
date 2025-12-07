"""

@author: Carolin Ector
Date: August 2025

"""

"""
Correlation analysis between growth and circadian parameters
Clock-Neuroblastoma Figure S4
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import pearsonr, linregress
import matplotlib.colors as mcolors
import itertools
import os
from scipy.stats import t
from matplotlib.colors import LinearSegmentedColormap

# Configuration
excel_file = 'neuroblastoma_circadian_data.xlsx'  # Update this to your Excel file name
sheet1_name = 'circ_values'  # Circadian parameters (update as needed)
sheet2_name = 'growth_values'  # Growth parameters (update as needed)

# Define enhanced color palette for cell lines (using colorbrewer-inspired colors)
cellline_colors = {
    'CHP212': '#1f77b4',   # Blue
    'SKNAS': '#ff7f0e',    # Orange  
    'SY5Y': '#2ca02c',     # Green
    'NGP': '#9467bd',      # Purple
    'GIMEN': '#8c564b',    # Brown
    'SKNSH': '#e377c2',    # Pink
    'CLBGA': '#7f7f7f',    # Gray
    'IMR5': '#bcbd22',     # Olive
    'SKNBE': '#17becf',    # Cyan
    'LAN5': '#aec7e8',     # Light Blue
    'SKNBE2': '#ffbb78',   # Light Orange
}

# Set global plotting parameters for better aesthetics
plt.rcParams.update({
    'font.size': 11,
    'font.family': 'Arial',
    'axes.linewidth': 1.2,
    'axes.spines.top': False,
    'axes.spines.right': False,
    'axes.grid': True,
    'grid.alpha': 0.3,
    'grid.linewidth': 0.8,
    'legend.frameon': True,
    'legend.fancybox': True,
    'legend.shadow': True,
    'legend.framealpha': 0.9,
    'legend.edgecolor': 'gray',
    'figure.dpi': 150,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'savefig.facecolor': 'white',
    'savefig.edgecolor': 'none'
})

# Exclude specific cell lines if needed
exclude_cells = []  # Add cell lines to exclude, e.g., ['SKNBE2']

# Create output directories
os.makedirs('growth_analysis_results/circadian_vs_growth', exist_ok=True)
os.makedirs('growth_analysis_results/circadian_vs_circadian', exist_ok=True)
os.makedirs('growth_analysis_results/heatmaps', exist_ok=True)

def load_data():
    """Load data from both sheets"""
    try:
        # Load circadian data (Sheet1)
        circadian_data = pd.read_excel(excel_file, sheet_name=sheet1_name, index_col=0)
        print(f"Loaded circadian data from {sheet1_name}: {circadian_data.shape}")
        print(f"Circadian parameters: {list(circadian_data.columns)}")
        
        # Load growth data (Sheet2)
        growth_data = pd.read_excel(excel_file, sheet_name=sheet2_name, index_col=0)
        print(f"Loaded growth data from {sheet2_name}: {growth_data.shape}")
        print(f"Growth parameters: {list(growth_data.columns)}")
        
        # Get numeric columns only
        circ_cols = [col for col in circadian_data.columns if pd.api.types.is_numeric_dtype(circadian_data[col])]
        growth_cols = [col for col in growth_data.columns if pd.api.types.is_numeric_dtype(growth_data[col])]
        
        circadian_data = circadian_data[circ_cols]
        growth_data = growth_data[growth_cols]
        
        print(f"Numeric circadian parameters: {circ_cols}")
        print(f"Numeric growth parameters: {growth_cols}")
        
        return circadian_data, growth_data, circ_cols, growth_cols
        
    except Exception as e:
        print(f"Error loading data: {e}")
        return None, None, [], []

def create_scatter_plot(x_data, y_data, x_label, y_label, cell_names, title, filename):
    """Create enhanced scatter plot with linear regression and 95% CIs"""
    mask = ~np.isnan(x_data) & ~np.isnan(y_data)
    x_clean = x_data[mask]
    y_clean = y_data[mask]
    cells_clean = np.array(cell_names)[mask]
    
    if len(x_clean) < 3:
        print(f"Skipping {title}: insufficient data points ({len(x_clean)})")
        return None, None, None
    
    fig, ax = plt.subplots(figsize=(10, 7))
    
    for xi, yi, cell in zip(x_clean, y_clean, cells_clean):
        color = cellline_colors.get(cell, '#2c2c2c')
        ax.scatter(xi, yi, color=color, s=120, label=cell, alpha=0.85, 
                   edgecolor='white', linewidth=2, zorder=3)

    # NEW: compute CI info
    ci = regression_ci(x_clean, y_clean, confidence=0.95)
    slope = ci['slope']
    intercept = ci['intercept']
    r = ci['r']
    p = ci['p']
    r2 = r**2

    # regression line
    ax.plot(ci['x_grid'], ci['y_grid'], color='#2c2c2c', linewidth=3, alpha=0.9,
            label='Linear fit', zorder=2)

    # 95% CI of the mean fit
    ax.fill_between(ci['x_grid'], ci['line_ci_lower'], ci['line_ci_upper'],
                    alpha=0.20, zorder=1, label='95% CI')

    # stats box, includes slope CI
    if p < 0.001:
        p_text = "p < 0.001"
    elif p < 0.01:
        p_text = f"p = {p:.3f}"
    else:
        p_text = f"p = {p:.3f}"

    stats_text = (
        f"R² = {r2:.3f}\n"
        f"Slope = {slope:.3f}\n"
        f"95% CI: [{ci['slope_ci_lower']:.3f}, {ci['slope_ci_upper']:.3f}]\n"
        f"{p_text}\n"
        f"n = {ci['n']}"
    )

    if r > 0:
        stats_x, stats_y, ha = 0.05, 0.95, 'left'
    else:
        stats_x, stats_y, ha = 0.95, 0.95, 'right'

    ax.text(stats_x, stats_y, stats_text, transform=ax.transAxes,
            va='top', ha=ha, fontsize=12, fontweight='bold',
            bbox=dict(boxstyle='round,pad=0.6', facecolor='white', alpha=0.9,
                      edgecolor='#cccccc', linewidth=1.5))

    ax.set_xlabel(x_label, fontsize=14, fontweight='bold', labelpad=10)
    ax.set_ylabel(y_label, fontsize=14, fontweight='bold', labelpad=10)
    ax.set_title(title, fontsize=16, fontweight='bold', pad=20)

    handles, labels = ax.get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    if 'Linear fit' in by_label and '95% CI' in by_label:
        # keep only cell labels plus CI label
        pass
    legend = ax.legend(by_label.values(), by_label.keys(),
                       bbox_to_anchor=(1.02, 1), loc='upper left',
                       fontsize=11, title='Cell Lines', title_fontsize=12,
                       markerscale=1.2, handletextpad=0.5, columnspacing=1.0)
    legend.get_title().set_fontweight('bold')

    x_range = np.max(x_clean) - np.min(x_clean)
    y_range = np.max(y_clean) - np.min(y_clean)
    ax.set_xlim(np.min(x_clean) - 0.05 * x_range, np.max(x_clean) + 0.05 * x_range)
    ax.set_ylim(np.min(y_clean) - 0.05 * y_range, np.max(y_clean) + 0.05 * y_range)

    ax.tick_params(axis='both', which='major', labelsize=11, width=1.2, length=6)
    ax.grid(True, alpha=0.3, linewidth=0.8)

    plt.tight_layout()
    plt.savefig(filename, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()

    return r, p, len(x_clean)


def calculate_correlations(data1, data2, cols1, cols2):
    """Calculate Pearson correlations between all parameter combinations"""
    n_params1 = len(cols1)
    n_params2 = len(cols2)
    
    # Initialize matrices
    pearson_r = np.full((n_params1, n_params2), np.nan)
    pearson_p = np.full((n_params1, n_params2), np.nan)
    sample_sizes = np.zeros((n_params1, n_params2))
    
    for i, param1 in enumerate(cols1):
        for j, param2 in enumerate(cols2):
            # Find common cell lines
            common_cells = data1.index.intersection(data2.index)
            common_cells = [cell for cell in common_cells if cell not in exclude_cells]
            
            if len(common_cells) < 3:
                continue
            
            # Get data for common cell lines
            x = data1.loc[common_cells, param1].values
            y = data2.loc[common_cells, param2].values
            
            # Remove NaN values
            mask = ~np.isnan(x) & ~np.isnan(y)
            if np.sum(mask) < 3:
                continue
            
            x_clean = x[mask]
            y_clean = y[mask]
            sample_sizes[i, j] = len(x_clean)
            
            # Calculate correlations
            try:
                pearson_r[i, j], pearson_p[i, j] = pearsonr(x_clean, y_clean)
            except:
                continue
    
    return pearson_r, pearson_p, sample_sizes

def pval_to_stars(p):
    """Convert p-values to significance stars"""
    if np.isnan(p):
        return ''
    elif p <= 0.0001:
        return '****'
    elif p <= 0.001:
        return '***'
    elif p <= 0.01:
        return '**'
    elif p <= 0.05:
        return '*'
    else:
        return ''

def regression_ci(x, y, confidence=0.95, grid_pts=200):
    n = len(x)
    if n < 3:
        return None

    # fit
    slope, intercept, r, p, stderr = linregress(x, y)

    # residual variance
    yhat = slope * x + intercept
    resid = y - yhat
    mse = np.sum(resid**2) / (n - 2)

    # design stats
    x_mean = np.mean(x)
    ss_x = np.sum((x - x_mean) ** 2)

    # t critical
    alpha = 1 - confidence
    t_val = t.ppf(1 - alpha/2, n - 2)

    # slope CI
    se_slope = np.sqrt(mse / ss_x)
    slope_ci_lower = slope - t_val * se_slope
    slope_ci_upper = slope + t_val * se_slope

    # grid across x range
    x_grid = np.linspace(np.min(x), np.max(x), grid_pts)
    y_grid = slope * x_grid + intercept

    # CI of the mean fit
    se_mean = np.sqrt(mse * (1/n + (x_grid - x_mean)**2 / ss_x))
    line_ci_lower = y_grid - t_val * se_mean
    line_ci_upper = y_grid + t_val * se_mean

    # prediction interval
    se_pred = np.sqrt(mse * (1 + 1/n + (x_grid - x_mean)**2 / ss_x))
    pred_ci_lower = y_grid - t_val * se_pred
    pred_ci_upper = y_grid + t_val * se_pred

    return dict(
        slope=slope,
        intercept=intercept,
        r=r,
        p=p,
        n=n,
        se_slope=se_slope,
        slope_ci_lower=slope_ci_lower,
        slope_ci_upper=slope_ci_upper,
        x_grid=x_grid,
        y_grid=y_grid,
        line_ci_lower=line_ci_lower,
        line_ci_upper=line_ci_upper,
        pred_ci_lower=pred_ci_lower,
        pred_ci_upper=pred_ci_upper
    )

def create_heatmap(corr_matrix, p_matrix, row_labels, col_labels, title, filename):
    """Create enhanced correlation heatmap with significance annotations and dynamic colors"""
    # Create annotation matrix with stars
    annot_matrix = np.empty(corr_matrix.shape, dtype=object)
    for i in range(corr_matrix.shape[0]):
        for j in range(corr_matrix.shape[1]):
            val = corr_matrix[i, j]
            stars = pval_to_stars(p_matrix[i, j])
            if np.isnan(val):
                annot_matrix[i, j] = ''
            else:
                annot_matrix[i, j] = f'{val:.2f}{stars if stars else ""}'
    
    # Dynamic limits based on data
    vmin = np.nanmin(corr_matrix)
    vmax = np.nanmax(corr_matrix)

    # Define green–white–magenta colormap
    cmap = LinearSegmentedColormap.from_list(
        "green_magenta",
        ["#a20056", "white", "#146730"],  # magenta → white → green
        N=256
    )
    
    # Plot
    fig, ax = plt.subplots(figsize=(max(10, len(col_labels) * 1.2),
                                    max(8, len(row_labels) * 0.6)))
    
    heatmap = sns.heatmap(corr_matrix, annot=annot_matrix, fmt='', 
                         xticklabels=col_labels, yticklabels=row_labels,
                         cmap=cmap, center=0, vmin=vmin, vmax=vmax,
                         annot_kws={'fontsize': 10, 'ha': 'center', 'fontweight': 'bold'},
                         cbar_kws={'label': 'Pearson r', 'shrink': 0.8, 'aspect': 30},
                         linewidths=0.5, linecolor='white',
                         square=True)
    
    # Style
    cbar = heatmap.collections[0].colorbar
    cbar.ax.tick_params(labelsize=11)
    cbar.set_label('Pearson r', fontsize=12, fontweight='bold', labelpad=15)

    ax.set_title(title, fontsize=16, fontweight='bold', pad=25)
    ax.set_xlabel('Parameters', fontsize=14, fontweight='bold', labelpad=15)
    ax.set_ylabel('Parameters', fontsize=14, fontweight='bold', labelpad=15)
    
    ax.set_xticklabels(col_labels, rotation=45, ha='right', fontsize=11)
    ax.set_yticklabels(row_labels, rotation=0, ha='right', fontsize=11)
    
    plt.tight_layout()
    plt.savefig(filename, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()

def main():
    print("Starting comprehensive correlation analysis...")
    
    # Load data
    circadian_data, growth_data, circ_cols, growth_cols = load_data()
    
    if circadian_data is None or growth_data is None:
        print("Failed to load data. Exiting.")
        return
    
    print(f"\\nAnalyzing {len(circ_cols)} circadian parameters and {len(growth_cols)} growth parameters")
    print(f"Cell lines in circadian data: {list(circadian_data.index)}")
    print(f"Cell lines in growth data: {list(growth_data.index)}")
    
    # Find common cell lines
    common_cells = circadian_data.index.intersection(growth_data.index)
    common_cells = [cell for cell in common_cells if cell not in exclude_cells]
    print(f"Common cell lines for analysis: {common_cells}")
    
    # 1. CIRCADIAN vs GROWTH CORRELATIONS
    print("\\n=== 1. Analyzing Circadian vs Growth correlations ===")
    
    # Create scatter plots for all circadian vs growth combinations
    circ_growth_results = []
    print(f"Creating {len(circ_cols) * len(growth_cols)} scatter plots for circadian vs growth...")
    
    for i, (circ_param, growth_param) in enumerate(itertools.product(circ_cols, growth_cols)):
        if (i + 1) % 10 == 0:  # Progress indicator
            print(f"  Progress: {i + 1}/{len(circ_cols) * len(growth_cols)} plots")
            
        x_data = circadian_data.loc[common_cells, circ_param].values
        y_data = growth_data.loc[common_cells, growth_param].values
        
        title = f'{circ_param} vs {growth_param}'
        # Clean filename for better organization
        safe_circ = circ_param.replace('/', '_').replace(' ', '_')
        safe_growth = growth_param.replace('/', '_').replace(' ', '_')
        filename = f'growth_analysis_results/circadian_vs_growth/{safe_circ}_vs_{safe_growth}.svg'
        
        r, p, n = create_scatter_plot(x_data, y_data, circ_param, growth_param, 
                                    common_cells, title, filename)
        
        if r is not None:
            circ_growth_results.append({
                'circadian_param': circ_param,
                'growth_param': growth_param,
                'pearson_r': r,
                'p_value': p,
                'n_samples': n
            })
    
    # Calculate correlation matrices for circadian vs growth
    pearson_r_cg, pearson_p_cg, sizes_cg = calculate_correlations(
        circadian_data, growth_data, circ_cols, growth_cols)
    
    # Create heatmaps for circadian vs growth
    create_heatmap(pearson_r_cg, pearson_p_cg, circ_cols, growth_cols,
                   'Circadian vs Growth - Pearson Correlations',
                   'growth_analysis_results/heatmaps/circadian_vs_growth_pearson.svg')
    
    # 2. CIRCADIAN vs CIRCADIAN CORRELATIONS
    print("\\n=== 2. Analyzing Circadian vs Circadian correlations ===")
    
    # Create scatter plots for all circadian vs circadian combinations
    circ_circ_results = []
    n_combinations = len(circ_cols) * (len(circ_cols) - 1) // 2  # Only upper triangle
    print(f"Creating {n_combinations} scatter plots for circadian vs circadian...")
    
    plot_count = 0
    for i, circ_param1 in enumerate(circ_cols):
        for j, circ_param2 in enumerate(circ_cols):
            if i >= j:  # Skip diagonal and duplicate combinations
                continue
            
            plot_count += 1
            if plot_count % 5 == 0:  # Progress indicator
                print(f"  Progress: {plot_count}/{n_combinations} plots")
                
            x_data = circadian_data.loc[common_cells, circ_param1].values
            y_data = circadian_data.loc[common_cells, circ_param2].values
            
            title = f'{circ_param1} vs {circ_param2}'
            # Clean filename for better organization
            safe_param1 = circ_param1.replace('/', '_').replace(' ', '_')
            safe_param2 = circ_param2.replace('/', '_').replace(' ', '_')
            filename = f'growth_analysis_results/circadian_vs_circadian/{safe_param1}_vs_{safe_param2}.svg'
            
            r, p, n = create_scatter_plot(x_data, y_data, circ_param1, circ_param2,
                                        common_cells, title, filename)
            
            if r is not None:
                circ_circ_results.append({
                    'param1': circ_param1,
                    'param2': circ_param2,
                    'pearson_r': r,
                    'p_value': p,
                    'n_samples': n
                })
    
    # Calculate correlation matrices for circadian vs circadian
    pearson_r_cc, pearson_p_cc, sizes_cc = calculate_correlations(
        circadian_data, circadian_data, circ_cols, circ_cols)
    
    # Create heatmaps for circadian vs circadian
    create_heatmap(pearson_r_cc, pearson_p_cc, circ_cols, circ_cols,
                   'Circadian vs Circadian - Pearson Correlations',
                   'growth_analysis_results/heatmaps/circadian_vs_circadian_pearson.svg')

    # 3. SAVE RESULTS TO EXCEL
    print("\\n=== 3. Saving results to Excel ===")
    
    with pd.ExcelWriter('correlation_analysis_results.xlsx') as writer:
        # Circadian vs Growth results
        if circ_growth_results:
            df_cg = pd.DataFrame(circ_growth_results)
            df_cg.to_excel(writer, sheet_name='Circadian_vs_Growth_Stats', index=False)
        
        # Circadian vs Circadian results
        if circ_circ_results:
            df_cc = pd.DataFrame(circ_circ_results)
            df_cc.to_excel(writer, sheet_name='Circadian_vs_Circadian_Stats', index=False)
        
        # Correlation matrices
        pd.DataFrame(pearson_r_cg, index=circ_cols, columns=growth_cols).to_excel(
            writer, sheet_name='Circ_Growth_Pearson_R')
        pd.DataFrame(pearson_p_cg, index=circ_cols, columns=growth_cols).to_excel(
            writer, sheet_name='Circ_Growth_Pearson_P')
        
        pd.DataFrame(pearson_r_cc, index=circ_cols, columns=circ_cols).to_excel(
            writer, sheet_name='Circ_Circ_Pearson_R')
        pd.DataFrame(pearson_p_cc, index=circ_cols, columns=circ_cols).to_excel(
            writer, sheet_name='Circ_Circ_Pearson_P')
    
    # 4. SUMMARY STATISTICS
    print("\\n=== 4. Summary Statistics ===")
    print(f"Total circadian vs growth correlations analyzed: {len(circ_growth_results)}")
    print(f"Total circadian vs circadian correlations analyzed: {len(circ_circ_results)}")
    
    if circ_growth_results:
        cg_df = pd.DataFrame(circ_growth_results)
        sig_cg = cg_df[cg_df['p_value'] <= 0.05]
        print(f"Significant circadian vs growth correlations (p≤0.05): {len(sig_cg)}")
        if len(sig_cg) > 0:
            print("Top significant circadian vs growth correlations:")
            print(sig_cg.nsmallest(5, 'p_value')[['circadian_param', 'growth_param', 'pearson_r', 'p_value']])
    
    if circ_circ_results:
        cc_df = pd.DataFrame(circ_circ_results)
        sig_cc = cc_df[cc_df['p_value'] <= 0.05]
        print(f"\\nSignificant circadian vs circadian correlations (p≤0.05): {len(sig_cc)}")
        if len(sig_cc) > 0:
            print("Top significant circadian vs circadian correlations:")
            print(sig_cc.nsmallest(5, 'p_value')[['param1', 'param2', 'pearson_r', 'p_value']])
    
    print("\\n=== Analysis Complete! ===")
    print("Check the following directories for results:")
    print("- growth_analysis_results/circadian_vs_growth/ - Individual scatter plots")
    print("- growth_analysis_results/circadian_vs_circadian/ - Individual scatter plots")
    print("- growth_analysis_results/heatmaps/ - Correlation heatmaps")
    print("- correlation_analysis_results.xlsx - All numerical results")

if __name__ == "__main__":
    main()