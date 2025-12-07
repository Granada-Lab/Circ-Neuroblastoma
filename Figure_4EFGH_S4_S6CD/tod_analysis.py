#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Time-of-Day (ToD) Response Analysis 
Author: Carolin Ector
Date: August 2025
Clock-Neuroblastoma Figure 4 & Figure S6

This script analyzes cell confluence response to drug treatments at different times of day.
It performs:
1. Data loading and preprocessing (outlier filtering, normalization)
2. Confluence curve visualization
3. Time-of-Day response calculation (max-min from interpolated curve)
4. Processing of two experiments

Key features:
- Processes two experiments (experiment1 and experiment2)
- Focuses solely on confluence measurements
- Calculates ToD response amplitudes without curve fitting
- Exports results per timepoint per experiment
"""

import os
import re
import logging
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Any
from dataclasses import dataclass
from collections import defaultdict

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.signal import savgol_filter
from scipy.interpolate import PchipInterpolator

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

@dataclass
class ExperimentConfig:
    """Configuration parameters for the ToD analysis experiment."""
    base_data_dir: str = "Data/time_of_day_raw_data"
    cell_lines: List[str] = None
    drugs: List[str] = None
    timepoints: List[str] = None
    treatment_times: Dict[str, int] = None
    n_measurements: int = 46  # Number of measurements post-treatment
    measurement_interval_hours: int = 2  # Time between measurements
    outlier_threshold: float = 3.0  # MAD threshold for outlier detection
    experiments: List[str] = None
    
    def __post_init__(self):
        if self.cell_lines is None:
            self.cell_lines = ["CHP212", "SY5Y", "GIMEN", "SKNSH", "SKNBE2"]
        
        if self.drugs is None:
            self.drugs = ["Lorlatinib", "Doxorubicin", "Cisplatin"]
        
        if self.timepoints is None:
            self.timepoints = ["0h", "4h", "8h", "16h", "20h", "24h"]
        
        if self.treatment_times is None:
            self.treatment_times = {
                "0h": 46, "4h": 46, "8h": 46, "UNRESET": 46,
                "16h": 60, "20h": 60, "24h": 60
            }
        
        if self.experiments is None:
            self.experiments = ["experiment1", "experiment2"]
        
        self.conditions = self.timepoints + ["UNRESET"]

class DataProcessor:
    """Handles data loading, cleaning, and preprocessing."""
    
    def __init__(self, config: ExperimentConfig):
        self.config = config
    
    def load_excel_data(self, filepath: str) -> Optional[pd.DataFrame]:
        """Load and preprocess Excel data file."""
        try:
            if not os.path.exists(filepath):
                logger.warning(f"File not found: {filepath}")
                return None
            
            df = pd.read_excel(filepath)
            
            # Standardize column names
            df.rename(columns={df.columns[0]: "Time"}, inplace=True)
            
            # Clean empty cells and convert to numeric
            df.replace(r'^\s*$', np.nan, regex=True, inplace=True)
            df[df.columns[1:]] = df[df.columns[1:]].apply(pd.to_numeric, errors="coerce")
            
            # Interpolate missing values
            df.interpolate(method="linear", axis=0, inplace=True)
            
            logger.info(f"Successfully loaded {filepath}: {df.shape}")
            return df
            
        except Exception as e:
            logger.error(f"Error loading {filepath}: {e}")
            return None
    
    def filter_outliers(self, data: pd.DataFrame, threshold: float = None) -> pd.DataFrame:
        """Filter outlier wells based on baseline confluence using MAD."""
        if threshold is None:
            threshold = self.config.outlier_threshold
        
        baseline = data.iloc[0]  # First measurement as baseline
        median_baseline = np.median(baseline)
        mad = np.median(np.abs(baseline - median_baseline))
        
        # Keep wells above threshold
        outlier_threshold = median_baseline - threshold * mad
        valid_wells = baseline[baseline > outlier_threshold].index
        
        logger.debug(f"Filtered {len(baseline) - len(valid_wells)} outlier wells")
        return data[valid_wells]
    
    def normalize_to_baseline(self, data: pd.DataFrame) -> pd.DataFrame:
        """Normalize data to the first timepoint (treatment time)."""
        return data.div(data.iloc[0])

class TodResponseCalculator:
    """Calculates Time-of-Day response metrics."""
    
    def __init__(self, config: ExperimentConfig):
        self.config = config
    
    def calculate_tod_response(self, timepoints: List[str], values: List[float], 
                             control_value: float) -> Tuple[float, Dict[str, float]]:
        """Calculate ToD response amplitude (max-min) and individual timepoint responses."""
        if len(timepoints) < 3 or np.isnan(control_value) or control_value <= 0:
            return np.nan, {}
        
        # Convert timepoints to hours and normalize values
        hours = np.array([int(tp.replace("h", "")) for tp in timepoints])
        normalized_values = np.array(values) / control_value
        
        # Remove NaN values
        valid_mask = ~np.isnan(normalized_values)
        if np.sum(valid_mask) < 3:
            return np.nan, {}
        
        hours_valid = hours[valid_mask]
        values_valid = normalized_values[valid_mask]
        timepoints_valid = [tp for i, tp in enumerate(timepoints) if valid_mask[i]]
        
        # Create individual timepoint responses dictionary
        timepoint_responses = {tp: val for tp, val in zip(timepoints_valid, values_valid)}
        
        # Interpolate and smooth for amplitude calculation
        interpolator = PchipInterpolator(hours_valid, values_valid)
        hours_dense = np.linspace(hours_valid.min(), hours_valid.max(), 200)
        values_dense = interpolator(hours_dense)
        values_smooth = savgol_filter(values_dense, window_length=21, polyorder=2, mode="interp")
        
        # Calculate ToD-MR amplitude (max - min)
        tod_amplitude = values_smooth.max() - values_smooth.min()
        
        return tod_amplitude, timepoint_responses

class PlotGenerator:
    """Generates plots for ToD analysis."""
    
    def __init__(self, config: ExperimentConfig, output_dir: str = "plots"):
        self.config = config
        self.output_dir = Path(output_dir)
        self._create_directories()
    
    def _create_directories(self):
        """Create necessary output directories."""
        subdirs = ["raw_curves", "normalized_curves"]
        for experiment in self.config.experiments:
            for subdir in subdirs:
                (self.output_dir / experiment / subdir).mkdir(parents=True, exist_ok=True)
    
    def plot_raw_curves(self, cell: str, drug: str, condition_data: Dict[str, Tuple], experiment: str):
        """Plot raw confluence curves for all conditions."""
        plt.figure(figsize=(12, 8))
        
        for condition, (time_vals, mean_vals, std_vals, n_wells, treatment_conf) in condition_data.items():
            label = f"{condition} ({treatment_conf:.1f}%) [n={n_wells}]"
            plt.plot(time_vals, mean_vals, label=label, linewidth=2)
            plt.fill_between(time_vals, mean_vals - std_vals, mean_vals + std_vals, alpha=0.3)
        
        plt.title(f"{experiment} - {cell} - {drug}: Raw Confluence Curves")
        plt.xlabel("Time (hours)")
        plt.ylabel("Confluence (%)")
        plt.legend()
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.savefig(self.output_dir / experiment / "raw_curves" / f"{cell}_{drug}_raw.svg")
        plt.close()
    
    def plot_normalized_curves(self, cell: str, drug: str, condition_data: Dict[str, Tuple], experiment: str):
        """Plot normalized confluence curves."""
        plt.figure(figsize=(12, 8))
        
        for condition, (time_vals, mean_vals, std_vals, n_wells, treatment_conf) in condition_data.items():
            label = f"{condition} ({treatment_conf:.1f}%) [n={n_wells}]"
            plt.plot(time_vals, mean_vals, label=label, linewidth=2)
            plt.fill_between(time_vals, mean_vals - std_vals, mean_vals + std_vals, alpha=0.3)
        
        plt.title(f"{experiment} - {cell} - {drug}: Normalized Confluence Curves")
        plt.xlabel("Time Post-Treatment (hours)")
        plt.ylabel("Normalized Confluence")
        plt.legend()
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.savefig(self.output_dir / experiment / "normalized_curves" / f"{cell}_{drug}_normalized.svg")
        plt.close()

class TodAnalyzer:
    """Main class orchestrating the ToD analysis workflow."""
    
    def __init__(self, config: ExperimentConfig = None):
        self.config = config or ExperimentConfig()
        self.data_processor = DataProcessor(self.config)
        self.tod_calculator = TodResponseCalculator(self.config)
        self.plot_generator = PlotGenerator(self.config)
        
        # Results storage for both experiments
        self.confluence_results = defaultdict(lambda: defaultdict(lambda: defaultdict(float)))  # [experiment][drug][cell]
        self.timepoint_responses = defaultdict(lambda: defaultdict(lambda: defaultdict(dict)))  # [experiment][drug][cell][timepoint]
    
    def analyze_single_experiment(self, cell: str, drug: str, experiment: str) -> bool:
        """Analyze a single cell line + drug combination for one experiment."""
        base = Path(self.config.base_data_dir)
        filepath = base / experiment / f"{cell}_{drug}.xlsx"
        logger.info(f"Analyzing {experiment} - {cell} - {drug} from {filepath}")
    
        df = self.data_processor.load_excel_data(str(filepath))
        if df is None:
            return False
    
        raw_data = {}
        normalized_data = {}
        confluence_values = {}
    
        for condition in self.config.conditions:
            try:
                result = self._process_condition(df, condition, cell, drug)
                if result is None:
                    continue
                raw_stats, norm_stats, confluence_value = result
                raw_data[condition] = raw_stats
                normalized_data[condition] = norm_stats
                confluence_values[condition] = confluence_value
            except Exception as e:
                logger.error(f"Error processing {condition} for {experiment}-{cell}-{drug}: {e}")
    
        if not confluence_values:
            logger.warning(f"No valid data for {experiment}-{cell}-{drug}")
            return False
    
        if raw_data:
            self.plot_generator.plot_raw_curves(cell, drug, raw_data, experiment)
        if normalized_data:
            self.plot_generator.plot_normalized_curves(cell, drug, normalized_data, experiment)
    
        self._calculate_tod_responses(cell, drug, experiment, confluence_values)
        return True

    def _process_condition(self, df: pd.DataFrame, condition: str, cell: str, drug: str):
        """Process a single experimental condition."""
        treatment_time = self.config.treatment_times.get(condition)
        if treatment_time is None:
            return None
        
        # Find columns for this condition
        pattern = rf"\b{re.escape(condition)}\b"
        columns = df.columns[df.columns.str.contains(pattern)]
        if len(columns) == 0:
            return None
        
        # Extract data from treatment time onwards
        data_subset = df.loc[df["Time"] >= treatment_time, ["Time"] + list(columns)].dropna()
        if data_subset.empty or len(data_subset) < self.config.n_measurements:
            return None
        
        # Take exactly n_measurements
        time_values = data_subset["Time"].values[:self.config.n_measurements]
        well_data = data_subset[columns].iloc[:self.config.n_measurements]
        
        # Create post-treatment time array
        time_post_treatment = np.arange(0, self.config.n_measurements * self.config.measurement_interval_hours, 
                                       self.config.measurement_interval_hours)
        
        # Filter outliers
        filtered_data = self.data_processor.filter_outliers(well_data)
        if filtered_data.empty:
            return None
        
        # Calculate statistics for raw data
        mean_raw = filtered_data.mean(axis=1)
        std_raw = filtered_data.std(axis=1)
        treatment_confluence = mean_raw.iloc[0]
        n_wells = len(filtered_data.columns)
        
        # Normalize data
        normalized_data = self.data_processor.normalize_to_baseline(filtered_data)
        mean_norm = normalized_data.mean(axis=1)
        std_norm = normalized_data.std(axis=1)
        
        # Calculate final confluence (average of last 3 measurements)
        final_confluence = mean_norm.iloc[-3:].mean()
        
        # Prepare return data
        raw_stats = (time_values, mean_raw, std_raw, n_wells, treatment_confluence)
        norm_stats = (time_post_treatment, mean_norm, std_norm, n_wells, treatment_confluence)
        
        return raw_stats, norm_stats, final_confluence
    
    def _calculate_tod_responses(self, cell: str, drug: str, experiment: str, confluence_values: Dict):
        """Calculate ToD responses for confluence."""
        # Confluence ToD response
        control_confluence = confluence_values.get("0h", np.nan)
        if not np.isnan(control_confluence) and control_confluence > 0:
            timepoint_values = [confluence_values.get(tp, np.nan) for tp in self.config.timepoints]
            
            amplitude, timepoint_responses = self.tod_calculator.calculate_tod_response(
                self.config.timepoints, timepoint_values, control_confluence
            )
            
            if not np.isnan(amplitude):
                self.confluence_results[experiment][drug][cell] = amplitude
                self.timepoint_responses[experiment][drug][cell] = timepoint_responses
    
    def run_full_analysis(self):
        logger.info("Starting ToD response analysis for both experiments...")
        data_dir = Path(self.config.base_data_dir)

        if not data_dir.exists():
            logger.error(f"Data directory not found: {data_dir}")
            logger.error(f"Current working directory: {Path.cwd()}")
            logger.error(f"Here is what I see: {list(Path('.').iterdir())}")
            return
        
        # Check if the time_of_day_raw_data folder exists
        data_dir = os.path.join("Data", "time_of_day_raw_data")
        if not os.path.exists(data_dir):
            logger.error(f"Data directory not found: {data_dir}")
            logger.error("Available folders in Data: " + str(os.listdir("Data")))
            return
        
        total_successful = 0
        total_analyses = len(self.config.experiments) * len(self.config.cell_lines) * len(self.config.drugs)
        
        # Analyze each experiment
        for experiment in self.config.experiments:
            logger.info(f"Processing {experiment}...")
            experiment_dir = os.path.join(data_dir, experiment)
            if not os.path.exists(experiment_dir):
                logger.warning(f"Experiment directory not found: {experiment_dir}")
                continue
                
            successful_analyses = 0
            
            # Analyze each cell line + drug combination
            for cell in self.config.cell_lines:
                for drug in self.config.drugs:
                    if self.analyze_single_experiment(cell, drug, experiment):
                        successful_analyses += 1
                        total_successful += 1
            
            logger.info(f"Completed {successful_analyses}/{len(self.config.cell_lines) * len(self.config.drugs)} analyses for {experiment}")
        
        logger.info(f"Total completed: {total_successful}/{total_analyses} analyses")
        
        # Export results
        self._export_results()
        
        logger.info("Analysis complete!")
    
    def _export_results(self):
        """Export results to Excel files."""
        logger.info("Exporting results to Excel...")
        
        output_dir = Path("tod_results")
        output_dir.mkdir(exist_ok=True)
        
        # Export confluence amplitude results for each experiment
        for experiment in self.config.experiments:
            if experiment in self.confluence_results:
                confluence_df = pd.DataFrame(self.confluence_results[experiment]).reindex(
                    index=self.config.cell_lines, columns=self.config.drugs
                )
                confluence_df.to_excel(output_dir / f"{experiment}_tod_confluence_amplitudes.xlsx")
        
        # Export individual timepoint responses for each experiment
        for experiment in self.config.experiments:
            if experiment in self.timepoint_responses:
                # Create a comprehensive DataFrame with all timepoints
                all_responses = []
                for drug in self.config.drugs:
                    for cell in self.config.cell_lines:
                        if cell in self.timepoint_responses[experiment][drug]:
                            for timepoint, response in self.timepoint_responses[experiment][drug][cell].items():
                                all_responses.append({
                                    'Experiment': experiment,
                                    'Drug': drug,
                                    'Cell_Line': cell,
                                    'Timepoint': timepoint,
                                    'Relative_Confluence': response
                                })
                
                if all_responses:
                    responses_df = pd.DataFrame(all_responses)
                    responses_df.to_excel(output_dir / f"{experiment}_timepoint_responses.xlsx", index=False)
        
        # Export combined results across experiments
        if len(self.config.experiments) > 1:
            # Combined amplitude comparison
            combined_amplitudes = []
            for experiment in self.config.experiments:
                if experiment in self.confluence_results:
                    for drug in self.config.drugs:
                        for cell in self.config.cell_lines:
                            if cell in self.confluence_results[experiment][drug]:
                                combined_amplitudes.append({
                                    'Experiment': experiment,
                                    'Drug': drug,
                                    'Cell_Line': cell,
                                    'ToD_Amplitude': self.confluence_results[experiment][drug][cell]
                                })
            
            if combined_amplitudes:
                combined_df = pd.DataFrame(combined_amplitudes)
                combined_df.to_excel(output_dir / "combined_experiments_amplitudes.xlsx", index=False)

def main():
    """Main execution function."""
    # Create configuration
    config = ExperimentConfig()
    
    # Initialize analyzer
    analyzer = TodAnalyzer(config)
    
    # Run analysis
    analyzer.run_full_analysis()

if __name__ == "__main__":
    main()