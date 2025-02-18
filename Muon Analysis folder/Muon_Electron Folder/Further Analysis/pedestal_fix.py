import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.stats as stats
from matplotlib.backends.backend_pdf import PdfPages

# File paths
pedestals_file = 'calibrations/pedestals.csv'  # Old pedestal file
source_file = 'analysis_files/run_20220425_fpga_run.csv'  # Source data file
output_pdf = 'plots/pedestal_plots.pdf'  # Output PDF file

# Read pedestal data
pedestals_data = pd.read_csv(pedestals_file)
# Read source data
source_data = pd.read_csv(source_file)

# Create a PDF file to save all plots
with PdfPages(output_pdf) as pdf:
    # Group data by layer and strip
    for (layer, strip), group in source_data.groupby(['layer', 'strip']):
        # Retrieve pedestal and std_dev for end 0
        pedestal_info_end0 = pedestals_data[
            (pedestals_data['layer'] == layer) &
            (pedestals_data['strip'] == strip) &
            (pedestals_data['end'] == 0)
        ]
        
        # Retrieve pedestal and std_dev for end 1
        pedestal_info_end1 = pedestals_data[
            (pedestals_data['layer'] == layer) &
            (pedestals_data['strip'] == strip) &
            (pedestals_data['end'] == 1)
        ]
        
        if pedestal_info_end0.empty or pedestal_info_end1.empty:
            print(f"No pedestal information for Layer {layer}, Strip {strip}")
            continue
        
        # Get old pedestal and std_dev (rounded to 3 decimal places)
        pedestal_old_end0 = round(pedestal_info_end0['pedestal'].values[0], 3)
        std_dev_old_end0 = round(pedestal_info_end0['std_dev'].values[0], 3)
        
        pedestal_old_end1 = round(pedestal_info_end1['pedestal'].values[0], 3)
        std_dev_old_end1 = round(pedestal_info_end1['std_dev'].values[0], 3)

        # Define full display range: pedestal ± 400
        display_range_min = min(pedestal_old_end0, pedestal_old_end1) - 400
        display_range_max = max(pedestal_old_end0, pedestal_old_end1) + 400

        # Define fit range: pedestal ± 5 * std_dev
        fit_range_min_end0 = pedestal_old_end0 - 5 * std_dev_old_end0
        fit_range_max_end0 = pedestal_old_end0 + 5 * std_dev_old_end0
        
        fit_range_min_end1 = pedestal_old_end1 - 5 * std_dev_old_end1
        fit_range_max_end1 = pedestal_old_end1 + 5 * std_dev_old_end1

        # Select ADC sum data
        adc_sum_end0 = group['adc_sum_end0']
        adc_sum_end1 = group['adc_sum_end1']

        # Select data within the display range
        display_data_end0 = adc_sum_end0[(adc_sum_end0 >= display_range_min) & (adc_sum_end0 <= display_range_max)]
        display_data_end1 = adc_sum_end1[(adc_sum_end1 >= display_range_min) & (adc_sum_end1 <= display_range_max)]

        # Select data within the fit range
        fit_data_end0 = adc_sum_end0[(adc_sum_end0 >= fit_range_min_end0) & (adc_sum_end0 <= fit_range_max_end0)]
        fit_data_end1 = adc_sum_end1[(adc_sum_end1 >= fit_range_min_end1) & (adc_sum_end1 <= fit_range_max_end1)]
        
        # Gaussian fit for end 0
        mean_end0, std_dev_revised_end0 = stats.norm.fit(fit_data_end0)
        mean_end0 = round(mean_end0, 3)
        std_dev_revised_end0 = round(std_dev_revised_end0, 3)
        
        # Gaussian fit for end 1
        mean_end1, std_dev_revised_end1 = stats.norm.fit(fit_data_end1)
        mean_end1 = round(mean_end1, 3)
        std_dev_revised_end1 = round(std_dev_revised_end1, 3)

        # Plot end 0 data
        fig1, ax1 = plt.subplots(figsize=(8, 6))
        bins = np.linspace(display_range_min, display_range_max, 100)  # Use 100 bins for full range
        ax1.hist(display_data_end0, bins=bins, density=True, histtype='step', color='gray', label='Data')

        # Fit curve for end 0
        x0 = np.linspace(fit_range_min_end0, fit_range_max_end0, 100)
        p0 = stats.norm.pdf(x0, mean_end0, std_dev_revised_end0)
        ax1.plot(x0, p0, 'k', linewidth=2)
        ax1.set_title(f'Layer: {layer}, Strip: {strip}, End: 0')
        ax1.set_xlabel('ADC Sum End 0')
        ax1.set_ylabel('Density')
        ax1.legend()

        # Add pedestal information
        fig1.suptitle(f'Pedestal Old: {pedestal_old_end0}, Pedestal Revised: {mean_end0}\nStd Dev Old: {std_dev_old_end0}, Std Dev Revised: {std_dev_revised_end0}')
        
        # Save end 0 plot to PDF
        pdf.savefig(fig1)
        plt.close(fig1)

        # Plot end 1 data
        fig2, ax2 = plt.subplots(figsize=(8, 6))
        ax2.hist(display_data_end1, bins=bins, density=True, histtype='step', color='gray', label='Data')

        # Fit curve for end 1
        x1 = np.linspace(fit_range_min_end1, fit_range_max_end1, 100)
        p1 = stats.norm.pdf(x1, mean_end1, std_dev_revised_end1)
        ax2.plot(x1, p1, 'k', linewidth=2)
        ax2.set_title(f'Layer: {layer}, Strip: {strip}, End: 1')
        ax2.set_xlabel('ADC Sum End 1')
        ax2.set_ylabel('Density')
        ax2.legend()

        # Add pedestal information
        fig2.suptitle(f'Pedestal Old: {pedestal_old_end1}, Pedestal Revised: {mean_end1}\nStd Dev Old: {std_dev_old_end1}, Std Dev Revised: {std_dev_revised_end1}')

        # Save end 1 plot to PDF
        pdf.savefig(fig2)
        plt.close(fig2)

print(f"PDF saved to {output_pdf}")
