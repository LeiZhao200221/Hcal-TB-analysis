import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import pylandau
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.ticker import MultipleLocator

class calculateMIPsManual:
    def __init__(self, data_file_name, pedestal_file_name, temp_mip_fit_cut_file, out_directory, plots_directory):
        self.data_file_name = data_file_name
        self.pedestal_file_name = pedestal_file_name
        self.temp_mip_fit_cut_file = temp_mip_fit_cut_file
        self.out_directory = out_directory
        self.plots_directory = plots_directory

        self.in_data = pd.read_csv(self.data_file_name)
        self.in_peds = pd.read_csv(self.pedestal_file_name)
        self.in_mip_fit_cut = pd.read_csv(self.temp_mip_fit_cut_file)

        if not os.path.exists(self.out_directory):
            os.makedirs(self.out_directory)
        if not os.path.exists(self.plots_directory):
            os.makedirs(self.plots_directory)

    def calculate_mip(self):
        pdf_path = os.path.join(self.plots_directory, 'all_plots.pdf')
        new_fit_cut = []
        mip_results = []

        with PdfPages(pdf_path) as pdf:
            for _, row in self.in_mip_fit_cut.iterrows():
                layer = row['layer']
                strip = row['strip']
                end = int(row['end'])
                multiplier = row['multiplier']
                fit_low = row['low']
                fit_upper = row['upper']

                # Get pedestal and standard deviation
                selection_ped = (self.in_peds['strip'] == strip) & (self.in_peds['end'] == end) & (self.in_peds['layer'] == layer)
                pedestal_ = self.in_peds[selection_ped]['pedestal'].iloc[0]
                std_dev = self.in_peds[selection_ped]['std_dev'].iloc[0]

                # Select data for the given layer and strip
                selection_data = (self.in_data['layer'] == layer) & (self.in_data['strip'] == strip)
                dataframe = self.in_data[selection_data].copy()
                adc_column = 'adc_sum_end' + str(end)

                if adc_column not in dataframe.columns:
                    print(f"Column {adc_column} not found in data for Layer: {layer}, Bar: {strip}, Side: {end}")
                    continue

                # Subtract pedestal from ADC sum
                dataframe[adc_column] = dataframe[adc_column].astype(float) - pedestal_

                # Filter based on multiplier * std_dev
                dataframe = dataframe[dataframe[adc_column] > multiplier * std_dev]

                # Create histogram
                counts, bins = np.histogram(dataframe[adc_column], bins=50, range=[0, 2000])
                y = counts
                x = (bins[:-1] + bins[1:]) / 2

                # Select the range for fitting
                new_y = y[(x >= fit_low) & (x <= fit_upper) & (y != 0)]
                new_x = x[(x >= fit_low) & (x <= fit_upper) & (y != 0)]

                if len(new_x) == 0:
                    print(f"No data in range for Layer: {layer}, Bar: {strip}, Side: {end}")
                    continue

                # Perform Langau fit
                try:
                    mpv = new_x[np.argmax(new_y)]
                    eta, sigma, A = 1, 1, max(new_y)
                    coeff, _ = curve_fit(pylandau.langau, new_x, new_y, sigma=np.sqrt(new_y), absolute_sigma=True, 
                                         p0=(mpv, eta, sigma, A), bounds=(1, 10000))
                    y_fit = pylandau.langau(new_x, *coeff)
                    chi_square = np.sum(((new_y - y_fit) / np.sqrt(new_y)) ** 2)
                except Exception as e:
                    print(f'Fit failed: {e}')
                    coeff = [-999, -999, -999, -999]
                    chi_square = np.inf

                # Save MIP results
                mip_results.append({
                    'layer': layer,
                    'strip': strip,
                    'end': end,
                    'multiplier': multiplier,
                    'mpv': coeff[0],
                    'eta': coeff[1],
                    'sigma': coeff[2],
                    'A': coeff[3],
                    'chi_square': chi_square
                })

                # Generate plot
                fig, ax = plt.subplots(figsize=(8, 8))
                if coeff[0] != -999:
                    ax.plot(new_x, pylandau.langau(new_x, *coeff), "-", label="Fit")
                ax.errorbar(x, y, yerr=np.sqrt(y), drawstyle='steps-mid', label="Data")
                ax.set_yscale('log')
                ax.set_xlabel(f'Sum ADC {end}')
                ax.set_ylabel('Events')
                ax.set_ylim(1e0, 1e4)
                ax.set_xlim(-50, 2000)
                ax.xaxis.set_major_locator(MultipleLocator(100))
                ax.tick_params(axis='x', rotation=45, labelsize=8)
                ax.legend()
                ax.set_title(f'Layer: {layer} Bar: {strip} Side: {end}')
                pdf.savefig(fig)
                plt.close()

                # Save new fit cut range
                new_fit_cut.append({
                    'layer': layer,
                    'strip': strip,
                    'end': end,
                    'multiplier': multiplier,
                    'low': round(fit_low, 2),
                    'upper': round(fit_upper, 2)
                })

        # Save new fit cut CSV
        new_fit_cut_df = pd.DataFrame(new_fit_cut)
        new_fit_cut_df.to_csv(self.temp_mip_fit_cut_file, index=False)

        # Save MIP results CSV
        mip_df = pd.DataFrame(mip_results)
        mip_df.to_csv(os.path.join(self.out_directory, 'miprevision.csv'), index=False)

if __name__ == "__main__":
    temp_mip_fit_cut_path = 'MIP_revision/mip_fit_cut_for_range.csv'
    data_file_path = 'analysis_files/run_20220425_fpga_run.csv'
    pedestal_file_path = 'calibrations/pedestals.csv'
    out_directory = 'MIP_revision/revision'
    plots_directory = 'MIP_revision/Plot'

    calc_mip_manual = calculateMIPsManual(data_file_path, pedestal_file_path, temp_mip_fit_cut_path, out_directory, plots_directory)
    calc_mip_manual.calculate_mip()
