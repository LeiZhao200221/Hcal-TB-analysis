import sys
import os

# Add the src directory to the system path
sys.path.append('src')

# Import necessary modules
from calculatePedestals import calculatePedestals
from calculateMIPs import calculateMIPs

# Step 1: Compute Pedestal
pedestal_calculator = calculatePedestals(
    root_file_name='data/ntuple_decoded_hcal_run_287_20220425_0738127.root',  # ROOT file
    out_directory='calibrations',  # Output directory for pedestal files
    plot_pedestals=False,  # Whether to generate pedestal plots
    plots_directory='plots/pedestals',  # Directory for pedestal plots
    do_one_bar=False  # Process only one bar (for debugging)
)

# Run pedestal calculation and generate CSV file
pedestal_calculator.get_pedestals()

# Step 2: Compute MIP
mip_calculator = calculateMIPs(
    data_file_name='analysis_files/run_20220425_fpga_run.csv',  # CSV file generated in step 1
    pedestal_file_name='calibrations/pedestals.csv',  # Pedestal CSV file generated in step 2
    mip_fit_cut_file_name='calibrations/mip_fit_cut.csv',  # MIP fit cut file
    out_directory='calibrations',  # Output directory for MIP files
    plot_mips=True,  # Whether to generate MIP plots
    plots_directory='plots/mips',  # Directory for MIP plots
    do_one_bar=False,  # Process only one bar (for debugging)
    calc_multipliers=False  # Whether to manually calculate multipliers
)

# Run MIP calculation and generate CSV file
mip_calculator.get_mips()
