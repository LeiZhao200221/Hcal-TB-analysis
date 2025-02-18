import shutil
import os
import sys

# Add the src directory to the system path
sys.path.append('src')

# Import the makeAnalysisFiles module
import makeAnalysisFiles

# Define the source file name and target directory
source_file_name = 'source.root'
destination_folder = 'data'

# Ensure the destination folder exists
os.makedirs(destination_folder, exist_ok=True)

# Move the file to the target directory
target_file_name = 'ntuple_decoded_hcal_run_287_20220425_0738127.root'
shutil.move(source_file_name, os.path.join(destination_folder, target_file_name))

print("File successfully moved to the target directory.")

# Initialize the makeAnalysisFiles class
analysis = makeAnalysisFiles.makeAnalysisFiles(
    root_file_name=os.path.join(destination_folder, target_file_name),  # ROOT file
    out_directory='analysis_files',  # Output directory for CSV files
    calibration_file='calibrations/toa_calibration_phase3.csv',  # TOA calibration file
    do_one_bar=False,  # Process only one bar (for debugging)
    do_alignment=False,  # Perform FPGA alignment
    alignment_threshold=20  # FPGA alignment threshold
)

# Run data processing and generate CSV files
analysis.create_dataframes()
