# 1. Decoding.py
## Description
This script moves any file from an arbitrary location to a designated folder and calls Joe’s `makeAnalysisFiles` to decode the NTuple data.

## Notes:
- The target file **must follow the format like**: `ntuple_decoded_hcal_run_287_20220425_0738127.root`, otherwise, Joe’s `makeAnalysisFiles` will not function correctly.
- `do_one_bar` and `do_alignment` are **disabled by default**, but users can modify them as needed.
- The calibration file corresponds to **Phase 3 Test Beam data** from **Run 280 to 315 (April 2022)**. For details on Phase 3, refer to the notebook.

# Pedestal and MIP Calculation (`pedestal_mip.py`)

## Description
This script calls Joe’s **pedestal analysis module** and **MIP analysis module** to compute pedestal values and MIP distributions from the NTuple data.

## Process:
1. **Compute Pedestal:**  
   - Uses `calculatePedestals` to analyze the ROOT file and generate pedestal values.
   - The pedestal results are saved as `calibrations/pedestals.csv`.

2. **Compute MIP:**  
   - Uses `calculateMIPs` to process the pedestal-corrected data.
   - Reads the pedestal file (`pedestals.csv`) from the previous step.
   - Outputs MIP analysis results to `calibrations/mip_fit_cut.csv`.

## Notes:
- The **pedestal file is required** for the MIP analysis.
- Users can modify parameters such as `do_one_bar` and `plot_mips` as needed.


# MIP Fit Cut (`mip_fit_cut.py`)

## Description
This script processes MIP data by refining the fit range and performing a Poisson-based Langau fit to eliminate multi-muon and electron contamination.

## Data Source
- The script requires an **initial fit range file**:  
  `MIP_revision/temp_mip_fit_cut.csv`
- This file was created after **running `pedestal_mip.py`** from the `First Glance` folder.
- The initial fit range (lower and upper bounds) was **manually identified** based on the first peak in the MIP distribution plots.

## Processing Steps
1. **Pedestal Subtraction**:  
   - Reads `calibrations/pedestals.csv` to correct ADC values.
2. **Filtering & Selection**:  
   - Uses the manually set `low` and `upper` fit range from `temp_mip_fit_cut.csv`.(Hope someone could develop a code to solve this)
   - Removes **multi-muon and electron contamination**.
3. **Poisson-Based Langau Fit**:  
   - Performs a refined fit within the selected range.
   - Updates `temp_mip_fit_cut.csv` with the new refined fit range.
4. **Result Output**:  
   - Saves the updated MIP analysis to `MIP_revision/miprevision.csv`.
   - Generates **histogram plots** of the fits in `MIP_revision/Plot/all_plots.pdf`.

## Notes
- Users can adjust `temp_mip_fit_cut.csv` manually if necessary.
- The refined fit ensures better isolation of the **single-muon MIP peak**.

