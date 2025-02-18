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
