# 1. Decoding.py
## Description
This script moves any file from an arbitrary location to a designated folder and calls Joe’s `makeAnalysisFiles` to decode the NTuple data.

## Notes:
- The target file **must follow the format like**: `ntuple_decoded_hcal_run_287_20220425_0738127.root`, otherwise, Joe’s `makeAnalysisFiles` will not function correctly.
- `do_one_bar` and `do_alignment` are **disabled by default**, but users can modify them as needed.
- The calibration file corresponds to **Phase 3 Test Beam data** from **Run 280 to 315 (April 2022)**. For details on Phase 3, refer to the notebook.

# pedestal_mip.py

## Description  
This script calls Joe’s `calculatePedestals` and `calculateMIPs` modules to compute pedestal values and MIP distributions from NTuple data.  

## Notes  
- The pedestal file **must be generated first** before running the MIP analysis.  
- The pedestal results are saved in `calibrations/pedestals.csv`.  
- The MIP analysis reads `pedestals.csv` and outputs the results to `calibrations/mip_fit_cut.csv`.  
- `do_one_bar` and `plot_mips` are **disabled by default**, but users can modify them as needed.  

---

# mip_fit_cut.py  

## Description  
This script refines the MIP fit range using a **Poisson-based Langau fit** to eliminate multi-muon and electron contamination.  

## Notes  
- The script requires an initial fit range file: `MIP_revision/temp_mip_fit_cut.csv`.  
- This file was manually created after running `pedestal_mip.py` from the **First Glance** folder.  
- The **first peak** in the MIP distribution plots was manually identified to set the `low` and `upper` fit range.  
- The refined fit updates `temp_mip_fit_cut.csv` and saves MIP analysis results in `MIP_revision/miprevision.csv`.  
- Histogram plots of the fits are saved in `MIP_revision/Plot/all_plots.pdf`.  

---

# pedestal_fix.py  

## Description  
This script refines pedestal values by performing a **second Gaussian fit** based on the initial standard deviation obtained from the **First Glance** step.  

## Notes  
- The new fit range is determined using the **standard deviation from the first fit**.  
- The refined pedestal values are saved in `plots/pedestal_plots.pdf`.  
- The script reads pedestal data from `calibrations/pedestals.csv`.  
- The original pedestal range was **too broad**; this refinement ensures better accuracy.  

---

# final_factor.py  

## Description  
This script **corrects pedestal values** by adding a **new ADC selection range** and computes the **muon-calibrated factor values**.  

## Notes  
- A new column **`new_range = pedestal + 15 * std_dev`** is added to `pedestals.csv`.  
- This corrected pedestal file is saved as `calibrations/pedestals_updated.csv`.  
- The muon calibration factor is computed as `factor = 612 / mpv` and saved in `/factor.csv`.  
