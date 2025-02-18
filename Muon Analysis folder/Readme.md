# 1. Decoding.py
## Description
This script moves any file from an arbitrary location to a designated folder and calls Joe’s `makeAnalysisFiles` to decode the NTuple data.

## Notes:
- The target file **must follow the format like**: `ntuple_decoded_hcal_run_287_20220425_0738127.root`, otherwise, Joe’s `makeAnalysisFiles` will not function correctly.
- `do_one_bar` and `do_alignment` are **disabled by default**, but users can modify them as needed.
- The calibration file corresponds to **Phase 3 Test Beam data** from **Run 280 to 315 (April 2022)**. For details on Phase 3, refer to the notebook.

# 2. pedestal_mip.py

## Description  
This script calls Joe’s `calculatePedestals` and `calculateMIPs` modules to compute pedestal values and MIP distributions from NTuple data.  

## Notes  
- The pedestal file **must be generated first** before running the MIP analysis.  
- The pedestal results are saved in `calibrations/pedestals.csv`.  
- The MIP analysis reads `pedestals.csv` and outputs the results to `calibrations/mip_fit_cut.csv`.  
- `do_one_bar` and `plot_mips` are **disabled by default**, but users can modify them as needed.  

---

# 3.mip_fit_cut.py  

## Description  
This script refines the MIP fit range using a **Poisson-based Langau fit** to eliminate multi-muon and electron contamination.  

## Notes  
- The script requires an initial fit range file: `MIP_revision/temp_mip_fit_cut.csv`.  
- This file was manually created after running `pedestal_mip.py` from the **First Glance** folder.  
- The **first peak** in the MIP distribution plots was manually identified to set the `low` and `upper` fit range.  
- The refined fit updates `temp_mip_fit_cut.csv` and saves MIP analysis results in `MIP_revision/miprevision.csv`.  
- Histogram plots of the fits are saved in `MIP_revision/Plot/all_plots.pdf`.  

---

# 4. pedestal_fix.py  

## Description  
This script refines pedestal values by performing another **Gaussian fit** based on the initial standard deviation obtained from the **First Glance** step.  

## Notes  
- The new fit range is determined using the **standard deviation from the first fit**.  
- The refined pedestal values are saved in `plots/pedestal_plots.pdf`.  
- The script reads pedestal data from `calibrations/pedestals.csv`.  
- The original pedestal range was **too broad**; this refinement ensures better accuracy.  

---

# 5.final_factor.py  

## Description  
This script **corrects pedestal values** by adding a **new ADC selection range** and computes the **muon-calibrated factor values**.  

## Notes  
- A new column **`new_range = pedestal + 15 * std_dev`** is added to `pedestals.csv`.  
- This corrected pedestal file is saved as `calibrations/pedestals_updated.csv`.  
- The muon calibration factor is computed as `factor = 612 / mpv` and saved in `/factor.csv`.

---
# 6.ADC_Selection.py

## Description  
This script filters events based on ADC sum values and removing TOT conditions using updated pedestal thresholds. Till now, there is no difference between electron and Muon analysis and this is the last step should be undertaken for First Glance part.

## Notes  
- Reads the pedestal file **`calibrations/pedestals1.csv`**, which contains `new_range` thresholds.  
- Reads the input event data from **`analysis_files/run_20220425_fpga_run.csv`**.  
- Filters events where:  
  - `adc_sum_end0` exceeds `new_range_end0` **or**  
  - `adc_sum_end1` exceeds `new_range_end1`, **and**  
  - Both `tot_end0` and `tot_end1` are **zero**.  
- Saves the filtered results to **`cleaned.csv`**.  
- The final output keeps only `pf_event`, `layer`, `strip`, `adc_sum_end0`, and `adc_sum_end1`.  

---
# 7.energy_calculation.py  (electron selection ends here)

## Description  
This script calculates the calibrated energy values by:  
1. **Subtracting pedestal values** from ADC sums.  
2. **Applying MIP calibration factors** to obtain corrected ADC values.  
3. **Converting ADC values to energy** using a scaling factor.  

## Notes  
- Reads input event data from **`cleaned.csv`**.  
- Uses pedestal values from **`calibrations/pedestals.csv`** to correct ADC sums.  
- Reads MIP calibration values from **`calibrations/mip.csv`**.  
- Uses muon calibration factors from **`factor.csv`**.  
- Filters out events where **energy values (`eng0` or `eng1`) are negative**.  Although ADC selection will be effectively removing all negative value here, just in case, an extra check is executed here.
- Saves the final calibrated energy data in **`cal_data.csv`** for further use. (like, Layer energy deposit Pattern/ Strip energy deposit pattern).
- There is a factor 612 which I select as standard for calibration. This value is close to the mean/median value of data in Run 287, 04/2022. It should be carefully selected by different value for different period. Both value here and 612 in section 5 should be changed when data from different period are introduced.
---
#8.muon_selection.py (a further selection should be done for muon)

## Description  
This script applies a **Muon-specific selection** to ensure event purity by:  
1. **Filtering events where only one strip per layer is hit** to eliminate electron contamination and multiple-muon events.  
2. **Applying an energy selection (1.84 - 7.83 MIP energy)** to retain only single-strip, single-hit, single-muon events.  

## Notes  
- Muons **should not interact significantly** in the HCal, so single-strip responses ensure **pure muon selection**.  
- Multiple muons hitting the same strip may deposit higher energy, so an energy window of **0.4 - 1.6 MIP (1.84 - 7.83 absolute energy)** is applied.  
- The script reads input data from **`cal_data.csv`**.  
- The final selected data is saved in **`cal_data_muon.csv`**.  

---
All Template for calibration factor/results I used/calculated could be found in calibration folder

