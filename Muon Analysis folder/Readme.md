# 1. Decoding.py
## Description
This script moves any file from an arbitrary location to a designated folder and calls Joe’s `makeAnalysisFiles` to decode the NTuple data.

## Notes:
- The target file **must follow the format like**: `ntuple_decoded_hcal_run_287_20220425_0738127.root`, otherwise, Joe’s `makeAnalysisFiles` will not function correctly.
- `do_one_bar` and `do_alignment` are **disabled by default**, but users can modify them as needed.
- The calibration file corresponds to **Phase 3 Test Beam data** from **Run 280 to 315 (April 2022)**. For details on Phase 3, refer to the notebook.

