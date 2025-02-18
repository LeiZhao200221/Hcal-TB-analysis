import pandas as pd
import os

# Define file paths
pedestals_path = "calibrations/pedestals.csv"
pre_final_path = "Datafortable/cleaned.csv"
mip_path = "calibrations/mip.csv"
factor_path = "factor.csv"
output_csv = "cal_data.csv"

# Load data
pedestals_df = pd.read_csv(pedestals_path)
pre_final_df = pd.read_csv(pre_final_path)
mip_df = pd.read_csv(mip_path)
factor_df = pd.read_csv(factor_path)

# Ensure ADC values are floats
pre_final_df["adc_sum_end0"] = pre_final_df["adc_sum_end0"].astype(float)
pre_final_df["adc_sum_end1"] = pre_final_df["adc_sum_end1"].astype(float)

# Convert pedestal and MIP data to dictionary for fast lookup
pedestals_dict = pedestals_df.set_index(["layer", "strip", "end"]).to_dict(orient="index")
mip_dict = mip_df.set_index(["layer", "strip", "end"]).to_dict(orient="index")

# Function to compute pedestal-subtracted ADC and MIP values
def calculate_mip(adc_sum, layer, strip, end):
    pedestal_key = (layer, strip, end)
    mip_key = (layer, strip, end)

    if pedestal_key in pedestals_dict and mip_key in mip_dict:
        pedestal = pedestals_dict[pedestal_key]["pedestal"]
        adc_minus_pedestal = adc_sum - pedestal
        mpv = mip_dict[mip_key]["mpv"]
        if mpv != 999:
            mip_value = round(adc_minus_pedestal / mpv, 2)
        else:
            mip_value = None
    else:
        adc_minus_pedestal = None
        mip_value = None

    return adc_minus_pedestal, mip_value

# Process data
output_data = []

for _, row in pre_final_df.iterrows():
    pf_event = row["pf_event"]
    layer = row["layer"]
    strip = row["strip"]

    adc_minus_pedestal_end0, mip_value_end0 = calculate_mip(row["adc_sum_end0"], layer, strip, 0)
    adc_minus_pedestal_end1, mip_value_end1 = calculate_mip(row["adc_sum_end1"], layer, strip, 1)

    output_data.append([pf_event, layer, strip, adc_minus_pedestal_end0, adc_minus_pedestal_end1])

# Convert to DataFrame
output_columns = ["pf_event", "layer", "strip", "adc_minus_pedestal_end0", "adc_minus_pedestal_end1"]
df_jian_suo = pd.DataFrame(output_data, columns=output_columns)

# Merge calibration factor
df_jian_suo = pd.merge(
    df_jian_suo, factor_df[factor_df["end"] == 0][["layer", "strip", "factor"]],
    on=["layer", "strip"], how="left"
).rename(columns={"factor": "factor_end0"})

df_jian_suo = pd.merge(
    df_jian_suo, factor_df[factor_df["end"] == 1][["layer", "strip", "factor"]],
    on=["layer", "strip"], how="left"
).rename(columns={"factor": "factor_end1"})

# Apply calibration factors
df_jian_suo["cali_adc_minus_pede_end0"] = df_jian_suo["adc_minus_pedestal_end0"] * df_jian_suo["factor_end0"]
df_jian_suo["cali_adc_minus_pede_end1"] = df_jian_suo["adc_minus_pedestal_end1"] * df_jian_suo["factor_end1"]

# Compute energy values
df_jian_suo["eng0"] = df_jian_suo["cali_adc_minus_pede_end0"] * (4.83143 / 612)
df_jian_suo["eng1"] = df_jian_suo["cali_adc_minus_pede_end1"] * (4.83143 / 612)

# Drop unnecessary columns
df_jian_suo = df_jian_suo.drop(columns=["factor_end0", "factor_end1"])

# Convert event, layer, strip to integer
df_jian_suo["pf_event"] = df_jian_suo["pf_event"].astype(int)
df_jian_suo["layer"] = df_jian_suo["layer"].astype(int)
df_jian_suo["strip"] = df_jian_suo["strip"].astype(int)

# Remove negative energy values
df_jian_suo = df_jian_suo[(df_jian_suo["eng0"] >= 0) & (df_jian_suo["eng1"] >= 0)]

# Save final energy data
df_jian_suo.to_csv(output_csv, index=False, float_format="%.6f")

print(f"Energy data saved to {output_csv}")
