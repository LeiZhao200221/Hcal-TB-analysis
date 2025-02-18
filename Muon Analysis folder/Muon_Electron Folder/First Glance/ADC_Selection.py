import pandas as pd

# Define file paths
input_file_path = "analysis_files/run_20220425_fpga_run.csv"
pedestal_file_path = "calibrations/pedestals.csv"
output_folder = "Datafortable"
output_file = f"{output_folder}/cleaned.csv"

# Read pedestal file
pedestal_df = pd.read_csv(pedestal_file_path)

# Split pedestal data into end0 and end1 parts
pedestal_end0_df = pedestal_df[pedestal_df["end"] == 0][["layer", "strip", "new_range"]].rename(
    columns={"new_range": "new_range_end0"}
)
pedestal_end1_df = pedestal_df[pedestal_df["end"] == 1][["layer", "strip", "new_range"]].rename(
    columns={"new_range": "new_range_end1"}
)

# Define chunk size
chunksize = 100000

# Initialize an empty list to store filtered data chunks
filtered_chunks = []

# Read and process data in chunks
for chunk in pd.read_csv(input_file_path, chunksize=chunksize):
    # Merge chunk with pedestal_end0_df and pedestal_end1_df
    merged_chunk = chunk.merge(pedestal_end0_df, on=["layer", "strip"], how="left")
    merged_chunk = merged_chunk.merge(pedestal_end1_df, on=["layer", "strip"], how="left")

    # Filter data based on new criteria
    filtered_chunk = merged_chunk[
        (merged_chunk["adc_sum_end0"] > merged_chunk["new_range_end0"])
        | (merged_chunk["adc_sum_end1"] > merged_chunk["new_range_end1"])
        & (merged_chunk["tot_end0"] == 0)
        & (merged_chunk["tot_end1"] == 0)
    ]

    # Append filtered chunk to list
    filtered_chunks.append(filtered_chunk)

# Concatenate all filtered data chunks
filtered_df = pd.concat(filtered_chunks)

# Sort the filtered data by pf_event, layer, and strip
sorted_df = filtered_df.sort_values(by=["pf_event", "layer", "strip"])

# Keep only the required columns
final_df = sorted_df[["pf_event", "layer", "strip", "adc_sum_end0", "adc_sum_end1"]]

# Save the filtered data to a new CSV file
final_df.to_csv(output_file, index=False)

print(f"Data has been saved to: {output_file}")
