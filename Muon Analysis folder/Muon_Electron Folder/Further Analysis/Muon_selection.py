import pandas as pd

# Define file paths
file_path = "cal_data.csv"
output_file_path = "cal_data_muon.csv"

# Load data
data = pd.read_csv(file_path)

# Step 1: Ensure only one strip per event per layer
layer_strip_counts = data.groupby(["pf_event", "layer"])["strip"].count()
valid_events = layer_strip_counts[layer_strip_counts == 1].index.get_level_values("pf_event")

# Filter the data to keep only valid events
filtered_data = data[data["pf_event"].isin(valid_events)]

# Step 2: Apply energy selection for layer 1
layer1_data = filtered_data[filtered_data["layer"] == 1]
valid_eng = layer1_data[(layer1_data["eng0"].between(1.84, 7.83)) & (layer1_data["eng1"].between(1.84, 7.83))]

# Get final valid event IDs
valid_final_events = valid_eng["pf_event"].unique()

# Select final valid data
final_filtered_data = filtered_data[filtered_data["pf_event"].isin(valid_final_events)]

# Save results
final_filtered_data.to_csv(output_file_path, index=False)

print(f"Filtered data saved to {output_file_path}")
