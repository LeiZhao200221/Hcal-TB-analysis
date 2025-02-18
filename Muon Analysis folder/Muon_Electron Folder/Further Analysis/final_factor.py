import pandas as pd

# File paths
mip_file_path = 'MIP_revision/revision/mip.csv'  # MIP file
pedestal_file_path = 'calibrations/pedestals.csv'  # Old pedestal file
output_factor_file = 'factor.csv'  # Output factor file
output_pedestal_file = 'calibrations/pedestals_updated.csv'  # Updated pedestal file

# Step 1: Compute factor and save new CSV
# Read MIP data
mip_data = pd.read_csv(mip_file_path)

# Calculate new factor column (rounded to 4 decimal places)
mip_data['factor'] = (612 / mip_data['mpv']).round(4)

# Keep only required columns
factor_output_data = mip_data[['layer', 'strip', 'end', 'mpv', 'factor']]

# Save factor results
factor_output_data.to_csv(output_factor_file, index=False)
print(f"Factor file saved to: {output_factor_file}")

# Step 2: Update pedestal file with new_range
# Read pedestal data
pedestal_data = pd.read_csv(pedestal_file_path)

# Compute new_range (pedestal + 15 * std_dev)
pedestal_data['new_range'] = pedestal_data['pedestal'] + 15 * pedestal_data['std_dev']

# Save updated pedestal data
pedestal_data.to_csv(output_pedestal_file, index=False)
print(f"Updated pedestal file saved to: {output_pedestal_file}")
