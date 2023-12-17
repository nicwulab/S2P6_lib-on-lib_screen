import pandas as pd
import numpy as np

# Read the TSV file into a DataFrame
input_file = 'result/barcode_count.tsv'  # Replace 'your_input_file.tsv' with the actual file path
df = pd.read_csv(input_file, sep='\t')

# Split the DataFrame into 10 smaller DataFrames
split_dfs = np.array_split(df, 10)

# Save each split DataFrame as a separate TSV file
#output_folder = 'output_folder/'  # Replace 'output_folder/' with the desired output folder path

for i, split_df in enumerate(split_dfs):
    output_file = f'split_{i + 1}.tsv'
    split_df.to_csv(output_file, sep='\t', index=False)