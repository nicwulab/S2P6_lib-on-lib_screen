import dask.dataframe as dd
from dask.diagnostics import ProgressBar
from multiprocessing import cpu_count
import pandas as pd
import time
import multiprocessing
import glob

def process_df(list_of_df, barcode_df):
    print("start process dfs")
    barcode_dict = barcode_df.set_index('barcode').T.to_dict('list')
    print(list_of_df)
    with multiprocessing.Pool() as pool:
        processed_dfs = pool.starmap(barcode_to_mut, [(df, barcode_dict) for df in list_of_df])
    combined_df = combine_dataframes(processed_dfs)
    return combined_df

def combine_dataframes(processed_dfs):
    print("start combine dfs")
    # Concatenate all DataFrames into a single DataFrame
    combined_df = pd.concat(processed_dfs, ignore_index=True)
    print(combined_df)
    # Group by 'SH_pep' and 'mut_ID' and sum other columns
    combined_df = combined_df.drop(columns=['muts'])  # Drop 'muts' column
    combined_df = combined_df.groupby(['SH_pep', 'mut_ID']).sum().reset_index()
    print(combined_df)
    return combined_df

def barcode_to_mut(count_df_name, barcode_dict):
  print("start " + count_df_name)
  count_df = pd.read_csv(count_df_name, sep='\t')
  count_df['SH_pep'] = ''
  count_df['mut_ID'] = ''
  for index, row in count_df.iterrows():
    barcode = count_df['muts'][index]
    if barcode in barcode_dict.keys():
      SH_pep = barcode_dict[barcode][1]
      mut_ID = barcode_dict[barcode][2]
      count_df.at[index, 'SH_pep'] = SH_pep
      count_df.at[index, 'mut_ID'] = mut_ID
  print(count_df)
  return count_df

    
def main():
  start_time = time.time()
  barcode_file = 'data/barcode_count_info.tsv' 
  outfile = 'result/mut_count.tsv'
  list_of_df = glob.glob('result/split_files/*.tsv')
  barcode_df = pd.read_csv(barcode_file, sep='\t')
  print('start processing..')
  outputdf = process_df(list_of_df, barcode_df)
  cols = ['SH_pep', 'mut_ID',
          'Rep1_t0_gate1_count', 'Rep1_t0_gate2_count', 'Rep1_t0_gate3_count', 'Rep1_t0_gate4_count', 
          'Rep1_t1_gate1_count', 'Rep1_t1_gate2_count', 'Rep1_t1_gate3_count', 'Rep1_t1_gate4_count', 
          'Rep1_exp_pos_count', 'Rep1_exp_neg_count', 'Rep2_t0_gate1_count', 'Rep2_t0_gate2_count', 
          'Rep2_t0_gate3_count', 'Rep2_t0_gate4_count','Rep2_t1_gate1_count', 'Rep2_t1_gate2_count', 
          'Rep2_t1_gate3_count', 'Rep2_t1_gate4_count','Rep2_exp_pos_count', 'Rep2_exp_neg_count']
  outputdf = outputdf[cols]
  outputdf.to_csv(outfile, sep="\t", index=False)
  total_time = time.time() - start_time
  print(f"Total processing time: {total_time:.2f} seconds")


if __name__ == "__main__":
  main()
