#!/usr/bin/python
import sys
from Bio import SeqIO
from collections import Counter
import pandas as pd
import numpy as np
from multiprocessing import Pool, cpu_count
import time
import glob

def read_cell_count(infile):
  df = pd.read_csv(infile, sep='\t', header=0)
  return dict(zip(df.iloc[:, 0], df.iloc[:, 1]))

def barcode_to_mut(count_df, barcode_df):
  barcode_dict = barcode_df.set_index('barcode').T.to_dict('list')
  count = 0
  for index, row in count_df.iterrows():
    barcode = count_df['muts'][index]
    if barcode in barcode_dict.keys():
      SH_pep = barcode_dict[barcode][0]
      mut_ID = barcode_dict[barcode][1]
      count_df['SH_pep'][index] = SH_pep
      count_df['mut_ID'][index] = mut_ID
  return count_df

def count_to_freq_col(df, colname):
    df[colname] = pd.to_numeric(df[colname], errors='coerce')
    new_col_name = colname[:-6] + '_freq'
    print('calculate freq for: ' + colname[:-6])
    df[new_col_name] = (df[colname] + 1) / (df[colname].sum() + len(df))
    return df

def apply_count_to_freq_parallel(df, columns):
    with Pool() as pool:
        results = pool.starmap(count_to_freq_col, [(df, col) for col in columns])
    return results

def count_to_total_freq(df, dict_cellcount):
    for rep in ['Rep1', 'Rep2']:
        for t in ['t0', 't1']:
            total_cell = sum([dict_cellcount[rep+'_'+t+'_gate'+str(gate)] for gate in [1,2,3,4]])
            print (rep, t, total_cell)
            df[rep+'_'+t+'_freq'] = df[rep+'_'+t+'_gate1_count']/df[rep+'_'+t+'_gate1_count'].sum()*dict_cellcount[rep+'_'+t+'_gate1'] + \
                                    df[rep+'_'+t+'_gate2_count']/df[rep+'_'+t+'_gate2_count'].sum()*dict_cellcount[rep+'_'+t+'_gate2'] + \
                                    df[rep+'_'+t+'_gate3_count']/df[rep+'_'+t+'_gate3_count'].sum()*dict_cellcount[rep+'_'+t+'_gate3'] + \
                                    df[rep+'_'+t+'_gate4_count']/df[rep+'_'+t+'_gate4_count'].sum()*dict_cellcount[rep+'_'+t+'_gate4']
            df[rep+'_'+t+'_freq'] = df[rep+'_'+t+'_freq']/float(total_cell)
        df[rep+'_totalfreq'] = (df[rep+'_t0_freq']+df[rep+'_t1_freq'])/2
    df['totalfreq'] = (df['Rep1_totalfreq']+df['Rep2_totalfreq'])/2
    return(df)

def read_count_data(df, dict_cellcount):
    colnames = [colname for colname in df if 'mut' not in colname and "SH_pep" not in colname and "freq" not in colname]
    results = apply_count_to_freq_parallel(df, colnames)
    for idx, result_df in enumerate(results):
      column_names = result_df.columns.tolist()
      freq_colname = column_names[-1]
      freq_column = result_df.iloc[:, -1]
      df[freq_colname] = freq_column
    df = count_to_total_freq(df, dict_cellcount)
    cols = ['SH_pep', 'mut_ID',
          'Rep1_t0_gate1_count', 'Rep1_t0_gate2_count', 'Rep1_t0_gate3_count', 'Rep1_t0_gate4_count', 
          'Rep1_t1_gate1_count', 'Rep1_t1_gate2_count', 'Rep1_t1_gate3_count', 'Rep1_t1_gate4_count', 
          'Rep1_exp_pos_count', 'Rep1_exp_neg_count', 'Rep2_t0_gate1_count', 'Rep2_t0_gate2_count', 
          'Rep2_t0_gate3_count', 'Rep2_t0_gate4_count','Rep2_t1_gate1_count', 'Rep2_t1_gate2_count', 
          'Rep2_t1_gate3_count', 'Rep2_t1_gate4_count','Rep2_exp_pos_count', 'Rep2_exp_neg_count',
          'Rep1_t0_gate1_freq', 'Rep1_t0_gate2_freq', 'Rep1_t0_gate3_freq', 'Rep1_t0_gate4_freq',
          'Rep1_t1_gate1_freq', 'Rep1_t1_gate2_freq', 'Rep1_t1_gate3_freq', 'Rep1_t1_gate4_freq', 
          'Rep1_exp_pos_freq', 'Rep1_exp_neg_freq', 'Rep1_t0_freq','Rep1_t1_freq','Rep1_totalfreq',
          'Rep2_t0_gate1_freq', 'Rep2_t0_gate2_freq', 'Rep2_t0_gate3_freq', 'Rep2_t0_gate4_freq',
          'Rep2_t1_gate1_freq', 'Rep2_t1_gate2_freq', 'Rep2_t1_gate3_freq', 'Rep2_t1_gate4_freq',
          'Rep2_exp_pos_freq', 'Rep2_exp_neg_freq', 'Rep2_t0_freq','Rep2_t1_freq','Rep2_totalfreq',
          'totalfreq']
    df = df[cols]
    return df

def main():
  start_time = time.time()
  barcode_file = 'data/barcode_count_info.tsv' 
  outfile = 'result/mut_freq.tsv'
  inputfile = 'result/mut_count.tsv'
  cellfile  = 'data/cell_count.tsv'

  dict_cellcount = read_cell_count(cellfile)

  input_df = pd.read_csv(inputfile, sep='\t')
  #barcode_df = pd.read_csv(barcode_file, sep='\t')

  freq_df = read_count_data(input_df, dict_cellcount)
  print('muts before filter: ', len(freq_df))
  freq_df.to_csv(outfile, sep="\t", index=False)
  total_time = time.time() - start_time
  print(f"Total processing time: {total_time:.2f} seconds")


if __name__ == "__main__":
  main()

