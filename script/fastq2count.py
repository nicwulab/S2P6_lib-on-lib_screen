#!/usr/bin/python
import operator
import glob
from Bio import SeqIO
from collections import Counter
import pandas as pd
import multiprocessing
import time

def hamming(str1, str2):
    assert len(str1) == len(str2)
    return sum(map(operator.ne, str1, str2))

def translation(seq):
  dnamap = {"TTT":"F", "TTC":"F", "TTA":"L", "TTG":"L",
    "TCT":"S", "TCC":"S", "TCA":"S", "TCG":"S",
    "TAT":"Y", "TAC":"Y", "TAA":"_", "TAG":"_",
    "TGT":"C", "TGC":"C", "TGA":"_", "TGG":"W",
    "CTT":"L", "CTC":"L", "CTA":"L", "CTG":"L",
    "CCT":"P", "CCC":"P", "CCA":"P", "CCG":"P",
    "CAT":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
    "CGT":"R", "CGC":"R", "CGA":"R", "CGG":"R",
    "ATT":"I", "ATC":"I", "ATA":"I", "ATG":"M",
    "ACT":"T", "ACC":"T", "ACA":"T", "ACG":"T",
    "AAT":"N", "AAC":"N", "AAA":"K", "AAG":"K",
    "AGT":"S", "AGC":"S", "AGA":"R", "AGG":"R",
    "GTT":"V", "GTC":"V", "GTA":"V", "GTG":"V",
    "GCT":"A", "GCC":"A", "GCA":"A", "GCG":"A",
    "GAT":"D", "GAC":"D", "GAA":"E", "GAG":"E",
    "GGT":"G", "GGC":"G", "GGA":"G", "GGG":"G",}
  pep = []
  i = 0
  while i < len(seq):
    codon = seq[i:i+3]
    aa = dnamap[codon]
    pep.append(aa)
    i = i + 3
  pep = ''.join(pep)
  return pep
  
def ProcessMultilib(Rfile, barcode_list):
  print ("Reading %s" % Rfile)
  records = SeqIO.parse(Rfile,"fastq")
  variants = [] 
  record_count = 0
  for record in records:
    record_count += 1
    Rseq  = record.seq
    Rroi = Rseq
    if ((hamming(Rroi[0:28],"TACGCTCTGCAGGCTAGTGCGTAATAAT") == 0) and (hamming(Rroi[-28:],"CTGATAACAACAGTGTAGATGTAACAAA") == 0) and (len(Rroi) == 108)): 
      # Only include those that have the correct forward primer sequence, correct reverse primer sequence, and the correct number of nucleotides between the primers
      Rroi = Rroi[33:43] + '-' + Rroi[-43:-33] # Trim forward and reverse primers
      #print(Rroi)
      if "N" in Rroi: continue
      #print(Rroi[4:6]+ Rroi[15:17])
      if Rroi[4:6] + Rroi[15:17]!= 'CTGA': continue
      if Rroi not in barcode_list: continue
      variants.append(Rroi)
    #if record_count == 10000: break
  return Counter(variants)

def process_file(fastq_file):
    barcode_file = 'data/barcode_count_info.tsv'
    barcode_df = pd.read_csv(barcode_file, sep='\t')
    barcode_dict = barcode_df.set_index('barcode').T.to_dict('list')
    barcode_list = barcode_dict.keys()
    count_dictionary = ProcessMultilib(fastq_file, barcode_list)
    file_names_df = pd.read_csv('data/file_names.tsv', sep = '\t')
    file_names_dict = file_names_df.set_index('File_name')['Sample_name'].to_dict()
    #print(file_names_dict)
    name_info =  file_names_dict[fastq_file[13:]]
    #print(name_info)
    return name_info, count_dictionary

def main():
  start_time = time.time()
   
  outfile = 'result/barcode_count.tsv'
  fastq_list = glob.glob('fastq_merged/*assembled.fastq')
  
  #print(fastq_list)
  count_df = pd.DataFrame()
  count = 0
  num_processes = 7


  with multiprocessing.Pool(processes=num_processes) as pool:
      results = pool.map(process_file, fastq_list)

  # Initialize an empty DataFrame
  count_df = None

  for name_info, count_dictionary in results:
      if count_df is None:
          count_df = pd.DataFrame.from_records(list(dict(count_dictionary).items()), columns=['muts', name_info + '_count'])
          count_df['muts'] = count_df['muts'].apply(lambda tup: ''.join(map(str, tup)))
      else:
          sub_df = pd.DataFrame.from_records(list(dict(count_dictionary).items()), columns=['muts', name_info + '_count'])
          sub_df['muts'] = sub_df['muts'].apply(lambda tup: ''.join(map(str, tup)))
          count_df = count_df.merge(sub_df, on='muts', how='outer')

  print(count_df)
  cols = ['muts', 
          'Rep1_t0_gate1_count', 'Rep1_t0_gate2_count', 'Rep1_t0_gate3_count', 'Rep1_t0_gate4_count', 
          'Rep1_t1_gate1_count', 'Rep1_t1_gate2_count', 'Rep1_t1_gate3_count', 'Rep1_t1_gate4_count', 
          'Rep1_exp_pos_count', 'Rep1_exp_neg_count', 'Rep2_t0_gate1_count', 'Rep2_t0_gate2_count', 
          'Rep2_t0_gate3_count', 'Rep2_t0_gate4_count','Rep2_t1_gate1_count', 'Rep2_t1_gate2_count', 
          'Rep2_t1_gate3_count', 'Rep2_t1_gate4_count','Rep2_exp_pos_count', 'Rep2_exp_neg_count']
  
  count_df = count_df[cols]
  count_df = count_df.fillna(0)
  count_df.to_csv(outfile, sep="\t", index = False)

  total_time = time.time() - start_time
  print(f"Total processing time: {total_time:.2f} seconds")

if __name__ == "__main__":
  main()
