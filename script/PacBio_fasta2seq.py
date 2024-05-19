#!/usr/bin/python
import string
from Bio import SeqIO
from collections import Counter, defaultdict

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
    "GGT":"G", "GGC":"G", "GGA":"G", "GGG":"G",
    "NNK":"X"}
  pep = []
  i = 0
  while i < len(seq):
    codon = seq[i:i+3]
    aa = dnamap[codon]
    pep.append(aa)
    i = i + 3
  pep = ''.join(pep)
  return pep

def rc(seq):
  seq = str(seq)
  complements = str.maketrans('acgtrymkbdhvACGTRYMKBDHV', 'tgcayrkmvhdbTGCAYRKMVHDB')
  rcseq = seq.translate(complements)[::-1]
  return rcseq

def read_to_seq(record, barcode_spacer):
  if barcode_spacer in record.seq:
    return str(record.seq)
  elif rc(barcode_spacer) in record.seq:
    return str(rc(record.seq))
  else:
    return 'FAIL'

def QC(seq, features):
  for feature in features:
    if seq.count(feature) != 1:
      return 'FAIL'
  return 'GOOD'

def seq_to_barcode(seq, barcode_spacer):
  seq_array = seq.rsplit(barcode_spacer)
  barcode1 = seq_array[0][-10::]
  barcode2 = seq_array[1][0:10]
  if barcode1[4:6] == 'CT' or barcode2[4:6] == 'GA': 
    return barcode1+'-'+barcode2
  else:
    return 'FAIL'
  
def seq_to_feature(seq, feature_5, feature_3, feature_len):
  if feature_5 == 'SH_helix':
    feature = seq.rsplit(feature_3)[0][-42::]
  else:
    feature = seq.rsplit(feature_5)[1].rsplit(feature_3)[0]
  if len(feature) == feature_len:
    return translation(feature)
  else:
    return 'FAIL'

def extract_info_from_fasta(fasta_file, barcode_spacer, AAQPA_linker, post_SH_linker, gs_linker, gibson, outfile):
  print ('writing: %s' % outfile)
  outfile = open(outfile, 'w')
  outfile.write("\t".join(['barcode', 'SH_pep', 'mutID'])+"\n")
  count_record = 0
  for record in SeqIO.parse(fasta_file, "fasta"):
    count_record += 1
    if count_record % 200000==0:
      print ('processed %s reads' % count_record)
    seq = read_to_seq(record, barcode_spacer)
    if seq == 'FAIL': continue
    QC_result = QC(seq, [barcode_spacer, AAQPA_linker, post_SH_linker, gs_linker, gibson])
    if QC_result == 'FAIL': continue
    barcode = seq_to_barcode(seq, barcode_spacer)
    if barcode == 'FAIL': continue
    SH_seq = seq_to_feature(seq, 'SH_helix', post_SH_linker, 42)
    if SH_seq == 'FAIL': continue
    VH_pep = seq_to_feature(seq, AAQPA_linker, gs_linker, 354)
    if VH_pep == 'FAIL': continue
    VL_pep = seq_to_feature(seq, gs_linker, gibson, 330)
    if VL_pep == 'FAIL': continue
    mutID = VH_pep[31:32]+VH_pep[47:48]+VH_pep[55:57]+VH_pep[69:70]+VH_pep[78:79]+'-'+VL_pep[3:5]+VL_pep[29:30]+VL_pep[31:32]
    outfile.write("\t".join(map(str,[barcode, SH_seq, mutID]))+"\n")
  outfile.close()

def main():
  outfile = 'data/barcode_SHpep_mutID.tsv'
  fasta_file = 'data/S2P6_lib_amplicon-for-rev.fasta'
  barcode_spacer = 'CTACGTTCAAGGCTATGCATCG'
  AAQPA_linker   = 'GCGGCCCAGCCGGCC'
  gs_linker = 'GGCGGAGGTGGGAGTGGAGGAGGCGGTTCTGGTTCAGGTGGCGGAGGTTCC'
  post_SH_linker = 'GAATTCGGCGGAGGTGGG'
  gibson = 'GGCCTCGGGGGCCTGTACCC'
  extract_info_from_fasta(fasta_file, barcode_spacer, AAQPA_linker, post_SH_linker, gs_linker, gibson, outfile)

if __name__ == "__main__":
  main()
