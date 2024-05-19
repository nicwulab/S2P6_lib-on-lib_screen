#!/usr/bin/python
import string
from Bio import SeqIO
from collections import Counter, defaultdict

def EFseq(readlist,Ecutoff,Scutoff):
  totalread = len(readlist)
  if totalread < Scutoff: return 'bad'
  realread = ''
  for j in range(0,len(readlist[0])):
    bases = {}
    for read in readlist:
      base = read[j]
      if base in bases: bases[base] += 1
      else: bases[base] = 1
    check = 'bad'
    bs = ''
    for b in bases.keys():
      if float(bases[b])/float(totalread) >= Ecutoff:
        bs += b
    if len(bs) != 1: return 'bad'
    else:
      realread += bs
      check = 'good'
    if check == 'bad':
      return 'bad'
  return realread

def QC_mutID(mutID):
  if mutID[0] not in ['Y','Q']: return 'fail'
  elif mutID[1] not in ['M','I']: return 'fail'
  elif mutID[2] not in ['G','V']: return 'fail'
  elif mutID[3] not in ['S','H']: return 'fail'
  elif mutID[4] not in ['M','L']: return 'fail'
  elif mutID[5] not in ['V','L']: return 'fail'
  elif mutID[7] not in ['L','M']: return 'fail'
  elif mutID[8] not in ['T','M']: return 'fail'
  elif mutID[9] not in ['S','R']: return 'fail'
  elif mutID[10] not in ['S','N']: return 'fail'
  else: return 'pass'

def clustering(infile, SH_pep_list):
  clusters = {}
  infile = open(infile, 'r')
  count_line = 0
  for line in infile.readlines():
    count_line += 1
    if count_line % 200000==0:
      print ('processed %s reads' % count_line)
    if 'barcode' in line:
      continue
    barcode, SH_pep, mutID = line.rstrip().rsplit("\t")
    if SH_pep not in SH_pep_list: 
      continue
    if QC_mutID(mutID) == 'fail': 
      continue
    try:
      clusters[barcode]['SH_pep'].append(SH_pep)
      clusters[barcode]['mutID'].append(mutID)
    except:
      clusters[barcode] = {'SH_pep':[SH_pep],
                           'mutID':[mutID]}
  return clusters
 
def compress_cluster(clusters, outfile):
  print ('writing: %s' % outfile)
  outfile = open(outfile, 'w')
  outfile.write("\t".join(['barcode', 'count_reads', 'SH_pep', 'mutID'])+"\n")
  for barcode in clusters.keys():
    count_reads = len(clusters[barcode]['SH_pep'])
    SH_pep_list = clusters[barcode]['SH_pep']
    mutID_list  = clusters[barcode]['mutID']
    SH_pep = EFseq(SH_pep_list,0.8,2)
    mutID  = EFseq(mutID_list,0.8,2)
    if SH_pep != 'bad' and mutID != 'bad':
      outfile.write("\t".join(map(str,[barcode, count_reads, SH_pep, mutID]))+"\n")
  outfile.close()

def main():
  outfile = 'data/barcode_count_info.tsv'
  infile  = 'data/barcode_SHpep_mutID.tsv'
  SH_pep_dict = SeqIO.to_dict(SeqIO.parse('Fasta/SH_pep_ref.fa', "fasta"))
  SH_pep_list = list(map(lambda x:str(x.seq), SH_pep_dict.values()))
  print (SH_pep_list)
  clusters = clustering(infile, SH_pep_list)
  compress_cluster(clusters, outfile)

if __name__ == "__main__":
  main()
