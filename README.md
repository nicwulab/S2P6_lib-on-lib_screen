## Yeast Library to Library Screening (Beta Coronavirus Stem Helix vs. S2P6)

### Merge Illumina seqeuncing reads
``python3 script/merge_reads.py``<br />
    - Input file:<br />
      - All .fastq files in [fastq/]<br />
    - Output files:<br />
      - merged files in [fastq/merged]<br />


### Counting unique barcode sequences and keeping the ones in PacBio reference file
``python3 script/fastq2count.py``<br />
    - Input file:<br />
      - All merged fastq files in [fastq/merged]<br />
      - [./data/file_names.tsv](./data/file_names.tsv)<br />
    - Output files:<br />
      - [./result/barcode_count.tsv](./result/barcode_count.tsv)<br />

### Split count file for faster processing
``python3 script/split_count_df.py``<br />
    - Input file:<br />
      - [./result/barcode_count.tsv](./result/barcode_count.tsv)<br />
    - Output files:<br />
      - [./result/split_files/split_*.tsv](./result/split_files)<br />

### Indentify SH pep and S2P6 mutant pair
``python3 script/mut_ID.py``<br />
    - Input file:<br />
      - [./result/split_files/split_*.tsv](./result/split_files)<br />
      - [./data/barcode_count_info.tsv](./data/barcode_count_info.tsv)<br />
    - Output files:<br />
      - [./result/mut_count.tsv](./result/mut_count.tsv)<br />

### Calculate the frequency, then apply the input frequncy filter and recalculate frequency
``python3 script/count2freq.py``<br />
    - Input file:<br />
      - [./result/mut_count.tsv](./result/mut_count.tsv)<br />
    - Output files:<br />
      - [./result/mut_freq.tsv](./result/mut_freq.tsv)<br />
      - [./result/mut_freq_1*10-5.tsv](./result/mut_freq_1*10-5.tsv)<br /> 
      - [./result/mut_freq_2*10-5.tsv](./result/mut_freq_2*10-5.tsv)<br /> 

### Calculate the expression score, expression pos/neg, and binding score
``python3 script/freq2score.py``<br />
    - Input file:<br />
      - [./result/mut_freq.tsv](./result/mut_freq.tsv)<br />
      - [./result/mut_freq_1*10-5.tsv](./result/mut_freq_1*10-5.tsv)<br /> 
      - [./result/mut_freq_2*10-5.tsv](./result/mut_freq_2*10-5.tsv)<br /> 
    - Output files:<br />
      - [./result/mut_scores.tsv](./result/mut_scores.tsv)<br />
      - [./result/mut_scores_1*10-5.tsv](./result/mut_scores_1*10-5.tsv)<br /> 
      - [./result/mut_scores_2*10-5.tsv](./result/mut_scores_2*10-5.tsv)<br /> 

### Plot correlation between expression score and binding scores
``python3 script/plot_replicate_qc.py``<br />
    - Input file:<br />
      - [./result/mut_scores.tsv](./result/mut_scores.tsv)<br />
      - [./result/mut_scores_1*10-5.tsv](./result/mut_scores_1*10-5.tsv)<br /> 
      - [./result/mut_scores_2*10-5.tsv](./result/mut_scores_2*10-5.tsv)<br /> 
    - Output files:<br />
      - [./graph/QC/Binding_score_correlation.html](./graph/QC/Binding_score_correlation.html)<br />
      - [./graph/QC/Binding_score_correlation_1*10-5.html](./graph/QC/Binding_score_correlation_1*10-5.html)<br />
      - [./graph/QC/Binding_score_correlation_2*10-5.html](./graph/QC/Binding_score_correlation_2*10-5.html)<br />
      - [./graph/QC/Expression_score_correlation.html](./graph/QC/Expression_score_correlation.html)<br />
      - [./graph/QC/Expression_score_correlation_1*10-5.html](./graph/QC/Expression_score_correlation_1*10-5.html)<br />
      - [./graph/QC/Expression_score_correlation_2*10-5.html](./graph/QC/Expression_score_correlation_2*10-5.html)<br />

      

      
