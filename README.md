## Library-on-library screen of betacoronavirus stem helix peptides vs S2P6 mutants

### Obejctive
Studying the evolutionary trajectories of S2P6 for breadth expansion using a library-on-library screen that involves 27 unique betacoronavirus stem helix peptides and 1,024 S2P6 variants.

### Dependencies
* [Python](https://www.python.org/) (version 3.9)
* [Biopython](https://github.com/biopython/biopython)
* [R](https://www.r-project.org/) (version 4.1)
* [PEAR](https://github.com/tseemann/PEAR)

### Input files
* [./Fasta/SH_pep_ref.fa](./Fasta/SH_pep_ref.fa): A list of stem helix peptides in this experiment
* [./data/file_names.tsv](./data/file_names.tsv): Filenames for the merged read files
* Raw PacBio seqeucing data in fasta format from NIH SRA database [BioProject PRJNAXXXXXX](https://www.ncbi.nlm.nih.gov/bioproject/PRJNAXXXXXX)
* Raw Illumina seqeucing data in fastq format from NIH SRA database [BioProject PRJNAXXXXXX](https://www.ncbi.nlm.nih.gov/bioproject/PRJNAXXXXXX)

### Linking barcodes to variants based on PacBio sequencing data
1. Identify the sequences of stem helix peptide, S2P6 variant, and barcode in each read   
``python3 script/PacBio_fasta2seq.py``<br />
    - Input file:<br />
      - Fasta file from the PacBio sequencing
    - Output file:<br />
      - data/barcode\_SHpep\_mutID.tsv

2. Filter barcodes with low read counts and perform error correction   
    - Input file: 
      - data/barcode\_SHpep\_mutID.tsv
      - [./Fasta/SH_pep_ref.fa](./Fasta/SH_pep_ref.fa)
    - Ouput file:
      - [./data/barcode_count_info.tsv](./data/barcode_count_info.tsv)<br />

### Analyze the bacode sequencing data for the library-on-library screen
1. Merging Illumina seqeuncing reads   
``python3 script/merge_reads.py``<br />
    - Input file:<br />
      - All .fastq files in [fastq/]<br />
    - Output files:<br />
      - merged files in [fastq/merged]<br />

2. Counting unique barcode sequences   
``python3 script/fastq2count.py``<br />
    - Input file:<br />
      - All merged fastq files in [fastq/merged]<br />
      - [./data/file_names.tsv](./data/file_names.tsv)<br />
    - Output files:<br />
      - [./result/barcode_count.tsv](./result/barcode_count.tsv)<br />

3. Splitting count file for faster processing   
``python3 script/split_count_df.py``<br />
    - Input file:<br />
      - [./result/barcode_count.tsv](./result/barcode_count.tsv)<br />
    - Output files:<br />
      - [./result/split_files/split_\*.tsv](./result/split_files)<br />

4. Indentifying pairs of stem helix peptide and S2P6 mutant   
``python3 script/mut_ID.py``<br />
    - Input file:<br />
      - [./result/split_files/split_\*.tsv](./result/split_files)<br />
      - [./data/barcode_count_info.tsv](./data/barcode_count_info.tsv)<br />
    - Output files:<br />
      - [./result/mut_count.tsv](./result/mut_count.tsv)<br />

5. Calculating the frequency of each variant   
``python3 script/count2freq.py``<br />
    - Input file:<br />
      - [./result/mut_count.tsv](./result/mut_count.tsv)<br />
    - Output files:<br />
      - [./result/mut_freq.tsv](./result/mut_freq.tsv)<br />

6. Calculate the expression scores and binding scores   
``python3 script/freq2score.py``<br />
    - Input file:<br />
      - [./result/mut_freq.tsv](./result/mut_freq.tsv)<br />
    - Output files:<br />
      - [./result/mut_scores.tsv](./result/mut_scores.tsv)<br />

### Plotting   
1. Plot correlation between expression scores and binding scores   
``python3 script/plot_replicate_qc.py``<br />
    - Input file:<br />
      - [./result/mut_scores.tsv](./result/mut_scores.tsv)<br />
    - Output files:<br />
      - [./graph/QC/Binding_score_correlation.html](./graph/QC/Binding_score_correlation.html)<br />
      - [./graph/QC/Binding_score_correlation_1\*10-5.html](./graph/QC/Binding_score_correlation_1*10-5.html)<br />
      - [./graph/QC/Binding_score_correlation_2\*10-5.html](./graph/QC/Binding_score_correlation_2*10-5.html)<br />
      - [./graph/QC/Expression_score_correlation.html](./graph/QC/Expression_score_correlation.html)<br />
      - [./graph/QC/Expression_score_correlation_1\*10-5.html](./graph/QC/Expression_score_correlation_1*10-5.html)<br />
      - [./graph/QC/Expression_score_correlation_2\*10-5.html](./graph/QC/Expression_score_correlation_2*10-5.html)<br />

2. Plot correlation between binding scores and effect of different frequency cutoffs   
``Rscript script/plot_QC.R``<br />
    - Input file:<br />
      - [./result/mut_scores.tsv](./result/mut_scores.tsv)<br />
    - Output files:<br />
      - [./graph/QC_cutoff_freq_vs_cor.png](./graph/QC_cutoff_freq_vs_cor.png)<br />
      - [./graph/QC_cutoff_freq_vs_variant_num.png](./graph/QC_cutoff_freq_vs_variant_num.png)<br />
      - [./graph/QC_replicate_cor.png](./graph/QC_replicate_cor.png)<br />
