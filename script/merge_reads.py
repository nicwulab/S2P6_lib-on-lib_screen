import os
import glob
from multiprocessing import Pool, cpu_count

def merge_files(R1_file):
    R2_file = R1_file.replace('_R1_', '_R2_')
    outfile = os.path.join('fastq_merged', os.path.basename(R1_file).replace('_R1_', 'merged'))
    assert(R1_file != R2_file)
    os.system(f'pear -f {R1_file} -r {R2_file} -o {outfile}')

def main():
    R1_files = glob.glob('fastq/*_R1_*')
    outfolder = 'fastq_merged/'
    
    # Create output folder if it doesn't exist
    os.makedirs(outfolder, exist_ok=True)

    # Define the number of cores to use for multiprocessing
    num_cores = 24  # Change this to the desired number of cores

    # Use multiprocessing Pool for concurrent processing with the specified number of cores
    with Pool(processes=num_cores) as pool:
        pool.map(merge_files, R1_files)

if __name__ == "__main__":
    main()