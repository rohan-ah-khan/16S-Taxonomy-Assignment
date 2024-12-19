import os
import gzip
import shutil
import subprocess


#get the raw dir path 
raw_dir = os.path.join(os.getcwd(), "raw")

if not os.path.exists(raw_dir):
    print(f'"raw" directory does not exist. Make subfolder which contains all raw fast.qz files')
else:
    file_list = os.listdir(raw_dir)
    gz_file_list = [k for k in file_list if 'fastq.gz' in k]
    if len(gz_file_list) % 2 != 0:
        print(f"The number of R1 and R2 files is not equal, please ensure there is a R1 and R2 file for each sample")

extract_dir = os.path.join(os.getcwd(), "extracted")
if not os.path.exists(extract_dir):
    os.makedirs(extract_dir)

for n in gz_file_list:
    raw_file = os.path.join(raw_dir, n)
    extracted_file = os.path.join(extract_dir, n[:-3])
    with gzip.open(raw_file, "rb") as gz:
        with open(extracted_file, "wb") as extracted:
            shutil.copyfileobj(gz, extracted)

extract_file_list = os.listdir(extract_dir)
fastq_file_list = [k for k in extract_file_list if 'fastq' in k]
if len(fastq_file_list) % 2 != 0:
    print(f"The number of R1 and R2 files is not equal, please ensure there is a R1 and R2 file for each sample")

r1_file_list = [k for k in fastq_file_list if '_L001_R1_001.fastq' in k]
sample_list = []
for n in r1_file_list:
    sample_list.append(n.removesuffix('_L001_R1_001.fastq'))    

with open("calist.txt", 'w') as file:
    for samples in sample_list:
        s = "".join(map(str, samples))
        file.write(s+'\n')

cutadapt_path = "cutadapt"  # Use the command directly, as it is in PATH
output_dir = os.path.join(os.getcwd(), "trimmed")
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# Primer sequences
primer_fwd = "CCTACGGGAGGCAGCAG"
primer_rev = "GACTACHVGGGTATCTAATCC"

for sample in sample_list:
    # Define input and output files
    r1_input = os.path.join(extract_dir, f"{sample}_L001_R1_001.fastq")
    r2_input = os.path.join(extract_dir, f"{sample}_L001_R2_001.fastq")
    r1_output = os.path.join(output_dir, f"{sample}_L001_R1_001_trimmed.fastq")
    r2_output = os.path.join(output_dir, f"{sample}_L001_R2_001_trimmed.fastq")
    
    # Construct the Cutadapt command
    cutadapt_cmd = [
        cutadapt_path,
        "-g", primer_fwd,
        "-G", primer_rev,
        "-o", r1_output,
        "-p", r2_output,
        r1_input,
        r2_input
    ]
    
    # Run the command
    try:
        print(f"Running Cutadapt for sample: {sample}")
        subprocess.run(cutadapt_cmd, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error occurred while running Cutadapt for sample {sample}: {e}")

print("Cutadapt processing complete.")



