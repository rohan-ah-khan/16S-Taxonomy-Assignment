import os
import gzip
import shutil
import subprocess


raw_dir = os.path.join(os.getcwd(), "raw")

if not os.path.exists(raw_dir):
    print(f'"raw" directory does not exist. Make subfolder which contains all raw fast.qz files')
    exit(1)
else:
    file_list = os.listdir(raw_dir)
    gz_file_list = [k for k in file_list if 'fastq.gz' in k]
    
    sample_names = {}
    for file in gz_file_list:
        if "_L001_R1_001.fastq.gz" in file:
            sample_name = file.replace("_L001_R1_001.fastq.gz", "")
            sample_names[sample_name] = sample_names.get(sample_name, 0) | 1  
        elif "_L001_R2_001.fastq.gz" in file:
            sample_name = file.replace("_L001_R2_001.fastq.gz", "")
            sample_names[sample_name] = sample_names.get(sample_name, 0) | 2  

    unmatched_samples = [name for name, value in sample_names.items() if value != 3]
    if unmatched_samples:
        print(f"The following samples are missing either R1 or R2 files: {unmatched_samples}")
        exit(1)  
    else:
        print(f"All samples have matching R1 and R2 files.")    

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
    exit(1)

r1_file_list = [k for k in fastq_file_list if '_L001_R1_001.fastq' in k]
sample_list = []
for n in r1_file_list:
    sample_list.append(n.removesuffix('_L001_R1_001.fastq'))    

with open("calist.txt", 'w') as file:
    for samples in sample_list:
        s = "".join(map(str, samples))
        file.write(s+'\n')

cutadapt_path = "cutadapt.exe"  
output_dir = os.path.join(os.getcwd(), "trimmed")
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

primer_fwd = "CCTACGGGAGGCAGCAG"
primer_rev = "GACTACHVGGGTATCTAATCC"

for sample in sample_list:
    r1_input = os.path.join(extract_dir, f"{sample}_L001_R1_001.fastq")
    r2_input = os.path.join(extract_dir, f"{sample}_L001_R2_001.fastq")
    r1_output = os.path.join(output_dir, f"{sample}_L001_R1_001_trimmed.fastq")
    r2_output = os.path.join(output_dir, f"{sample}_L001_R2_001_trimmed.fastq")
    
    cutadapt_cmd = [
        cutadapt_path,
        "-g", primer_fwd,
        "-G", primer_rev,
        "-o", r1_output,
        "-p", r2_output,
        r1_input,
        r2_input
    ]
    
    try:
        print(f"Running Cutadapt for sample: {sample}")
        subprocess.run(cutadapt_cmd, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error occurred while running Cutadapt for sample {sample}: {e}")
        exit(1)

print("Cutadapt processing complete.")



