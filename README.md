# README: Unzip and Trim FASTQ Files with Cutadapt

## Overview
This repository contains `unzip_trim.py` which is a Python script to process paired-end sequencing data from Illumina platforms. The scripts are designed to:

1. Extract `.fastq.gz` files that are in a `raw` directory.
2. Trim primers from paired-end reads using **Cutadapt**.

The repository is organized into:
- **MacOS Directory:** Contains `unzip_trim.py` script for macOS.
- **Windows Directory:** Contains `unzip_trim.py` script for Windows.

## Requirements
- **Python**
- **Cutadapt** ([Documentation](https://cutadapt.readthedocs.io/en/stable/index.html))

For macOS, you need to install Cutadapt beforehand. Use the following command:
```bash
sudo apt install cutadapt
```

### Install Python on macOS
If Python is not installed, you can install it using **Homebrew**:
```bash
brew install python
```

## Setting Up Your Working Directory
To use the script, you need a working directory (WD) with the following structure:

```
<working_directory>/
├── raw/
│   └── <raw fastq.gz files from Illumina output>
├── unzip_trim.py
```

### Input Files
The `raw` directory must contain paired-end `.fastq.gz` files in the standard Illumina naming format:
```
<samplename>_S<samplenumber>_L001_R<1|2>_001.fastq.gz
```

For example:
```
sample1_S1_L001_R1_001.fastq.gz
sample1_S1_L001_R2_001.fastq.gz
```

### Output
The script will:
1. Extract `.fastq.gz` files into a directory named `extracted/`.
2. Trim primers and save trimmed files in a directory named `trimmed/`.

## Primer Sequences
Currently 16s rRNA 341F 785R primers are being used. To use the script with a specific primer set, update the following lines in `unzip_trim.py` with the desired primer sequences:
```python
primer_fwd = "CCTACGGGAGGCAGCAG"
primer_rev = "GACTACHVGGGTATCTAATCC"
```

## Running the Script

### macOS
1. Install Python and Cutadapt as described above.
2. Navigate to the directory containing the `unzip_trim.py` script.
3. Run the script:
   ```bash
   python3 unzip_trim.py
   ```

### Windows
1. Ensure Python and Cutadapt are installed.
2. Navigate to the directory containing the `unzip_trim.py` script.
3. Run the script:
   ```bash
   python unzip_trim.py
   ```

## Demo Data
A `raw` directory with demo `.fastq.gz` files from BioProject [PRJNA607006](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA607006) is provided for testing. These demo files are not associated with this project.

## Notes
- The script assumes the standard Illumina naming convention for input files.
- Ensure `cutadapt` is in your system's PATH.

## Acknowledgements
Demo data is sourced from BioProject [PRJNA607006](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA607006). The author is not associated with this BioProject.
