# README: DADA2 Workflow for ASV Inference

## Introduction

This document outlines a workflow using the DADA2 package to process
trimmed FASTQ files. The pipeline performs quality filtering and
trimming, error learning, and taxonomy assignment using the SILVA
database. This workflow is adapted from
`https://benjjneb.github.io/dada2/tutorial.html` The workflow assumes
the trimmed files are in the `trimmed` directory created by the Python
script.

## Prerequisites

Before running this workflow, ensure the required libraries are
installed. The `dada2` and `seqinr` packages are necessary. If they are
not installed, use `install.packages("seqinr")` beforehand. For `dada2`
you can use the following to install
``` r
if (!requireNamespace("BiocManager", quietly = TRUE)) 
install.packages("BiocManager") 
BiocManager::install("dada2", version = "3.20")

# Version `3.20` might change based on your version of R. You may also need to update your packages using
# `BiocManager::install(version = '3.20')` before you can use the previous code

```


``` r
# Completely clear the workspace
rm(list=ls(all=T))

# Load required libraries
library(dada2)
library(seqinr)
```

## Working Directory Setup and Data Import

We will set the working directory to match the one used in the Python
script and specify the location of the trimmed FASTQ files. If following
the Python script, this directory should contain the `trimmed` folder
the `raw` folder and the Python script. Change working directory to
match yours.

``` r
#setwd("/MacOS")
path <- "trimmed"
list.files(path)
```

## Step 1: Examine Quality Profile

Make sure to examine the quality to the forward and reverse reads to
decide the values for `filterAndTrim()` in Step 2. Look at Dada2
documentation for information of how to interpret plots and decide
values `https://benjjneb.github.io/dada2/tutorial.html`

``` r
fnFs <- list.files(path, pattern = "_R1_001_trimmed.fastq", full.names = TRUE)
fnRs <- list.files(path, pattern = "_R2_001_trimmed.fastq", full.names = TRUE)

sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

#Examine quality profile of the forward reads of the first 20 samples
plotQualityProfile(fnFs[1:20])

#Examine quality profile of the reverse reads
plotQualityProfile(fnRs[1:20])
```

## Step 2: Filtering and Trimming

Quality filter and trim the sequences using Dada2’s `filterAndTrim()`
function. The filtered files will be save in a new directory name
`filtered`. Look at Dada2 documentation for information of how to decide
values for `filterAndTrim()` the following work well for 341F 785R
primer sequences that have already gone through `cutadapt`. The qc3.csv
file will give more details on the filtering and trimming process.

``` r
filtpath.f <- file.path(path, "filtered/fwd")
filtpath.r <- file.path(path, "filtered/rev")
filtFs <- file.path(filtpath.f, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filtpath.r, paste0(sample.names, "_R_filt.fastq.gz"))

if(length(fnFs) != length(fnRs)) stop("Forward and reverse files do not match.")

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(250,220), maxN=0, maxEE=c(2,6), truncQ=2, rm.phix=T, compress=T, multithread=F)

write.csv(out, file= "qc3.csv")

names(filtFs) <- sample.names
names(filtRs) <- sample.names
```

## Step 3: Learn Error Rates

Estimate error rates from the filtered data. These error models are used
in the next step for ASV inference. Visualize the estimated error rates
and how much do they diverge, look at Dada2 documentation for
information on how to interpret the plot

``` r
set.seed(100)
errF <- learnErrors(filtFs, nbases=1e6, multithread=TRUE, randomize = TRUE)
errR <- learnErrors(filtRs, nbases=1e6, multithread=TRUE, randomize = TRUE)
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)
```

## Step 4: Sample Inference and Pair-End Merging

Infer ASVs for each sample individually using DADA2. Merge forward and
reverse reads for each sample. This step processes samples one by one to
minimize memory usage for large datasets.

``` r
mergers <- vector("list", length(sample.names))
names(mergers) <- sample.names

head(sample.names)

for(sam in sample.names) {
  cat("Processing:", sam, "\n")
  derepF <- derepFastq(filtFs[[sam]])
  ddF <- dada(derepF, err=errF, multithread=TRUE)
  derepR <- derepFastq(filtRs[[sam]])
  ddR <- dada(derepR, err=errR, multithread=TRUE)
  merger <- mergePairs(ddF, derepF, ddR, derepR)
  mergers[[sam]] <- merger
}
```

## Step 5: Inspect Results

Inspect the results of ASV inference and merging to ensure the pipeline
has worked correctly. Only the first sample is inspected here, but you
should inspect all samples before removing the intermediates.

``` r
# Inspect the dada-class objects for the first sample
ddF[[1]]

# Inspect the mergers data.frame for the first sample
head(mergers[[1]])

# Remove intermediate files to save memory
rm(derepF); rm(derepR)
```

## Step 6: Sequence Table Construction

``` r
seqtab <- makeSequenceTable(mergers)

# Number of sequence variants
dim(seqtab)

# Size distribution of sequence variants
table(nchar(getSequences(seqtab)))

# Save sequence table
saveRDS(seqtab, "seqtab.rds")
```

## Step 7: Chimera Removal

Search for and remove chimeric sequences from the sequence table.

``` r
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)

# Check the impact of chimera removal
dim(seqtab.nochim)

sum(seqtab.nochim)/sum(seqtab)

# Size distribution after chimera removal
table(nchar(getSequences(seqtab.nochim)))
```

## Step 8: Sequence Filtering

Filter sequence variants by length, retaining only those of biologically
relevant sizes. Here 400 to 431 bp sequences are being kept, but you
will want to change based on the size distribution you get after chimera
removal. Save the ASV table

``` r
# Filter sequence variants between 400 and 431 bp
seqtab2.nochim <- seqtab.nochim[, nchar(colnames(seqtab.nochim)) %in% seq(400, 431)]
dim(seqtab2.nochim)
table(nchar(getSequences(seqtab2.nochim)))

# Save the sequence variant table
seqnum <- paste0("ASV", seq(ncol(seqtab2.nochim))) 
uniqueSeqs <- as.list(colnames(seqtab2.nochim))   
seqtab2.nochim.transposed <- t(seqtab2.nochim)  
rownames(seqtab2.nochim.transposed) <- as.character(seqnum) 
write.csv(seqtab2.nochim.transposed, file="ASV_abundance.csv") 
```

## Step 9: Write ASV Sequences to FASTA

Export ASV sequences to a FASTA file for downstream analysis.

``` r
write.fasta(uniqueSeqs, seqnum, "Asvs.fasta")
```

## Step 10: Taxonomy Assignment

Assign taxonomy to the sequence variants using the SILVA database.
Replace missing classifications with “Unclassified”. This chunk requires
that the silva database is in the working directory. You can install
from `https://zenodo.org/records/4587955` original database can be found
here
`https://www.arb-silva.de/no_cache/download/archive/current/Exports/`

``` r
# Assign taxonomy
tax <- assignTaxonomy(seqtab2.nochim, "silva_nr99_v138.1_train_set.fa.gz", multithread=TRUE, minBoot=60)
head(tax)

# Create a copy of taxonomy and update missing classifications
tax1 <- as.matrix(tax)
rownames(tax1) <- as.character(seqnum)
y <- which(is.na(tax1) == TRUE) 
tax1[y] <- "Unclassified"       
write.csv(tax1, file="ASV_tax.csv") 
```

## Step 11: Save Workspace

Save the entire R environment to an .Rdata file for reproducibility.

``` r
save.image(file="dada2.Rdata")
```
