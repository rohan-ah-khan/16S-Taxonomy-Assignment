DADA2 Workflow for ASV Inference
================

## Introduction

This document outlines a workflow using the DADA2 package to process
trimmed FASTQ files. The pipeline performs quality filtering and
trimming, error learning, and taxonomy assignment using the SILVA
database. This workflow is adapted from
`https://benjjneb.github.io/dada2/tutorial.html` The workflow assumes
the trimmed files are in the `trimmed` directory created by the Python
script.

------------------------------------------------------------------------

## Prerequisites

Before running this workflow, ensure the required libraries are
installed. The `dada2` and `seqinr` packages are necessary. If they are
not installed, use `install.packages("seqinr")` beforehand. For `dada2`
you can use the following to install
`if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager") BiocManager::install("dada2", version = "3.20")`
Version `3.20` might change based on your version of R. You may also
need to update your packages using
`BiocManager::install(version = '3.20')` before you can use the previous
code

``` r
# Completely clear the workspace
rm(list=ls(all=T))

# Load required libraries
library(dada2)
```

    ## Loading required package: Rcpp

``` r
library(seqinr)
```

------------------------------------------------------------------------

## Working Directory Setup and Data Import

We will set the working directory to match the one used in the Python
script and specify the location of the trimmed FASTQ files. If following
the Python script, this directory should contain the `trimmed` folder
the `raw` folder and the Python script. Change working directory to
match yours.

``` r
#setwd("/Desktop/Github/16S-Taxonomy-Assignment/MacOS")
path <- "trimmed"
list.files(path)
```

    ## [1] "filtered"                                 
    ## [2] "SAMN14119777_S5_L001_R1_001_trimmed.fastq"
    ## [3] "SAMN14119777_S5_L001_R2_001_trimmed.fastq"
    ## [4] "SAMN14119778_S4_L001_R1_001_trimmed.fastq"
    ## [5] "SAMN14119778_S4_L001_R2_001_trimmed.fastq"

------------------------------------------------------------------------

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
```

![](Dada2Workflow_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

``` r
#Examine quality profile of the reverse reads
plotQualityProfile(fnRs[1:20])
```

![](Dada2Workflow_files/figure-gfm/unnamed-chunk-3-2.png)<!-- -->

------------------------------------------------------------------------

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

------------------------------------------------------------------------

## Step 3: Learn Error Rates

Estimate error rates from the filtered data. These error models are used
in the next step for ASV inference. Visualize the estimated error rates
and how much do they diverge, look at Dada2 documentation for
information on how to interpret the plot

``` r
set.seed(100)
errF <- learnErrors(filtFs, nbases=1e6, multithread=TRUE, randomize = TRUE)
```

    ## 5295250 total bases in 21181 reads from 1 samples will be used for learning the error rates.

``` r
errR <- learnErrors(filtRs, nbases=1e6, multithread=TRUE, randomize = TRUE)
```

    ## 4659820 total bases in 21181 reads from 1 samples will be used for learning the error rates.

``` r
plotErrors(errF, nominalQ=TRUE)
```

    ## Warning in scale_y_log10(): log-10 transformation introduced infinite values.
    ## log-10 transformation introduced infinite values.

![](Dada2Workflow_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

``` r
plotErrors(errR, nominalQ=TRUE)
```

    ## Warning in scale_y_log10(): log-10 transformation introduced infinite values.

## ![](Dada2Workflow_files/figure-gfm/unnamed-chunk-5-2.png)<!-- -->

## Step 4: Sample Inference and Pair-End Merging

Infer ASVs for each sample individually using DADA2. Merge forward and
reverse reads for each sample. This step processes samples one by one to
minimize memory usage for large datasets.

``` r
mergers <- vector("list", length(sample.names))
names(mergers) <- sample.names

head(sample.names)
```

    ## [1] "SAMN14119777" "SAMN14119778"

``` r
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

    ## Processing: SAMN14119777 
    ## Sample 1 - 47505 reads in 17924 unique sequences.
    ## Sample 1 - 47505 reads in 32382 unique sequences.
    ## Processing: SAMN14119778 
    ## Sample 1 - 21181 reads in 8303 unique sequences.
    ## Sample 1 - 21181 reads in 14621 unique sequences.

------------------------------------------------------------------------

## Step 5: Inspect Results

Inspect the results of ASV inference and merging to ensure the pipeline
has worked correctly. Only the first sample is inspected here, but you
should inspect all samples before removing the intermediates.

``` r
# Inspect the dada-class objects for the first sample
ddF[[1]]
```

    ## TGAGGAATATTGGTCAATGGGCGCGAGCCTGAACCAGCCAAGTCGCGTGAGGGAAGACGGTCCTATGGATTGTAAACCTCTTTTGTCGGGGAGCAAAAGACGCCACGCGTGGCGTTCCGAGAGTACCCGAAGAAAAAGCATCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGCAGGCGGAAGATCAAGTCAGCGGTAAAATTG 
    ##                                                                                                                                                                                                                                                        721 
    ## CCTACGGGGGGCTGCAGTGAGGAATATTGGTCAATGGGCGCGAGCCTGAACCAGCCAAGTCGCGTGAGGGAAGACGGTCCTATGGATTGTAAACCTCTTTTGTCGGGGAGCAAAAGACGCCACGCGTGGCGTTCCGAGAGTACCCGAAGAAAAAGCATCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGCAGGCGGAAGATCA 
    ##                                                                                                                                                                                                                                                        649 
    ## TCGAGAATCATTCACAATGGGGGAAACCCTGATGGTGCGACGCCGCGTGGGGGAATGAAGGTCTTCGGATTGTAAACCCCTGTCATGTGGGAGCAAATTAAAAAGATAGTACCACAAGAGGAAGAGACGGCTAACTCTGTGCCAGCAGCCGCGGTAATACAGAGGTCTCAAGCGTTGTTCGGAATCACTGGGCGTAAAGCGTGCGTAGGCTGTTTCGTAAGTCGTGTGTGAAAGGCGCGGGCTCAACCCG 
    ##                                                                                                                                                                                                                                                        333 
    ## CCTACGGGCGGCTGCAGTGAGGAATATTGGTCAATGGGCGCGAGCCTGAACCAGCCAAGTCGCGTGAGGGAAGACGGTCCTATGGATTGTAAACCTCTTTTGTCGGGGAGCAAAAGACGCCACGCGTGGCGTTCCGAGAGTACCCGAAGAAAAAGCATCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGCAGGCGGAAGATCA 
    ##                                                                                                                                                                                                                                                        364 
    ## TGAGGAATATTGGTCAATGGCCGGAAGGCTGAACCAGCCAAGTCGCGTGAGGGAATAAGGCCCTACGGGTCGTAAACCTCTTTTGCCAGGGAGCAATGCCGTCCACGTGTGGACGGAAGGAGAGTACCTGGAGAAAAAGCATCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGTAGGCGGCCTGCCAAGTCAGCGGTAAAATT 
    ##                                                                                                                                                                                                                                                        348 
    ## CCTACGGGGGGCTGCAGTCGAGAATCATTCACAATGGGGGAAACCCTGATGGTGCGACGCCGCGTGGGGGAATGAAGGTCTTCGGATTGTAAACCCCTGTCATGTGGGAGCAAATTAAAAAGATAGTACCACAAGAGGAAGAGACGGCTAACTCTGTGCCAGCAGCCGCGGTAATACAGAGGTCTCAAGCGTTGTTCGGAATCACTGGGCGTAAAGCGTGCGTAGGCTGTTTCGTAAGTCGTGTGTGAAA 
    ##                                                                                                                                                                                                                                                        266 
    ## CCTACGGGGGGCTGCAGTGAGGAATATTGGTCAATGGACGAGAGTCTGAACCAGCCAAGTAGCGTGAAGGATGACTGCCCTATGGGTTGTAAACTTCTTTTATATGGGAATAAAATGTTCCACGTGTGGGATTTTGTATGTACCATATGAATAAGGATCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGATCCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGAGCGTAGGTGGATTGTTA 
    ##                                                                                                                                                                                                                                                        338 
    ## TGAGGAATATTGGTCAATGGACGAGAGTCTGAACCAGCCAAGTAGCGTGAAGGATGACTGCCCTATGGGTTGTAAACTTCTTTTATATGGGAATAAAATGTTCCACGTGTGGGATTTTGTATGTACCATATGAATAAGGATCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGATCCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGAGCGTAGGTGGATTGTTAAGTCAGTTGTGAAAGTT 
    ##                                                                                                                                                                                                                                                        340 
    ## CCTACGGGGGGCTGCAGTGAGGAATATTGGTCAATGGCCGGAAGGCTGAACCAGCCAAGTCGCGTGAGGGAATAAGGCCCTACGGGTCGTAAACCTCTTTTGCCAGGGAGCAATGCCGTCCACGTGTGGACGGAAGGAGAGTACCTGGAGAAAAAGCATCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGTAGGCGGCCTGCC 
    ##                                                                                                                                                                                                                                                        278 
    ## GACTACTGGGGTATCTAATCCTGTTCGATACCCACACTTTCGTGCATGAGCGTCAGTTGAGCGCCGGTATGCTGCCTTCGCAATCGGAGTTCTGCGTGATATCTATGCATTTCACCGCTACACCACGCATTCCGCATACTTCTCGCTCACTCAAGAAAACCAGTTTCAACGGCTCGAAGAGGTTGAGCCTCTCAATTTTACCGCTGACTTGATCTTCCGCCTGCGCACCCTTTAAACCCAATAAATCCGG 
    ##                                                                                                                                                                                                                                                        251 
    ## CCTACGGGCGGCTGCAGTCGAGAATCATTCACAATGGGGGAAACCCTGATGGTGCGACGCCGCGTGGGGGAATGAAGGTCTTCGGATTGTAAACCCCTGTCATGTGGGAGCAAATTAAAAAGATAGTACCACAAGAGGAAGAGACGGCTAACTCTGTGCCAGCAGCCGCGGTAATACAGAGGTCTCAAGCGTTGTTCGGAATCACTGGGCGTAAAGCGTGCGTAGGCTGTTTCGTAAGTCGTGTGTGAAA 
    ##                                                                                                                                                                                                                                                        192 
    ## GACTACTGGGGTATCTAATCCTGTTTGATACCCACACTTTCGAGCCTCAGTGTCAGTTGCAGTCCAGTGAGCTGCCTTCGCAATCGGAGTTCTTCGTGATATCTAAGCATTTCACCGCTACACCACGAATTCCGCCCACCTCTACTGTACTCAAGACTGCCAGTTTCAACTGCAATTTTACGGTTGAGCCGCAAACTTTCACAACTGACTTAACAATCCACCTACGCTCCCTTTAAACCCAATAAATCCG 
    ##                                                                                                                                                                                                                                                        190 
    ## CCTACGGGGGGCTGCAGTGAGGAATATTGGTCAATGGACGCAAGTCTGAACCAGCCATGCCGCGTGCAGGAAGACGGCTCTATGAGTTGTAAACTGCTTTTGTACGAGGGTAAACGCAGATACGTGTATCTGTCTGAAAGTATCGTACGAATAAGGATCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGATTCAAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGTAGGCGGTCGGATA 
    ##                                                                                                                                                                                                                                                        196 
    ## TGAGGAATATTGGTCAATGGGCGCGAGCCTGAACCAGCCAAGTAGCGTGCAGGACGACGGCCCTATGGGTTGTAAACTGCTTTTATACGGGGATAAAGTATGCCACGTGTGGTTTATTGCAGGTACCGTATGAATAAGGACCGGCTAATTCCGTGCCAGCAGCCGCGGTAATACGGAAGGTCCGGGCGTTATCCGGATTTATTGGGTTTAAAGGGAGCGTAGGCCGGGTTTTAAGCGTGCCGTGAAATGT 
    ##                                                                                                                                                                                                                                                        195 
    ## TGAGGAATATTGGTCAATGGACGCAAGTCTGAACCAGCCATGCCGCGTGCAGGAAGACGGCTCTATGAGTTGTAAACTGCTTTTGTACGAGGGTAAACGCAGATACGTGTATCTGTCTGAAAGTATCGTACGAATAAGGATCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGATTCAAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGTAGGCGGTCGGATAAGTTAGAGGTGAAATCC 
    ##                                                                                                                                                                                                                                                        188 
    ## GACTACTGGGGTATCTAATCCCTTTCGCTCCCCTGGCCTTCGTGCCTCAGCGTCAGTTAATGTCCAGGAACCCGCCTTCGCCACGAGTGTTCCTCTCGATATCTACGCATTTCACTGCTACACCGAGAATTCCGGTTCCCCCTCCATTACTCTAGTCTCGCAGTATCATGTGCCGTCCGCGGGTTGAGCCCGCGCCTTTCACACACGACTTACGAAACAGCCTACGCACGCTTTACGCCCAGTGATTCCG 
    ##                                                                                                                                                                                                                                                        152 
    ## CCTACGGGGGGCTGCAGTGAGGAATATTGGTCAATGGGCGCGAGCCTGAACCAGCCAAGTAGCGTGCAGGACGACGGCCCTATGGGTTGTAAACTGCTTTTATACGGGGATAAAGTATGCCACGTGTGGTTTATTGCAGGTACCGTATGAATAAGGACCGGCTAATTCCGTGCCAGCAGCCGCGGTAATACGGAAGGTCCGGGCGTTATCCGGATTTATTGGGTTTAAAGGGAGCGTAGGCCGGGTTTTA 
    ##                                                                                                                                                                                                                                                        157 
    ## CCTACGGGGGGCTGCAGTGAGGAATATTGGTCAATGGCCGGAAGGCTGAACCAGCCAAGTCGCGTGAGGGAATAAGGCCCTACGGGTCGTAAACCTCTTTTGTCAGGGAGCAAGGCCGCCCACGTGTGGACGGAAGGAGAGTACCTGGAGAAAAAGCATCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGTAGGCGGCCTGCC 
    ##                                                                                                                                                                                                                                                        108 
    ## GACTACTGGGGTATCTAATCCTGTTTGATCCCCGCACTTTCGTGCCTCAGCGTCAGTAGGGCGCCGGTATGCTGCCTTCGCAATCGGGGTTCTGCGTGATATCTATGCATTTCACCGCTACACCACGCATTCCGCATACTTCTCGCCCACTCGAGCCCGGCAGTTTCAACGGCTGTACGGGGTTGAGCCCCGCAATTTTACCGCTGACTTGGCAGGCCGCCTACGCACCCTTTAAACCCAATAAATCCGG 
    ##                                                                                                                                                                                                                                                        135 
    ## TGAGGAATATTGGTCAATGGCCGGAAGGCTGAACCAGCCAAGTCGCGTGAGGGAATAAGGCCCTACGGGTCGTAAACCTCTTTTGTCAGGGAGCAAGGCCGCCCACGTGTGGACGGAAGGAGAGTACCTGGAGAAAAAGCATCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGTAGGCGGCCTGCCAAGTCAGCGGTAAAATT 
    ##                                                                                                                                                                                                                                                        101 
    ## GACTACTGGGGTATCTAATCCTGTTTGCTCCCCACGCTTTCGTGCCTCAACGTCAGAAGTAGCTTGGTAAGCTGCCTTCGCAATCGGTGTTCTGTATGATCTCTAAGCATTTCACCGCTACACCATACATTCCGCCTACCGCAACTACTCTCTAGCCGAACAGTATCAGAGGCAATTCCGAAGTTGAGCCTCGGGATTTCACCTCTAACTTATCCGACCGCCTACGCACCCTTTAAACCCAATAAATCCG 
    ##                                                                                                                                                                                                                                                        117 
    ## GACTACTGGGGTATCTAATCCTGTTCGATACCCGCGCCTTCGAGCTTCAGCGTCAGTAGCGCTGCCGTATGCTGCCTTCGCAATCGGGGTTCTTCGTGATATCTAAGCATTTCACCGCTACACCACGAATTCCGCATACGTTCCGCGCACTCAAGGACTCCAGTTCGCGCCGCAGTGTCAAGGTTGAGCCCCGACATTTCACGGCACGCTTAAAACCCGGCCTACGCTCCCTTTAAACCCAATAAATCCG 
    ##                                                                                                                                                                                                                                                        123 
    ## TGAGGAATATTGGTCAATGGGAGAGATCCTGAACCAGCCAAGCCGCGTGAGGGAAGACGGCCCTATGGGTTGTAAACCTCTTTTGTCGGAGAACAAAACCCGGGACGTGTCCCGGACTGCGTGTATCCGAAGAAAAAGCATCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGTAGGCGGTCCGTTAAGTCAGCGGTAAAATTG 
    ##                                                                                                                                                                                                                                                         79 
    ## TGAGGAATATTGGTCAATGGGCGAGAGCCTGAACCAGCCAAGTCGCGTGAGGGAGTACTGCCCTATGGGTTGTAAACCTCTTTTGTCGGGGAGCAAAAGCCGGACGTGTCCGTCTGTGAGAGTACCCGAAGAAAAAGCATCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGCAGGCGGATTTTTAAGTCAGCGGTCAAATCGT 
    ##                                                                                                                                                                                                                                                         63 
    ## TGGGGAATATTGCACAATGGGCGAAAGCCTGATGCAGCAACGCCGCGTGAAGGAAGACGGTTTTCGGATTGTAAACTTCTGTTCTTAGTGAAGAATAATGACGGTAACTAAGGAGCAAGCCACGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGTGGCAAGCGTTGTCCGGAATTACTGGGTGTAAAGGGAGCGTAGGCGGGAAGCCAAGTCAGCTGTGAAAACTACGGGCTTAACCTGTAGAC 
    ##                                                                                                                                                                                                                                                         56 
    ## TCGGGAATATTGCGCAATGGAGGAAACTCTGACGCAGTGACGCCGCGTATAGGAAGAAGGTTTTCGGATTGTAAACTATTGTCGTTAGGGAAGATAAAAGACAGTACCTAAGGAGGAAGCTCCGGCTAACTATGTGCCAGCAGCCGCGGTAATACATAGGGAGCAAGCGTTATCCGGAATTATTGGGTGTAAAGGGTGCGTAGACGGAGGAACAAGTTAGTTGTGAAATCCCTCGGCTTAACTGAGGAAC 
    ##                                                                                                                                                                                                                                                         63 
    ## CCTACGGGGGGCTGCAGTCGGGAATATTGCGCAATGGAGGAAACTCTGACGCAGTGACGCCGCGTATAGGAAGAAGGTTTTCGGATTGTAAACTATTGTCGTTAGGGAAGATAAAAGACAGTACCTAAGGAGGAAGCTCCGGCTAACTATGTGCCAGCAGCCGCGGTAATACATAGGGAGCAAGCGTTATCCGGAATTATTGGGTGTAAAGGGTGCGTAGACGGAGGAACAAGTTAGTTGTGAAATCCCT 
    ##                                                                                                                                                                                                                                                         64 
    ## TGAGGAATATTGGTCAATGGGCGGGAGCCTGAACCAGCCAAGTCGCGTGAGGGACGACGGTCCTATGGATTGTAAACCTCTTTTGTCGGGGAACAAAGGCGCCCACGGGTGGGCGGATGAGTGTACCCGAAGAAAAAGCATCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGTAGGCGGACGGTCAAGTCAGCGGTAAAAATG 
    ##                                                                                                                                                                                                                                                         50 
    ## TGAGGAATATTGGTCAATGGCCGAGAGGCTGAACCAGCCAAGTCGCGTGAGGGAAGACGGCCCTACGGGTTGTAAACCTCTTTTGTCGGGGAGCAAAGAGCGGGATGCGTCCCGCCTCGAGAGTACCCGAAGAAAAAGCATCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGCAGGCGGAAGATCAAGTCAGCGGTAAAATCG 
    ##                                                                                                                                                                                                                                                         40 
    ## TGGGGAATATTGCACAATGGGGGAAACCCTGATGCAGCGACGCCGCGTGAAGGAAGAAGTATTTCGGTATGTAAACTTCTATCAGCAGGGAAGAAAATGACGGTACCTGAGTAAGAAGCCCCGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGGGGCAAGCGTTATCCGGATTTACTGGGTGTAAAGGGAGCGTAGACGGTAATGTAAGTCTGATGTGAAAGCCCGGGGCTCAACCCCGGGACT 
    ##                                                                                                                                                                                                                                                         53 
    ## TGAGGAATATTGGTCAATGGGCGGGAGCCTGAACCAGCCAAGTCGCGTGAGGGAAGACGGTCCTATGGATTGTAAACCTCTTTAGGCGGGGAGCAATGCCGGGCACGCGTGCCCGGAGGGAGAGTACCCGCAGAATAAGCATCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGTAGGCGGATTTTTAAGTCAGCGGTAAAATG 
    ##                                                                                                                                                                                                                                                         49 
    ## CCTACGGGGGGCTGCAGTGAGGAATATTGGTCAATGGGCGGGAGCCTGAACCAGCCAAGTCGCGTGAGGGAAGACGGTCCTATGGATTGTAAACCTCTTTAGGCGGGGAGCAATGCCGGGCACGCGTGCCCGGAGGGAGAGTACCCGCAGAATAAGCATCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGTAGGCGGATTTTT 
    ##                                                                                                                                                                                                                                                         77 
    ## TGGGGAATATTGGGCAATGGAGGAAACTCTGACCCAGCAACGCCGCGTGAATGATGAAGGTCTTCGGATTGTAAAGTTCTTTTCTAAGGGAAGAAGAAAGTGACGGTACCTTAGGAATAAGCCTCGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGAGGCAAGCGTTATCCGGAATGACTGGGCGTAAAGGGTGCGTAGGTGGTTCAGCAAGTTAGTAGCGTAACTCCGGGGCTCAACTCCGGC 
    ##                                                                                                                                                                                                                                                         49 
    ## CCTACGGGGGGCTGCAGTGGGGAATATTGCACAATGGGCGAAAGCCTGATGCAGCAACGCCGCGTGAAGGAAGACGGTTTTCGGATTGTAAACTTCTGTTCTTAGTGAAGAATAATGACGGTAACTAAGGAGCAAGCCACGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGTGGCAAGCGTTGTCCGGAATTACTGGGTGTAAAGGGAGCGTAGGCGGGAAGCCAAGTCAGCTGTGAAAACTAC 
    ##                                                                                                                                                                                                                                                         52 
    ## CCTACGGGGGGCTGCAGTGAGGAATATTGGTCAATGGGAGAGATCCTGAACCAGCCAAGCCGCGTGAGGGAAGACGGCCCTATGGGTTGTAAACCTCTTTTGTCGGAGAACAAAACCCGGGACGTGTCCCGGACTGCGTGTATCCGAAGAAAAAGCATCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGTAGGCGGTCCGTTA 
    ##                                                                                                                                                                                                                                                        139 
    ## CCTACGGGGGGCTGCAGTGAGGAATATTGGTCAATGGCCGAGAGGCTGAACCAGCCAAGTCGCGTGAGGGAAGACGGCCCTACGGGTTGTAAACCTCTTTTGTCGGGGAGCAAAGAGCGGGATGCGTCCCGCCTCGAGAGTACCCGAAGAAAAAGCATCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGCAGGCGGAAGATCA 
    ##                                                                                                                                                                                                                                                         91 
    ## GACTACTGGGGTATCTAATCCTGTTCGATACCCACGCTTTCGTGCATGAGCGTCAGTTGGGCGCCGGTACGCTGCCTTCGCAATCGGAGTTCTGCGTGATATCTATGCATTTCACCGCTACACCACGCATTCCGCGTACTTCTCACCCACTCAAGAAAACCAGTTTCAACGGCTGAAAGAGGTTGAGCCTCTCGATTTTACCGCTGACTTGATCTTCCGCCTGCGCACCCTTTAAACCCAATAAATCCGG 
    ##                                                                                                                                                                                                                                                        201 
    ## CCTACGGGGGGCTGCAGTGAGGAATATTGGTCAATGGGCGGGAGCCTGAACCAGCCAAGTCGCGTGAGGGACGACGGTCCTATGGATTGTAAACCTCTTTTGTCGGGGAACAAAGGCGCCCACGGGTGGGCGGATGAGTGTACCCGAAGAAAAAGCATCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGTAGGCGGACGGTCA 
    ##                                                                                                                                                                                                                                                         65 
    ## TGAGGAATATTGGTCAATGGCCGGAAGGCTGAACCAGCCAAGCCGCGTGAGGGAGGAAGGCGCAGAGCGTCGCAGACCTCTTTTGCCGGGGGACAAAAGGCCGGACTCGTCCGGTCCTGAGGGTACCCGGAGAAAAAGCATCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGATGCGAGCGTTATCCGGATTCATTGGGTTTAAAGGGTGCGCAGGCGGGCCTTTAAGCCGGCGGTGAAATCG 
    ##                                                                                                                                                                                                                                                         36 
    ## TGGGGAATATTGCACAATGGGGGAAACCCTGATGCAGCGATGCCGCGTGGAGGAAGAAGGTTTTCGGATTGTAAACTCCTTTTATCAGGGACGATAATGACGGTACCTGAAGAAAAAGCTCCGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGGAGCGAGCGTTGTCCGGAATTACTGGGTGTAAAGGGAGCGTAGGCGGGACGGCAAGTCAGATGTGAAAACTGAGGGCTCAACCTTCAGACT 
    ##                                                                                                                                                                                                                                                         32 
    ## GACTACCGGGGTATCTAATCCTGTTCGATCCCCACGCTTTCGTGCCTCAGCGTCAGTCTGGCGCCGGTACGCTGCCTTCGCAATCGGAGTTCTGCGCGATATCTATGCATTTCACCGCTACACCGCGCATTCCGCGTACTTCTCGCCAACTCAAGTCTGCCAGTTTCAACGGCTCGACGGGGTTGAGCCCCGCAATTTTACCGCTGACTTAACGGACCGCCTACGCACCCTTTAAACCCAATAAATCCGG 
    ##                                                                                                                                                                                                                                                        322 
    ## CCTACGGGGGGCTGCAGTGAGGAATATTGGTCAATGGGCGAGAGCCTGAACCAGCCAAGTCGCGTGAGGGAGTACTGCCCTATGGGTTGTAAACCTCTTTTGTCGGGGAGCAAAAGCCGGACGTGTCCGTCTGTGAGAGTACCCGAAGAAAAAGCATCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGCAGGCGGATTTTTAA 
    ##                                                                                                                                                                                                                                                         70 
    ## CCTACGGGGGGCTGCAGTGAGGAATATTGGTCAATGGGCGGGAGCCTGAACCAGCCAAGTCGCGTGAGGGAAGAAGGTACTGAGTATTGTAAACCTCTTTTGCCAGGGGACAAAGGCGGCCACGGGTGGCCGTAAGAGGGTACCTGGAGAAAAAGCATCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGCAGGCGGCGTGTCA 
    ##                                                                                                                                                                                                                                                         60 
    ## TGGGGAATTTTGGACAATGGGCGCAAGCCTGATCCAGCTATTCCGCGTGTGGGATGACGGCCCTCGGGTTGTAAACCACTTTAGTAGAGAACGAAAAGGAGCTATCGAATAAATAGCTCTGCTGACGGTACTCTAAGAATAAGCACCGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGGTGCGAGCGTTAATCGGAATTACTGGGCGTAAAGGGTGTGCAGGCGGCTGAGTAAGACAGATGTGA 
    ##                                                                                                                                                                                                                                                         44 
    ## GACTACTGGGGTATCTAATCCTGTTCGATACCCGCACTTTCGTGCCTCAGCGTCAGTTGAGCGCCGGTAAGCTGCCTTCGCAATCGGCGTTCTGCGTGATATCTATGCATTTCACCGCTACACCACGCATTCCGCTTACTTCTCACTCACTCAAGGGCGTCAGTTTCAACGGCACGACGGGGTTGAGCCCCGAAATTTTACCGCTGACTTGACACTCCGCCTACGCACCCTTTAAACCCAATAAATCCGG 
    ##                                                                                                                                                                                                                                                        136 
    ## TGAGGGATATTGGTCAATGGGGGAAACCCTGAACCAGCAACGCCGCGTGAGGGAAGACGGTTTTCGGATTGTAAACCTCTGTCCTCTGTGAAGATAATGACGGTAGCAGAGGAGGAAGCTCCGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGGAGCGAGCGTTGTCCGGATTTACTGGGTGTAAAGGGTGCGTAGGCGGGTAAGCAAGTCAGAAGTGAAATCCATGGGCTTAACCCATGAACT 
    ##                                                                                                                                                                                                                                                         29 
    ## CCTACGGGGGGCTGCAGTGGGGAATATTGCACAATGGGGGAAACCCTGATGCAGCGACGCCGCGTGAAGGAAGAAGTATTTCGGTATGTAAACTTCTATCAGCAGGGAAGAAAATGACGGTACCTGAGTAAGAAGCCCCGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGGGGCAAGCGTTATCCGGATTTACTGGGTGTAAAGGGAGCGTAGACGGTAATGTAAGTCTGATGTGAAAGCCCGG 
    ##                                                                                                                                                                                                                                                         60 
    ## CCTACGGGGGGCTGCAGTGGGGAATATTGGGCAATGGAGGAAACTCTGACCCAGCAACGCCGCGTGAATGATGAAGGTCTTCGGATTGTAAAGTTCTTTTCTAAGGGAAGAAGAAAGTGACGGTACCTTAGGAATAAGCCTCGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGAGGCAAGCGTTATCCGGAATGACTGGGCGTAAAGGGTGCGTAGGTGGTTCAGCAAGTTAGTAGCGTAACTC 
    ##                                                                                                                                                                                                                                                         53 
    ## GACTACTCGGGTATCTAATCCTGTTTGCTCCCCACACTTTCGAGCCTCAGCGTCAGTTAAAGCCCAGTTGGCCGCCTTCGCCACCGGTGTTCCTCCGAATATCTACGCATTTCACCGCTACACTCGGAATTCCGCCAACCTCTACTTCACTCAAGAAAGCCAGTTTCAACTGCAGTCTACAGGTTAAGCCCGTAGTTTTCACAGCTGACTTGGCTTCCCGCCTACGCTCCCTTTACACCCAGTAATTCCG 
    ##                                                                                                                                                                                                                                                         54 
    ## TGAGGAATATTGGTCAATGGCCGGAAGGCTGAACCAGCCAAGTCGCGTGAGGGAATAAGGCCCTATGGGTCGTAAACCTCTTTTGTCAGGGAGCAAAGGCGTCCACGTGTGGACGAGAGGAGAGTACCTGAAGAAAAAGCATCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGTAGGCGGAGTGTCAAGTCAGCGGTAAAATT 
    ##                                                                                                                                                                                                                                                         26 
    ## CCTACGGGGGGCTGCAGTGAGGAATATTGGTCAATGGCCGGAAGGCTGAACCAGCCAAGTCGCGTGAGGGAATAAGGCCCTATGGGTCGTAAACCTCTTTTGTCAGGGAGCAAAGGCGTCCACGTGTGGACGAGAGGAGAGTACCTGAAGAAAAAGCATCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGTAGGCGGAGTGTC 
    ##                                                                                                                                                                                                                                                         67 
    ## TGAGGAATATTGGTCAATGGGCGGGAGCCTGAACCAGCCAAGCCGCGTGAAGGAAGAAGGCGCCAAGCGTCGTAAACTTCTTTTGTCGGGGAACAAAGGGCGCCACGTGTGGCGTTGTGAGTGTACCCGAAGAAAAAGCATCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGTAGGCGGTCCGGAAAGTCAGCGGTAAAAGCC 
    ##                                                                                                                                                                                                                                                         26 
    ## GACTACTGGGGTATCTAATCCTGTTTGCTCCCCACGCTTTCGAGCCTCAACGTCAGTTTCCGTCCAGTAAGCCGCCTTCGCCACTGGTGTTCCTCCTAATATCTACGCATTTCACCGCTACACTAGGAATTCCGCTTACCTCTCCGGTACTCTAGTTACATAGTTTCCAATGCAGTCCCGGGGTTGAGCCCCGGGCTTTCACATCAGACTTACATTACCGTCTACGCTCCCTTTACACCCAGTAAATCCG 
    ##                                                                                                                                                                                                                                                         94 
    ## GACTACCGGGGTATCTAATCCTGTTCGATACCCACGCTTTCGTGCCTGAGCGTCAGTTGTGCACCGGTATGCTGCCTTCGCAATCGGGGTTCTGCGTGATATCTATGCATTTCACCGCTACACCACGCATTCCGCATACCTCTCGCACACTCTAGATCCCCAGTTTCAACGGCTGGATGGGGTTGAGCCCCACGATTTGACCGCTGACTTAAAAATCCGCCTGCGCACCCTTTAAACCCAATAAATCCGG 
    ##                                                                                                                                                                                                                                                        131 
    ## CCTACGGGGGGCTGCAGTGAGGGATATTGGTCAATGGGGGAAACCCTGAACCAGCAACGCCGCGTGAGGGAAGACGGTTTTCGGATTGTAAACCTCTGTCCTCTGTGAAGATAATGACGGTAGCAGAGGAGGAAGCTCCGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGGAGCGAGCGTTGTCCGGATTTACTGGGTGTAAAGGGTGCGTAGGCGGGTAAGCAAGTCAGAAGTGAAATCCATG 
    ##                                                                                                                                                                                                                                                         47 
    ## CCTACGGGGGGCTGCAGTGGGGAATATTGCACAATGGGGGAAACCCTGATGCAGCGATGCCGCGTGGAGGAAGAAGGTTTTCGGATTGTAAACTCCTTTTATCAGGGACGATAATGACGGTACCTGAAGAAAAAGCTCCGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGGAGCGAGCGTTGTCCGGAATTACTGGGTGTAAAGGGAGCGTAGGCGGGACGGCAAGTCAGATGTGAAAACTGAG 
    ##                                                                                                                                                                                                                                                         43 
    ## GACTACCGGGGTATCTAATCCTGTTCGATCCCCACGCTTTCGTGCCTCAGCGTCAGTCGGGCGCCGGTATGCTGCCTTCGCAATCGGAGTTCTGCGTGATATCTATGCATTTCACCGCTACACCACGCATTCCGCATACTTCTCGCCCACTCAAGATCCCCAGTTTCAACGGCCGGCCCGGGTTGAGCCCGGACATTTTACCGCTGACTTAAAAATCCGCCTACGCACCCTTTAAACCCAATAAATCCGG 
    ##                                                                                                                                                                                                                                                        120 
    ## TGAGGAATATTGGTCAATGGGCGGGAGCCTGAACCAGCCAAGTCGCGTGAGGGAAGAAGGTACTGAGTATTGTAAACCTCTTTTGCCAGGGGACAAAGGCGGCCACGGGTGGCCGTAAGAGGGTACCTGGAGAAAAAGCATCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGCAGGCGGCGTGTCAAGTCAGCGGTAAAATGT 
    ##                                                                                                                                                                                                                                                         26 
    ## TGGGGGATATTGGACAATGGGGGAAACCCTGATCCAGCGACGCCGCGTGAGTGAAGAAGTATTTCGGTATGTAAAGCTCTGTCAGCAGGGAAGAAAGAAATGACGGTACCTGACCAAGAAGCCCCGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGGGGCAAGCGTTATCCGGAATTACTGGGTGTAAAGGGAGCGTAGACGGTGATGCAAGTCTGGAGTGAAAGGCGGGGGCCCAACCCCCGG 
    ##                                                                                                                                                                                                                                                         23 
    ## CCTACGGGGGGCTGCAGTGGGGGATATTGGACAATGGGGGAAACCCTGATCCAGCGATGCCGCGTGAGGGAAGAAGGTTTTCGGATTGTAAACCTCTGTGGAGGGGGGCGATAATGACGGTACCCCTTTAGGAAGCCACGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGTGGCGAGCGTTGTCCGGAATTACTGGGTGTAAAGGGAGTGTAGGCGGGAAAGCAAGTCAGAAGTGAAAACTATG 
    ##                                                                                                                                                                                                                                                         30 
    ## GACTACTGGGGTATCTAATCCTGTTTGCTCCCCACGCTTTCGTGTATGAGCGTCAGTTGCGTCCCAGGAGGCTGCCTTCGCCATCGGTGTTCTTCCACATATCTACGCATTTCACTGCTACACGTGGAATTCCACCTCCCTCTGACGCACTCTAGCCAAGCAGTCACATATGCGGTTCCCAGGTTAAGCCCGGGGATTTCACATCTGTCTTACTCAGCCGCCTGCACACCCTTTACGCCCAGTAATTCCG 
    ##                                                                                                                                                                                                                                                         92 
    ## TGAGGAATATTGGTCAATGGGCGAGAGCCTGAACCAGCCAAGTCGCGTGAAGGATGAAGGATCTATGGTTTGTAAACTTCTTTTATATGGGAATAAAGTGAGGAACGTGTTCCTTTTTGTATGTACCATATGAATAAGCATCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGTAGGTGGTTAATTAAGTCAGCGGTGAAAGTT 
    ##                                                                                                                                                                                                                                                         20 
    ## GACTACTGGGGTATCTAATCCTGTTTGCTCCCCACACTTTCGTGCCTCAACGTCAGTTATCGTCCAGAAAGTCGCCTTCGCCACTGGTGTTCCTCCTAATATCTACGCATTTCACCGCTACACTAGGAATTCCACTTTCCTCTCCGACACTCAAGGGGAATAGTTTTAGTTGCAGTTCCTCAGTTAAGCCGAGGGATTTCACAACTAACTTGTTCCTCCGTCTACGCACCCTTTACACCCAATAATTCCG 
    ##                                                                                                                                                                                                                                                        147 
    ## CCTACGGGGGGCTGCAGTGAGGAATATTGGTCAATGGGCGGGAGCCTGAACCAGCCAAGTCGCGTGAGGGATGACGGTCCTACGGATTGTAAACCTCTTTTGCCGGGGAGCAAGCGTGCGTTCGTGAACGCGCGTCGAGAGTACCCGGAGAAAAAGCATCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGTAGGCGGGATGCC 
    ##                                                                                                                                                                                                                                                         50 
    ## CCTACGGGGGGCTGCAGTGAGGAATATTGGTCAATGGGCGCAAGCCTGAACCAGCCAAGTCGCGTGAGGGAAGACGGTCCTATGGATTGTAAACCTCTTTTGTCGGGGAGCAAAGCGCGGTACGAGTACCGCGAAGGAGAGTACCCGAAGAAAAAGCATCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGTAGGCGGATTTGT 
    ##                                                                                                                                                                                                                                                         53 
    ## CCTACGGGGGGCTGCAGTGGGGAATTTTGGACAATGGGCGCAAGCCTGATCCAGCTATTCCGCGTGTGGGATGACGGCCCTCGGGTTGTAAACCACTTTAGTAGAGAACGAAAAGGAGCTATCGAATAAATAGCTCTGCTGACGGTACTCTAAGAATAAGCACCGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGGTGCGAGCGTTAATCGGAATTACTGGGCGTAAAGGGTGTGCAGGCGGCT 
    ##                                                                                                                                                                                                                                                         41 
    ## CCTACGGGGGGCTGCAGTGGGGAATATTGCACAATGGGGGAAACCCTGATGCAGCAACGCCGCGTGAGTGAAGAAGTATTTCGGTATGTAAAGCTCTATCAGCAGGGAAGAAGAAAGACGGTACCTGACTAAGAAGCCCCGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGGGGCAAGCGTTATCCGGATTTACTGGGTGTAAAGGGAGCGTAGACGGCTGTGCAAGTCTGAAGTGAAAGCCCG 
    ##                                                                                                                                                                                                                                                         70 
    ## TCGGGAATATTGCGCAATGGAGGCAACTCTGACGCAGTGACGCCGCGTATAGGAAGAAGGTTTTCGGATTGTAAACTATTGTCCTCAGGGAAGAAAAAAGACTGTACCTGAGAAGAAAGCTCCGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGGAGCAAGCGTTATCCGGAATTATTGGGTGTAAAGGGTGCGTAGACGGGATAACAAGTTGGTTGTGAAACCCCTCGGCTCAACTGAGGAAC 
    ##                                                                                                                                                                                                                                                         16 
    ## TGAGGAATATTGGTCAATGGGCGCAAGCCTGAACCAGCCAAGTCGCGTGAGGGAAGACGGTCCTATGGATTGTAAACCTCTTTTGTCGGGGAGCAAAGCGCGGTACGAGTACCGCGAAGGAGAGTACCCGAAGAAAAAGCATCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGTAGGCGGATTTGTAAGTCAGCGGTAAAAAT 
    ##                                                                                                                                                                                                                                                         33 
    ## TGAGGAATATTGGTCAATGGCCGGGAGGCTGAACCAGCCAAGTCGCGTGAGGGAAGACGGCCCTACGGGTTGTAAACCTCTTTTGTCAGGGAGCAAGAGGCGGGTCGGGACCTGCTGTGAGAGTACCTGAAGAAAAAGCATCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGCAGGCGGAGGCACAAGTCAGCGGTAAAATCG 
    ##                                                                                                                                                                                                                                                         16 
    ## TGGGGAATATTGCACAATGGGGGAAACCCTGATGCAGCGACGCCGCGTGAGTGAAGAAGTATTTCGGTATGTAAAGCTCTATCAGCAGGGAAGAGAATGACGGTACCTGAGTAAGAAGCCCCGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGGGGCAAGCGTTATCCGGATTTACTGGGTGTAAAGGGAGCGTAGACGGCGAAGAAAGTCTGAAGTGAAAGCCCGCGGCTCAACCGCGGAACT 
    ##                                                                                                                                                                                                                                                         15 
    ## CCTACGGGGGGCTGCAGTGAGGAATATTGGTCAATGGCCGGAAGGCTGAACCAGCCAAGCCGCGTGAGGGAGGAAGGCGCAGAGCGTCGCAGACCTCTTTTGCCGGGGGACAAAAGGCCGGACTCGTCCGGTCCTGAGGGTACCCGGAGAAAAAGCATCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGATGCGAGCGTTATCCGGATTCATTGGGTTTAAAGGGTGCGCAGGCGGGCCTTTA 
    ##                                                                                                                                                                                                                                                         32 
    ## GACTACCCGGGTATCTAATCCTGTTTGCTCCCCACGCTTTCGTGCCTCAGTGTCAGTTACAGTCCAGCAACTCGCCTTCGCCACCGGTGTTCTTCCTAATATCTACGCATTCCACCGCTACACTAGGAATTCCAGTTGCCCCTCCTGCACTCAAGTCCAACAGTTTTAGTAGTAGTGCCGGAGTTGAGCCCCGGAGTTACGCTACTAACTTGCTGAACCACCTACGCACCCTTTACGCCCAGTCATTCCG 
    ##                                                                                                                                                                                                                                                         64 
    ## TGGGGAATATTGCGCAATGGGGGCAACCCTGACGCAGCAACGCCGCGTGATTGAAGAAGGCCTTCGGGTTGTAAAGATCTTTAATCGGGGACGAATTTTGACGGTACCCGAAGAATAAGCTCCGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGGAGCAAGCGTTATCCGGATTTACTGGGTGTAAAGGGCGTGTAGGCGGGCATGCAAGCCAGAAGTGAAATCTGGGGGCTTAACCCCCAAAC 
    ##                                                                                                                                                                                                                                                         11 
    ## TGAGGAATATTGGTCAATGGACGAAAGTCTGAACCAGCCATGCCGCGTGCAGGATGACGGCTCTATGAGTTGTAAACTGCTTTTGTATTAGGGTAAACACCCCTACGTGTAGGGGCTTGAAAGTATAATACGAATAAGGATCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGATCCAAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGTAGGCGGTTTGATAAGTTAGAGGTGAAATAC 
    ##                                                                                                                                                                                                                                                         16 
    ## GACTACTGGGGTATCTAATCCTGTTTGCTACCCACGCTTTCGTGCCTCAACGTCAGATATAGTTTGGTAAGCTGCCTTCGCAATCGGTGTTCTGTATGATCTCTAAGCATTTCACCGCTACACCATACATTCCGCCTACCGCAACTACTCTCTAGTTCAACAGTATTAGAGGCAGTTCCGGTGTTGAGCACCGGTATTTCACCTCTAACTTATCAAACCGCCTACGCACCCTTTAAACCCAATAAATCCG 
    ##                                                                                                                                                                                                                                                         50 
    ## CCTACGGGGGGCTGCAGTGAGGAATATTGGTCAATGGGCGGGAGCCTGAACCAGCCAAGCCGCGTGAAGGAAGAAGGCGCCAAGCGTCGTAAACTTCTTTTGTCGGGGAACAAAGGGCGCCACGTGTGGCGTTGTGAGTGTACCCGAAGAAAAAGCATCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGTAGGCGGTCCGGAA 
    ##                                                                                                                                                                                                                                                         24 
    ## CCTACGGGGGGCTGCAGTGAGGAATATTGGTCAATGGGCGAGAGCCTGAACCAGCCAAGTCGCGTGAAGGATGAAGGATCTATGGTTTGTAAACTTCTTTTATATGGGAATAAAGTGAGGAACGTGTTCCTTTTTGTATGTACCATATGAATAAGCATCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGTAGGTGGTTAATTA 
    ##                                                                                                                                                                                                                                                         25 
    ## CCTACGGGGGGCTGCAGTGAGGAATATTGGTCAATGGACGAAAGTCTGAACCAGCCATGCCGCGTGCAGGATGACGGCTCTATGAGTTGTAAACTGCTTTTGTATTAGGGTAAACACCCCTACGTGTAGGGGCTTGAAAGTATAATACGAATAAGGATCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGATCCAAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGTAGGCGGTTTGATA 
    ##                                                                                                                                                                                                                                                         30 
    ## GACTACCGGGGTATCTAATCCTGTTCGATACCCACGCTTTCGTGCCTCAGCGTCAGTAGGGCGCCGGTATGCTGCCTTCGCAATCGGAGTTCTGCGTGATATCTATGCATTTCACCGCTACACCACGCATTCCGCATACTTCTCGCCCACTCAAGAACGCCAGTTTCAACGGCTCGACGGCGTTGAGCACCGCTTTTTTACCGCTGACTTGGCATCCCGCCTACGCACCCTTTAAACCCAATAAATCCGG 
    ##                                                                                                                                                                                                                                                        107 
    ## GACTACTGGGGTATCTAATCCTGTTCGATACCCACGCTTTCGTGCCTCAGCGTCAGTGTGGAGCCGGTATGCTGCCTTCGCAATCGGGGTTCTTCGTGATATCTATGCATTTCACCGCTACACCACGAATTCCGCATACTTCTCTCCAACTCGAGCACGCCAGTTTCAACGGCCGGCCGGGGTTGAGCCCCGACATTTTACCGCTGACTTGACACGCCGCCTGCGCACCCTTTAAACCCAATAAATCCGG 
    ##                                                                                                                                                                                                                                                         73 
    ## TGGGGAATATTGCACAATGGGGGGAACCCTGATGCAGCGACGCCGCGTGAGTGAAGAAGTATTTCGGTACGTAAAGCTCTATCAGCAGGGAAGAAAAAAGGACGGTACCTGAGTAAGAAGCCCCGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGGGGCAAGCGTTATCCGGAATTACTGGGTGTAAAGGGAGCGTAGACGGCAATGCAAGTCCGGAGTGGAATGCGGCAGCTCAACTGCCGAA 
    ##                                                                                                                                                                                                                                                         10 
    ## TAGGGAATCTTTCACAATGGGCGAAAGCCTGATGGAGCAACGCCGCGTGCAGGATGAAGGCCTTCGGGTTGTAAACTGCTTTTATAAGCGAGAAATATGATGGTAACTTATGAATAAGGATCGGCTAACTACGTGCCAGCAGCCGCGGTCATACGTAGGATCCGAGCATTATCCGGAGTGACTGGGTGTAAAGAGTTGCGTAGGTGGCATAATAAGTAGCTAGTGAAATCTGGTGGCTCAACCATTCAGA 
    ##                                                                                                                                                                                                                                                         11 
    ## GACTACAGGGGTATCTAATCCTGTTCGATCCCCGCGCTTTCGTGCCTCAGCGTCAGTGGCGGCACGGCACGCTGCCTTCGCAATCGGGGTTCTGCGCGATATCTATGCATTTCACCGCTACACCGCGCATTCCGCATGCCTCCGCCGCACTCAAGGCCCCCAGTTCCGACGGCTCGGCGGGGTTGAGCCCCGCGATTTCACCGCCGGCTTAAAGGCCCGCCTGCGCACCCTTTAAACCCAATGAATCCGG 
    ##                                                                                                                                                                                                                                                         80 
    ## TGGGGGATATTGGACAATGGGGGAAACCCTGATCCAGCGATGCCGCGTGAGGGAAGAAGGTTTTCGGATTGTAAACCTCTGTGGAGGGGGGCGATAATGACGGTACCCCTTTAGGAAGCCACGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGTGGCGAGCGTTGTCCGGAATTACTGGGTGTAAAGGGAGTGTAGGCGGGAAAGCAAGTCAGAAGTGAAAACTATGGGCTTAACCCATAGCCT 
    ##                                                                                                                                                                                                                                                         19 
    ## GACTACTGGGGTATCTAATCCTGTTTGCTCCCCACGCTTTCGAGCCTCAGCGTCAGTAAAAGCCCAGCAAGCCGCCTTCGCCACTGGTGTTCCTCCTAATATCTACGCATTTCACCGCTACACTAGGAATTCCGCTTGCCTCTACTTCACTCAAGAAAAACAGTTTCAAAAGCAGGCTATGGGTTAAGCCCATAGTTTTCACTTCTGACTTGCTTTCCCGCCTACACTCCCTTTACACCCAGTAATTCCG 
    ##                                                                                                                                                                                                                                                         51 
    ## GACTACTGGGGTATCTAATCCTGTTTGCTCCCCACGCTTTCGAGCCTCAACGTCAGTTTCTGTCCAGCAGGCCGCCTTCGCCACTGGTGTTCCTCCTAATATCTACGCATTTCACCGCTACACTAGGAATTCCGCCTGCCTCTCCAGTACTCAAGCTATACAGTTTCCAAAGCAGTTCCGCGGTTGAGCCACGGGCTTTCACTTCAGACTTGCACAGCCGTCTACGCTCCCTTTACACCCAGTAAATCCG 
    ##                                                                                                                                                                                                                                                         95 
    ## GACTACAGGGGTATCTAATCCTGTTCGATACCCACGCTTTCGTGCATGAGCGTCAGTTGCGCGCCGGTAGGCTGCCTTCGCGATCGGAGTTCTGCGTGATATCTATGCATTTCACCGCTACACCACGCATTCCGCCTACCCTTCGCGCACTCAAGGGACACAGTTTCAACGGCGGAGCGGGGTTGAGCCCCGCGATTTTACCGCTGACTTGTGCCTCCGCCTGCGCACCCTTTAAACCCAATAAATCCGG 
    ##                                                                                                                                                                                                                                                         60 
    ## GACTACTGGGGTATCTAATCCTGTTTGCTACCCACGCTTTCGTGCTTCAGCGTCAGTTAAAGCCCAGTAGGCCGCCTTCGCCACTGGTGTTCCTCCCGATCTCTACGCATTTCACCGCTACACCGGGAATTCCGCCTACCTCTACTTCACTCAAGATAAACAGTTTCAAAAGCAGTTCATGGGTTAAGCCCATGGATTTCACTTCTGACTTGCTTACCCGCCTACGCACCCTTTACACCCAGTAAATCCG 
    ##                                                                                                                                                                                                                                                         65 
    ## CCTACGGGCGGCTGCAGTGGGGAATATTGGGCAATGGGCGCAAGCCTGACCCAGCAACGCCGCGTGAAGGAAGAAGGCTTTCGGGTTGTAAACTTCTTTTGACAGGGAAGAGTAGAAGACGGTACCTGTCGAATAAGCCACGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGTGGCAAGCGTTGTCCGGATTTACTGGGTGTAAAGGGCGTGTAGCCGGGTTGACAAGTCAGATGTGAAATCCT 
    ##                                                                                                                                                                                                                                                         15 
    ## GACTACCGGGGTATCTAATCCTGTTTGATCCCCACGCTTTCGTGCTTCAGTGTCAGTTATAGTTTAGTAAGCTGCCTTCGCAATCGGAGTTCTGCGTGATATCTATGCATTTCACCGCTACACCACGCATTCCGCCTACCTCAAATATACTCAAGTCAACCAGTTTCAACGGCAATTTTATGGTTGAGCCACAAACTTTCACCGCTGACTTAATTAACCACCTACGCACCCTTTAAACCCAATAAATCCG 
    ##                                                                                                                                                                                                                                                         35 
    ## TGGGGAATATTGGGCAATGGAGGCAACTCTGACCCAGCAATGCCGCGTGAATGAAGAAGGTCTTCGGATTGTAAAGTTCTTTAATCAGGGACGAAGAAGATGACGGTACCTGAGGAATAAGCCACGGCAAACTATGTGCCAGCAGCCGCGGTAATACATAGGTGGCAAGCGTTGTCCGGAATGACTGGGCGTAAAGGGAGCGTAGGCGGTCTGTTAAGTTGGAAGTGAAATCCCGGGGCTCAACCCCGGA 
    ##                                                                                                                                                                                                                                                          6 
    ## CCTACGGGCGGCTGCAGTGGGGAATATTGCGCAATGGGGGCAACCCTGACGCAGCAACGCCGCGTGATTGAAGAAGGCCTTCGGGTTGTAAAGATCTTTAATCGGGGACGAATTTTGACGGTACCCGAAGAATAAGCTCCGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGGAGCAAGCGTTATCCGGATTTACTGGGTGTAAAGGGCGTGTAGGCGGGCATGCAAGCCAGAAGTGAAATCTGG 
    ##                                                                                                                                                                                                                                                         10 
    ## CCTACGGGCGGCTGCAGTAGGGAATCTTTCACAATGGGCGAAAGCCTGATGGAGCAACGCCGCGTGCAGGATGAAGGCCTTCGGGTTGTAAACTGCTTTTATAAGCGAGAAATATGATGGTAACTTATGAATAAGGATCGGCTAACTACGTGCCAGCAGCCGCGGTCATACGTAGGATCCGAGCATTATCCGGAGTGACTGGGTGTAAAGAGTTGCGTAGGTGGCATAATAAGTAGCTAGTGAAATCTGG 
    ##                                                                                                                                                                                                                                                         12 
    ## GACTACAGGGGTATCTAAGCCCGTTTGCTCCCTACGCTTTCGTGCCTTAGTGTCAGAACCGTCCCAGTAACCTGCCTACGCTATTGGTGTTCTTTCTAATATCTACGGATTTCACTCCTACACTAGAAATTCCAGTTACCCCTAACGGTCTCGAGCTTAACAGTTTAGCTAATAGTCTGAATGGTTGAGCCACCAGATTTCACTAGCTACTTATTATGCCACCTACGCAACTCTTTACACCCAGTCACTC 
    ##                                                                                                                                                                                                                                                         29 
    ## CCTACGGGGGGCTGCAGTGAGGAATATTGGGCAATGGAGGCAACTCTGACCCAGCCATGCCGCGTGAGTGAAGAAGGTTTTCGGATTGTAAAGCTCTTTCGGATGTGACGATGATGACGGTAGCATCTAAAGAAGCCCCGGCTAACTTCGTGCCAGCAGCCGCGGTAATACGAAGGGGGCAAGCGTTGTTCGGAATTACTGGGCGTAAAGGGAGTGTAGGCGGTTTAGTAAGATAGTGGTGAAATCCCAG 
    ##                                                                                                                                                                                                                                                          6 
    ## TGGGGAATATTGGGCAATGGGCGCAAGCCTGACCCAGCAACGCCGCGTGAAGGAAGAAGGCCCTCGGGTTGTAAACTTCTTTTATTCGGGACGAAACAAATGACGGTACCGAATGAATAAGCCACGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGTGGCAAGCGTTATCCGGATTTACTGGGTGTAAAGGGCGTGTAGGCGGGAGAGCAAGTCAGACGTGAAATTCCAGGGCTCAACCCTGGA 
    ##                                                                                                                                                                                                                                                         18 
    ## TGGGGAATATTGCACAATGGGCGAAAGCCTGATGCAGCAACGCCGCGTGAAGGAAGAAGGCCTTCGGGTCGTAAACTTCTGTCCTTGGGGAAGAAGAACTGACGGTACCCAAGGAGGAAGCCCCGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGGGGCAAGCGTTATCCGGAATTATTGGGCGTAAAGAGTACGTAGGTGGCAACCTAAGCGCAGGGTTTAAGGCAATGGCTCAACCATTGTT 
    ##                                                                                                                                                                                                                                                          5 
    ## CCTACGGGGGGCTGCAGTGGGGAATATTGCGCAATGGACGAAAGTCTGACGCAGCGACGCCGCGTGAGGGATGAAGGTTTTCGGATCGTAAACCTCTGTCAGGAGGGAAGAAAGCGCATGGCGCCAATCAGCCATGCGTTGACGGTACCTCCAAAGGAAGCGCCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGGCGCGAGCGTTAATCGGAATTACTGGGCGTAAAGCGCGCGTAGGCTGTT 
    ##                                                                                                                                                                                                                                                          5 
    ## GACTACCCGGGTATCTAATCCCGTTCGCTACCCACGCTTTCGCGCCTCAGCGTCAGGTACAGTCCAGAAAGTCGCCTTCGCCACTGGTGTTCCTCCTAATATCTACGCATTTCACCGCTACACTAGGAATTCCACTTTCCTCTCCTGCCCTCAAGACACCCAGTTTCAAGTGCAGTCTCCAGGTTGAGCCCGGAGTTTTCACACCTGACTTAAGCGCCCACCTACGCGCCCTTTACGCCCAGTAATTCCG 
    ##                                                                                                                                                                                                                                                          4 
    ## GACTACTCGGGTATCTAATCCTGTTTGATCCCCATGCTTTCGAGCCTCAGCGTCAATGTACAGCCAGATAGTCGCCTTCGCCACTGGTGTTCTACCCAATATCTACGAATTTCACCTCTACACTGGGTATTCCACTATCCTCTCTGTAATTCTAGTCTAATAGTTATAATGGCAATTCCAAAGTTAAGCTCTGGGATTTCACCACTATCTTACTAAACCGCCTACACTCCCTTTACGCCCAGTAATTCCG 
    ##                                                                                                                                                                                                                                                         11 
    ## CCTACGGGCGGCTGCAGTGGGGAATCTTCCGCAATGGGCGCAAGCCTGACGGAGCAATGCCGCGTGAACGAAGAAGGTCTTCGGATTGTAAAGTTCTGTCCTTATCGAAGAGAGGGTATAGAGTGAAAAATGATATACTAGGACGGTAGATGAGGAGGAAGCCCCGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGGGGCGAGCGTTGTCCGGAATCACTGGGCGTAAAGGGCGCGTAGGCGGT 
    ##                                                                                                                                                                                                                                                          2 
    ## GACTACTAGGGTATCTAATCCTGTTTGCTCCCCACACTTTCGTGCCTCAGTGTCAGTAACAGTCCAGAAAGCCGCCTTCGCCACTGGTGTTCCTCCTAATATCTACGCATTTCACCGCTACACTAGGAATTCCGCTTTCCTCTCCTGCACTCAAGTAACCCAGTTCGCAGGGCGAACAATGGTTGAGCCATTGCCTTAAACCCCGCGCTTAAGTAACCACCTACGCACTCTTTACGCCCAATAATTCCGG 
    ##                                                                                                                                                                                                                                                         12 
    ## GACTACAAGGGTATCTAATCCTGTTTGCTCCCCACGCTTTCGCGCCTCAGCGTCAATACCGGTCCAGGTGGCCGCCTTCGCCACTGATGTTCCTCCAGATATCTACGGATTTCACTCCTACACCTGGAATTCCGCCACCCTCTCCCGGATTCAAGTCGCGCAGTATCAAGGGCTGTTCCACGGTTGAGCCGTGGGATTTCACCCCCAACTTACACAACAGCCTACGCGCGCTTTACGCCCAGTAATTCCG 
    ##                                                                                                                                                                                                                                                          8 
    ## CCTACGGGGGGCTGCAGTAGGGAATTTTCGGCAATGGGGGGAACCCTGACCGAGCAACGCCGCGTGAAGGAAGAAGGAATTCGTTCTGTAAACTTCTGTTATAAAGGAAGAAAGACGGATGGAGGAAATGACATCCGAGTGACGGTACTTTATGAGAAAGCCACGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGTGGCAAGCGTTATCCGGAATTATTGGGCGTAAAGAGGGAGCAGGCGGCA 
    ##                                                                                                                                                                                                                                                          3 
    ## GACTACAAGGGTATCTAATCCTGTTCGATACCCACACTTTCGTGCATGAGCGTCAGTTGAGCGCCGGTATGCTGCCTTCGCAATCGGAGTTCTGCGTGATATCTATGCATTTCACCGCTACACCACGCATTCCGCATACTTCTCGCTCACTCAAGAAAACCAGTTTCAACGGCTCGAAGAGGTTGAGCCTCTCAATTTTACCGCTGACTTGATCTTCCGCCTGCGCACCCTTTAAACCCAATAAATCCGG 
    ##                                                                                                                                                                                                                                                        108 
    ## GACTACCAGGGTATCTAATCCTGTTCGATACCCACACTTTCGTGCATGAGCGTCAGTTGAGCGCCGGTATGCTGCCTTCGCAATCGGAGTTCTGCGTGATATCTATGCATTTCACCGCTACACCACGCATTCCGCATACTTCTCGCTCACTCAAGAAAACCAGTTTCAACGGCTCGAAGAGGTTGAGCCTCTCAATTTTACCGCTGACTTGATCTTCCGCCTGCGCACCCTTTAAACCCAATAAATCCGG 
    ##                                                                                                                                                                                                                                                        120 
    ## GACTACAGGGGTATCTAATCCTGTTCGATACCCACACTTTCGTGCATGAGCGTCAGTTGAGCGCCGGTATGCTGCCTTCGCAATCGGAGTTCTGCGTGATATCTATGCATTTCACCGCTACACCACGCATTCCGCATACTTCTCGCTCACTCAAGAAAACCAGTTTCAACGGCTCGAAGAGGTTGAGCCTCTCAATTTTACCGCTGACTTGATCTTCCGCCTGCGCACCCTTTAAACCCAATAAATCCGG 
    ##                                                                                                                                                                                                                                                        200 
    ## TGAGGAATATTGGTCAATGGGCGCGAGCCTGAACCAGCCAAGTCGCGTGAGGGATGACGGTCCTACGGATTGTAAACCTCTTTTGCCGGGGAGCAAGCGTGCGTTCGTGAACGCGCGTCGAGAGTACCCGGAGAAAAAGCATCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGTAGGCGGGATGCCAAGTCAGCGGTAAAAAA 
    ##                                                                                                                                                                                                                                                         26 
    ## GACTACTCGGGTATCTAATCCTGTTTGCTCCCCACGCTTTCGAGCCTCAACGTCAGTTACCGTCCAGTAAGCCGCCTTCGCCACTGGTGTTCCTCCTAATATCTACGCATTTCACCGCTACACTAGGAATTCCGCTTACCTCTCCGGTACTCCAGCAAAAAAGTTTCCAAAGCAGTTCCGCGGTTGAGCCGCGGGCTTTCACTTCAGACTTTCTTCGCCGTCTACGCTCCCTTTACACCCAGTAAATCCG 
    ##                                                                                                                                                                                                                                                         40 
    ## CCTACGGGCGGCTGCAGTGAGGAATATTGGTCAATGGACGAGAGTCTGAACCAGCCAAGTAGCGTGAAGGATGACTGCCCTATGGGTTGTAAACTTCTTTTATATGGGAATAAAATGTTCCACGTGTGGGATTTTGTATGTACCATATGAATAAGGATCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGATCCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGAGCGTAGGTGGATTGTTA 
    ##                                                                                                                                                                                                                                                        186 
    ## TGGGGAATATTGCACAATGGGGGGAACCCTGATGCAGCAATGCCGCGTGGGTGAAGAAGTACCCCGGTATGTAAAGCCCTGTCAGCAGGGAAGAAAATGACGGTACCTGACTAAGAAGCCCCGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGGGGCGAGCGTTATCCGGATTTACTGGGTGTAAAGGGAGCGCAGGCGGCGTGGCAAGTCAGATGTGAAAACCCGGGGCCCAACCCCGGGACT 
    ##                                                                                                                                                                                                                                                         21 
    ## GACTACAGGGGTATCTAATCCTGTTTGATCCCCGCACTTTCGTGCCTCAGCGTCAGTCGGGCGCCGGTATGCTGCCTTCGCAATCGGGGTTCTGCGTGATATCTATGCATTTCACCGCTACACCACGCATTCCGCATACTTCTCGCCCACTCGAGCCCGGCAGTTTCAACGGCTGTACGGGGTTGAGCCCCGCAATTTTACCGCTGACTTGGCAGGCCGCCTACGCACCCTTTAAACCCAATAAATCCGG 
    ##                                                                                                                                                                                                                                                         97 
    ## GACTACCCGGGTATCTAATCCCTTTCGCTCCCCTGGCCTTCGTGCCTCAGCGTCAGTTAATGTCCAGGAACCCGCCTTCGCCACGAGTGTTCCTCTCGATATCTACGCATTTCACTGCTACACCGAGAATTCCGGTTCCCCCTCCATTACTCTAGTCTCGCAGTATCATGTGCCGTCCGCGGGTTGAGCCCGCGCCTTTCACACACGACTTACGAAACAGCCTACGCACGCTTTACGCCCAGTGATTCCG 
    ##                                                                                                                                                                                                                                                        101 
    ## TGGGGAATATTGCACAATGGGGGAAACCCTGATGCAGCGACGCCGCGTGAGCGAAGAAGTATTTCGGTATGTAAAGCTCTATCAGCAGGAAAGAAAATGACGGTACCTGACTAAGAAGCTCCGGCTAAATACGTGCCAGCAGCCGCGGTAATACGTATGGAGCAAGCGTTATCCGGATTTACTGGGTGTAAAGGGAGCGCAGGCGGTCTGGCAAGTCTGGTGTGAAAGGCCGGGGCTCAACCCCGGGACT 
    ##                                                                                                                                                                                                                                                         12 
    ## CCTACGGGCGGCTGCAGTGAGGAATATTGGTCAATGGCCGGAAGGCTGAACCAGCCAAGTCGCGTGAGGGAATAAGGCCCTACGGGTCGTAAACCTCTTTTGCCAGGGAGCAATGCCGTCCACGTGTGGACGGAAGGAGAGTACCTGGAGAAAAAGCATCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGTAGGCGGCCTGCC 
    ##                                                                                                                                                                                                                                                        161 
    ## GACTACCAGGGTATCTAATCCTGTTTGATACCCACACTTTCGAGCCTCAGTGTCAGTTGCAGTCCAGTGAGCTGCCTTCGCAATCGGAGTTCTTCGTGATATCTAAGCATTTCACCGCTACACCACGAATTCCGCCCACCTCTACTGTACTCAAGACTGCCAGTTTCAACTGCAATTTTACGGTTGAGCCGCAAACTTTCACAACTGACTTAACAATCCACCTACGCTCCCTTTAAACCCAATAAATCCG 
    ##                                                                                                                                                                                                                                                        122 
    ## GACTACACGGGTATCTAATCCTGTTTGATACCCACACTTTCGAGCCTCAGTGTCAGTTGCAGTCCAGTGAGCTGCCTTCGCAATCGGAGTTCTTCGTGATATCTAAGCATTTCACCGCTACACCACGAATTCCGCCCACCTCTACTGTACTCAAGACTGCCAGTTTCAACTGCAATTTTACGGTTGAGCCGCAAACTTTCACAACTGACTTAACAATCCACCTACGCTCCCTTTAAACCCAATAAATCCG 
    ##                                                                                                                                                                                                                                                         84 
    ## GACTACCAGGGTATCTAATCCTGTTTGCTCCCCACACTTTCGTGCCTCAACGTCAGTTACAGTCCAGAAAGTCGCCTTCGCCACCGGTGTTCTTCCTAATATCTACGCATTTCACCGCTACACTAGGAATTCCACTTTCCTCTCCTGCACTCAAGAATAATAGTTTTGGTTGCAGTTCCTCAGTTGAGCCGAGGGGTTTCACAACCAACTTGTTATCCCGTCTACGCACCCTTTACACCCAATAATTCCG 
    ##                                                                                                                                                                                                                                                         40 
    ## GACTACAAGGGTATCTAATCCCTTTCGCTCCCCTGGCCTTCGTGCCTCAGCGTCAGTTAATGTCCAGGAACCCGCCTTCGCCACGAGTGTTCCTCTCGATATCTACGCATTTCACTGCTACACCGAGAATTCCGGTTCCCCCTCCATTACTCTAGTCTCGCAGTATCATGTGCCGTCCGCGGGTTGAGCCCGCGCCTTTCACACACGACTTACGAAACAGCCTACGCACGCTTTACGCCCAGTGATTCCG 
    ##                                                                                                                                                                                                                                                        107 
    ## GACTACCGGGGTATCTAATCCTGTTCGATACCCACACTTTCGTGCATGAGCGTCAGTTGAGCGCCGGTATGCTGCCTTCGCAATCGGAGTTCTGCGTGATATCTATGCATTTCACCGCTACACCACGCATTCCGCATACTTCTCGCTCACTCAAGAAAACCAGTTTCAACGGCTCGAAGAGGTTGAGCCTCTCAATTTTACCGCTGACTTGATCTTCCGCCTGCGCACCCTTTAAACCCAATAAATCCGG 
    ##                                                                                                                                                                                                                                                        195 
    ## CCTACGGGCGGCTGCAGTGGGGAATATTGCACAATGGGGGAAACCCTGATGCAGCGACGCCGCGTGAGCGAAGAAGTATTTCGGTATGTAAAGCTCTATCAGCAGGAAAGAAAATGACGGTACCTGACTAAGAAGCTCCGGCTAAATACGTGCCAGCAGCCGCGGTAATACGTATGGAGCAAGCGTTATCCGGATTTACTGGGTGTAAAGGGAGCGCAGGCGGTCTGGCAAGTCTGGTGTGAAAGGCCGG 
    ##                                                                                                                                                                                                                                                         16 
    ## GACTACCGGGGTATCTAATCCTGTTTGATACCCACACTTTCGAGCCTCAGTGTCAGTTGCAGTCCAGTGAGCTGCCTTCGCAATCGGAGTTCTTCGTGATATCTAAGCATTTCACCGCTACACCACGAATTCCGCCCACCTCTACTGTACTCAAGACTGCCAGTTTCAACTGCAATTTTACGGTTGAGCCGCAAACTTTCACAACTGACTTAACAATCCACCTACGCTCCCTTTAAACCCAATAAATCCG 
    ##                                                                                                                                                                                                                                                        171 
    ## TGAGGAATATTGGTCAATGGGCGCGAGCCTGAACCAGCCAAGTCGCGTGAGGGAGGACGGCCCTACGGGTTGTAAACCTCTTTTGCCGGGGAGCAACGGGCGCCACGTGTGGCGCCACTGAGAGTACCCGGAGAAAAAGCATCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGTAGGCGGATCGTTAAGTCAGTGGTCAAATT 
    ##                                                                                                                                                                                                                                                          9 
    ## TGGGGAATATTGCACAATGGGGGAAACCCTGATGCAGCAACGCCGCGTGAGTGAAGAAGTATTTCGGTATGTAAAGCTCTATCAGCAGGGAAGAGAATGACGGTACCTGACTAAGAAGCCCCGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGGGGCGAGCGTTATCCGGATTTACTGGGTGTAAAGGGAGCGTAGACGGTTCTGCAAGTCTGAAGTGAAAGCCCGTGGCTTAACCGCGGAACG 
    ##                                                                                                                                                                                                                                                         14 
    ## GACTACTGGGGTATCTAATCCTGTTCGATACCCACGCTTTCGTGCCTCAGCGTCAGTCGGGCGCCGGTATGCTGCCTTCGCAATCGGAGTTCTGCGTGATATCTATGCATTTCACCGCTACACCACGCATTCCGCATACTTCTCGCCCACTCAAGACTGCCAGTTTCAACGGCTGTACGGCGTTGAGCACCGCATTTTTACCGCTGACTTGACCGTCCGCCTACGCACCCTTTAAACCCAATAAATCCGG 
    ##                                                                                                                                                                                                                                                        151 
    ## CCTACGGGGGGCTGCAGTGAGGAATATTGGTCAATGGGCGGAAGCCTGAACCAGCCAAGTCGCGTGAGGGACTAAGGCCCTACGGGTCGTAAACCTCTTTTGCCGGGGATCAGTGCCCAGCTCGCGAGCTGGGAGGGAGCGTACCCGGAGAAAAAGCATCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGTAGGCGGGCCTTC 
    ##                                                                                                                                                                                                                                                         18 
    ## GACTACAGGGGTATCTAATCCTGTTTGCTCCCCACGCTTTCGAGCCTCAGCGTCAGTAAAAGCCCAGTAAGCCGCCTTCGCCACCGATGTTCCTCCTAATATCTACGCATTTCACCGCTACACTAGGAATTCCGCTTACCTCTACTTCACTCAAGAACAGCAGTTTCAAATGCAGTCTGAAGGTTGAGCCCTCAGTTTTCACATCTGACTTGCCGTCCCGCCTACGCTCCCTTTACACCCAGTAATTCCG 
    ##                                                                                                                                                                                                                                                         55 
    ## GACTACAGGGGTATCTAATCCTGTTTGATACCCACACTTTCGAGCCTCAGTGTCAGTTGCAGTCCAGTGAGCTGCCTTCGCAATCGGAGTTCTTCGTGATATCTAAGCATTTCACCGCTACACCACGAATTCCGCCCACCTCTACTGTACTCAAGACTGCCAGTTTCAACTGCAATTTTACGGTTGAGCCGCAAACTTTCACAACTGACTTAACAATCCACCTACGCTCCCTTTAAACCCAATAAATCCG 
    ##                                                                                                                                                                                                                                                        170 
    ## GACTACAGGGGTATCTAATCCCTTTCGCTCCCCTGGCCTTCGTGCCTCAGCGTCAGTTAATGTCCAGGAACCCGCCTTCGCCACGAGTGTTCCTCTCGATATCTACGCATTTCACTGCTACACCGAGAATTCCGGTTCCCCCTCCATTACTCTAGTCTCGCAGTATCATGTGCCGTCCGCGGGTTGAGCCCGCGCCTTTCACACACGACTTACGAAACAGCCTACGCACGCTTTACGCCCAGTGATTCCG 
    ##                                                                                                                                                                                                                                                        119 
    ## CCTACGGGGGGCTGCAGTGGGGAATATTGCACAATGGGGGGAACCCTGATGCAGCAATGCCGCGTGGGTGAAGAAGTACCCCGGTATGTAAAGCCCTGTCAGCAGGGAAGAAAATGACGGTACCTGACTAAGAAGCCCCGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGGGGCGAGCGTTATCCGGATTTACTGGGTGTAAAGGGAGCGCAGGCGGCGTGGCAAGTCAGATGTGAAAACCCGG 
    ##                                                                                                                                                                                                                                                         20 
    ## GACTACTGGGGTATCTAAGCCTGTTTGCTCCCCACGCTTTCGCGCCTCAGCGTCAGTTGCTGTCCAGCAGACCGCCTTCGCCACTGGTGTTCCTCCTAATATCTACGCATTTCACCGCTACACTAGGAATTCCGTCTGCCTCTCCAGTACTCAAGAAAAACAGTTTCAAATGCAGGCCACAGGTTGAGCCCATGGTTTTCACATCTGACTTGCTTTCCCGCCTACACGCCCTTTACACCCAGTAAATCCG 
    ##                                                                                                                                                                                                                                                         17 
    ## GACTACTAGGGTATCTAATCCTGTTCGATACCCACACTTTCGTGCATGAGCGTCAGTTGAGCGCCGGTATGCTGCCTTCGCAATCGGAGTTCTGCGTGATATCTATGCATTTCACCGCTACACCACGCATTCCGCATACTTCTCGCTCACTCAAGAAAACCAGTTTCAACGGCTCGAAGAGGTTGAGCCTCTCAATTTTACCGCTGACTTGATCTTCCGCCTGCGCACCCTTTAAACCCAATAAATCCGG 
    ##                                                                                                                                                                                                                                                        137 
    ## TGAGGAATATTGGTCAATGGACGAGAGTCTGAACCAGCCAAGTAGCGTGAAGGATGACTGCCCTATGGGTTGTAAACTTCTTTTATATGGGAATAAAACAGGGTATGCATACCCTCTTGTATGTACCATATGAATAAGGATCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGATCCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGAGCGTAGGTGGATTGTTAAGTCAGTTGTGAAAGTT 
    ##                                                                                                                                                                                                                                                          8 
    ## CCTACGGGTGGCTGCAGTGAGGAATATTGGTCAATGGGCGCGAGCCTGAACCAGCCAAGTCGCGTGAGGGAAGACGGTCCTATGGATTGTAAACCTCTTTTGTCGGGGAGCAAAAGACGCCACGCGTGGCGTTCCGAGAGTACCCGAAGAAAAAGCATCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGCAGGCGGAAGATCA 
    ##                                                                                                                                                                                                                                                        135 
    ## TGGGGAATATTGGGCAATGGAGGAAACTCTGACCCAGCAACGCCGCGTGAAGGATGAAGGTTTTCGGATTGTAAACTTCTGTTGCGAGGGATGAAGAATGACGGTACCTCGTGAGGAAGCCCCGGCAAACTACGTGCCAGCAGCCGCGGTAATACGTAGGGGGCAAGCGTTGTCCGGAATGACTGGGCGTAAAGGGAGTGTAGGCGGCCTTTTAAGTTATATGTGAAAGCCCCCGGCTTAACCGAGGAAC 
    ##                                                                                                                                                                                                                                                          4 
    ## GACTACTCGGGTATCTAATCCTGTTCGATACCCACACTTTCGTGCATGAGCGTCAGTTGAGCGCCGGTATGCTGCCTTCGCAATCGGAGTTCTGCGTGATATCTATGCATTTCACCGCTACACCACGCATTCCGCATACTTCTCGCTCACTCAAGAAAACCAGTTTCAACGGCTCGAAGAGGTTGAGCCTCTCAATTTTACCGCTGACTTGATCTTCCGCCTGCGCACCCTTTAAACCCAATAAATCCGG 
    ##                                                                                                                                                                                                                                                        136 
    ## TGAGGAATATTGGTCAATGGGCGGAAGCCTGAACCAGCCAAGTCGCGTGAGGGATTAAGGCCCTGCGGGTCGTAAACCTCTTTTGCCGGGGAGCAGTGGTCCGGACGCGTCCGGGCCGGGAGAGTACCCGGAGAAAAAGCATCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGTAGGCGGGCATTTAAGCCGGCGGTCAAATG 
    ##                                                                                                                                                                                                                                                          7 
    ## CCTACGGGGGGCTGCAGTGAGGAATATTGGTCAATGGGCGGAAGCCTGAACCAGCCAAGTCGCGTGAGGGATTAAGGCCCTACGGGTCGTAAACCTCTTTTGTCAGGGAGCAAATGCGCCCACGTGTGGGCGAGTGGAGAGTACCTGAAGAAAAAGCATCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGTAGGCGGGTTGTT 
    ##                                                                                                                                                                                                                                                         14 
    ## GACTACCGGGGTATCTAAGCCTGTTTGCTCCCCACGCTTTCGAGCCTCAACGTCAGTTACAGTCCAGCAGGCCGCCTTCGCCACTGGTGTTCCTCCTAATATCTACGCATTTCACCGCTACACTAGGAATTCCGCCTGCCTCTCCTGCACTCCAGTTACACAGTTTCCAGAGCAGTCCGGGGGTTGGGCCCCCGCCTTTCACTCCAGACTTGCATTACCGTCTACGCTCCCTTTACACCCAGTAATTCCG 
    ##                                                                                                                                                                                                                                                         32 
    ## GACTACCGGGGTATCTAATCCCTTTCGCTCCCCTGGCCTTCGTGCCTCAGCGTCAGTTAATGTCCAGGAACCCGCCTTCGCCACGAGTGTTCCTCTCGATATCTACGCATTTCACTGCTACACCGAGAATTCCGGTTCCCCCTCCATTACTCTAGTCTCGCAGTATCATGTGCCGTCCGCGGGTTGAGCCCGCGCCTTTCACACACGACTTACGAAACAGCCTACGCACGCTTTACGCCCAGTGATTCCG 
    ##                                                                                                                                                                                                                                                        124 
    ## CCTACGGGGGGCTGCAGTGAGGAATATTGGTCAATGGACGAGAGTCTGAACCAGCCAAGTAGCGTGAAGGATGACTGCCCTATGGGTTGTAAACTTCTTTTATATGGGAATAAAACAGGGTATGCATACCCTCTTGTATGTACCATATGAATAAGGATCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGATCCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGAGCGTAGGTGGATTGTTA 
    ##                                                                                                                                                                                                                                                         18 
    ## CCTACGGGCGGCTGCAGTGGGGAATATTGCACAATGGGCGAAAGCCTGATGCAGCAACGCCGCGTGAAGGAAGAAGGCCTTCGGGTCGTAAACTTCTGTCCTTGGGGAAGAAGAACTGACGGTACCCAAGGAGGAAGCCCCGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGGGGCAAGCGTTATCCGGAATTATTGGGCGTAAAGAGTACGTAGGTGGCAACCTAAGCGCAGGGTTTAAGGCA 
    ##                                                                                                                                                                                                                                                         10 
    ## CCTACGGGTGGCTGCAGTGAGGAATATTGGTCAATGGCCGGGAGGCTGAACCAGCCAAGTCGCGTGAGGGAAGACGGCCCTACGGGTTGTAAACCTCTTTTGTCAGGGAGCAAGAGGCGGGTCGGGACCTGCTGTGAGAGTACCTGAAGAAAAAGCATCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGCAGGCGGAGGCACA 
    ##                                                                                                                                                                                                                                                         34 
    ## GACTACTAGGGTATCTAATCCTGTTTGCTACCCATGCTTTCGAGCCTCAGCGTCAGTTAGTGCCCAGTAGGCCGCCTTCGCCACTGGTGTTCCTCCCGATATCTACGCATTCCACCGCTACACCGGGAATTCCGCCTACCTCTACACCACTCAAGACTGACAGTTTTGAAAGCAGGTTACGGGTTGAGCCCGTAGTTTTCACTCCCAACTTGTCAATCCGCCTGCGCTCCCTTTACACCCAGTAATTCCG 
    ##                                                                                                                                                                                                                                                          9 
    ## GACTACAGGGGTATCTAATCCTGTTTGCTCCCCACGCTTTCGTGCCTCAACGTCAGATATAGATTGGTAAGCTGCCTTCGCAATCGGTGTTCTGTATGATCTCTAAGCATTTCACCGCTACACCATACATTCCGCCTACCGCATCTACTCTCTAGATCAACAGTATTAGAGGCAGTTCCCGAGTTGAGCTCGGACATTTCACCTCTAACTTATCAATCCGCCTACGCACCCTTTAAACCCAATAAATCCG 
    ##                                                                                                                                                                                                                                                         33 
    ## GACTACCCGGGTATCTAATCCTGTTTGATCCCCGCACTTTCGTGCCTCAGCGTCAGTCGGGCGCCGGTATGCTGCCTTCGCAATCGGGGTTCTGCGTGATATCTATGCATTTCACCGCTACACCACGCATTCCGCATACTTCTCGCCCACTCGAGCCCGGCAGTTTCAACGGCTGTACGGGGTTGAGCCCCGCAATTTTACCGCTGACTTGGCAGGCCGCCTACGCACCCTTTAAACCCAATAAATCCGG 
    ##                                                                                                                                                                                                                                                         53 
    ## CCTACGGGCGGCTGCAGTGAGGAATATTGGTCAATGGACGCAAGTCTGAACCAGCCATGCCGCGTGCAGGAAGACGGCTCTATGAGTTGTAAACTGCTTTTGTACGAGGGTAAACGCAGATACGTGTATCTGTCTGAAAGTATCGTACGAATAAGGATCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGATTCAAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGTAGGCGGTCGGATA 
    ##                                                                                                                                                                                                                                                        142 
    ## CCTACGGGGGGCTGCAGTGGGGAATATTGGGCAATGGGCGCAAGCCTGACCCAGCAACGCCGCGTGAAGGAAGAAGGCCCTCGGGTTGTAAACTTCTTTTATTCGGGACGAAACAAATGACGGTACCGAATGAATAAGCCACGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGTGGCAAGCGTTATCCGGATTTACTGGGTGTAAAGGGCGTGTAGGCGGGAGAGCAAGTCAGACGTGAAATTC 
    ##                                                                                                                                                                                                                                                         14 
    ## TGAGGAATATTGGTCAATGGGCGGGAGCCTGAACCAGCCAAGCCGCGTGAGGGAATAAGGCCCTATGGGTCGTAAACCTCTTTTGTCGGGGAACAAAAGCGGGGACGCGTCCCCGTCTGCGTGTACCCGAAGAAAAAGCATCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGTAGGCGGTATGGTAAGTCAGCGGTAAAAGCC 
    ##                                                                                                                                                                                                                                                          8 
    ## GACTACTGGGGTATCTAAGCCTGTTCGATACCCACACTTTCGTGCATGAGCGTCAGTTGAGCGCCGGTATGCTGCCTTCGCAATCGGAGTTCTGCGTGATATCTATGCATTTCACCGCTACACCACGCATTCCGCATACTTCTCGCTCACTCAAGAAAACCAGTTTCAACGGCTCGAAGAGGTTGAGCCTCTCAATTTTACCGCTGACTTGATCTTCCGCCTGCGCACCCTTTAAACCCAATAAATCCGG 
    ##                                                                                                                                                                                                                                                        117 
    ## CCTACGGGGGGCTGCAGTGGGGGATATTGGACAATGGGGGAAACCCTGATCCAGCGACGCCGCGTGAGTGAAGAAGTATTTCGGTATGTAAAGCTCTGTCAGCAGGGAAGAAAGAAATGACGGTACCTGACCAAGAAGCCCCGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGGGGCAAGCGTTATCCGGAATTACTGGGTGTAAAGGGAGCGTAGACGGTGATGCAAGTCTGGAGTGAAAGGC 
    ##                                                                                                                                                                                                                                                         41 
    ## GACTACAGGGGTATCTAATCCTGTTTGCTCCCCACACTTTCGTGCCTCAGCGTCAGTTACAGTCCAGAAAGTCGCCTTCGCCACTGATGTTCCTCCTAATATCTACGCATTCCACCGCTACACTAGGAATTCCACTTTCCTCTACTGCACTCAAGTATAGCAGTTTTAGAAGCAGCACCTGAGTTAAGCTCTAGGTATTTCACTTCTAACTTGCTACACCGCCTACGCACCCTTTACGCCCAGTCATTCC 
    ##                                                                                                                                                                                                                                                         26 
    ## TGGGGGATATTGCACAATGGGGGAAACCCTGATGCAGCGACGCCGCGTGAGTGAAGAAGTATTTCGGTACGTAAAGCTCTATCAGCAGGGAAGAATAATGACGGTACCTGACTAAGAAGCCCCGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGGGGCAAGCGTTATCCGGATTTACTGGGTGTAAAGGGAGCGTAGACGGCTGTGCAAGTCTGAAGTGAAATGCCGGGGCCCAACCCCGGAAC 
    ##                                                                                                                                                                                                                                                         26 
    ## CCTACGGGCGGCTGCAGTGGGGAATATTGCACAATGGGGGAAACCCTGATGCAGCGACGCCGCGTGAGTGAAGAAGTATTTCGGTATGTAAAGCTCTATCAGCAGGAACGAGACAAGACGGTACCTGACTAAGAAGCCCCGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGGGGCAAGCGTTATCCGGAATCACTGGGTGTAAAGGGAGCGTAGACGGCTGAGCAAGCCTGAAGTGAAAGGCGG 
    ##                                                                                                                                                                                                                                                         24 
    ## CCTACGGGCGGCTGCAGTGGGGAATATTGGGCAATGGAGGAAACTCTGACCCAGCAACGCCGCGTGAAGGATGAAGGTTTTCGGATTGTAAACTTCTGTTGCGAGGGATGAAGAATGACGGTACCTCGTGAGGAAGCCCCGGCAAACTACGTGCCAGCAGCCGCGGTAATACGTAGGGGGCAAGCGTTGTCCGGAATGACTGGGCGTAAAGGGAGTGTAGGCGGCCTTTTAAGTTATATGTGAAAGCCCC 
    ##                                                                                                                                                                                                                                                          6 
    ## CCTACGGGGGGCTGCAGTGAGGAATATTGGTCAATGGGCGCGAGCCTGAACCAGCCAAGTCGCGTGAGGGAGGACGGCCCTACGGGTTGTAAACCTCTTTTGCCGGGGAGCAACGGGCGCCACGTGTGGCGCCACTGAGAGTACCCGGAGAAAAAGCATCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGTAGGCGGATCGTT 
    ##                                                                                                                                                                                                                                                         11 
    ## GACTACCCGGGTATCTAATCCTGTTCGATACCCACACTTTCGTGCATGAGCGTCAGTTGAGCGCCGGTATGCTGCCTTCGCAATCGGAGTTCTGCGTGATATCTATGCATTTCACCGCTACACCACGCATTCCGCATACTTCTCGCTCACTCAAGAAAACCAGTTTCAACGGCTCGAAGAGGTTGAGCCTCTCAATTTTACCGCTGACTTGATCTTCCGCCTGCGCACCCTTTAAACCCAATAAATCCGG 
    ##                                                                                                                                                                                                                                                        103 
    ## GACTACCCGGGTATCTAATCCTGTTCGATACCCACGCTTTCGTGCCTCAGCGTCAGTCCGGCGCCGGTACGCTGCCTTCGCAATCGGAGTTCTGCGCGATATCTATGCATTTCACCGCTACACCGCGCATTCCGCGTACTTCTCGCCGTCTCTAGTCCGACAGTTTCAACGGCTCGCCGGGGTTGAGCCCCGGGCTTTTACCGCTGACTTTCCGGACCGCCTACGCACCCTTTAAACCCAATAAATCCGG 
    ##                                                                                                                                                                                                                                                         57 
    ## GACTACTGGGGTATCTAATCCTGTTCGATCCCCACGCTTTCGTGCCTCAGCGTCAGTTGAGCGCCGGTATGCTGCCTTCGCAATCGGAGTTCTGCGTGATATCTATGCATTTCACCGCTACACCACGCATTCCGCATACTTCTCGCTCACTCAAGACCCGCAGTTTCAACGGCGATACGGCGTTGAGCACCGCATTTTTACCGCTGACTTACAAATCCGCCTACGCACCCTTTAAACCCAATAAATCCGG 
    ##                                                                                                                                                                                                                                                         94 
    ## TGAGGAATATTGGTCAATGGGCGGAAGCCTGAACCAGCCAAGTCGCGTGAGGGATTAAGGCCCTACGGGTCGTAAACCTCTTTTGTCAGGGAGCAAATGCGCCCACGTGTGGGCGAGTGGAGAGTACCTGAAGAAAAAGCATCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGTAGGCGGGTTGTTAAGTCAGCGGTAAAATT 
    ##                                                                                                                                                                                                                                                          9 
    ## GACTACTAGGGTATCTAATCCTGTTTGATCCCCGCACTTTCGTGCCTCAGCGTCAGTCGGGCGCCGGTATGCTGCCTTCGCAATCGGGGTTCTGCGTGATATCTATGCATTTCACCGCTACACCACGCATTCCGCATACTTCTCGCCCACTCGAGCCCGGCAGTTTCAACGGCTGTACGGGGTTGAGCCCCGCAATTTTACCGCTGACTTGGCAGGCCGCCTACGCACCCTTTAAACCCAATAAATCCGG 
    ##                                                                                                                                                                                                                                                         63 
    ## CCTACGGGGGGCTGCAGTGGGGGATATTGCACAATGGGGGAAACCCTGATGCAGCGACGCCGCGTGAGGGAAGACGGTTTTCGGATTGTAAACCTCTGTCTTTGGGGACGAATTGAGACGGTACCCAAGGAGGAAGCTCCGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGGAGCAAGCGTTGTCCGGAATTACTGGGTGTAAAGGGAGCGTAGGCGGGATAGTAAGTTGATTGTTAAAACTAT 
    ##                                                                                                                                                                                                                                                          5 
    ## GACTACTAGGGTATCTAATCCTGTTTGATACCCACACTTTCGAGCCTCAGTGTCAGTTGCAGTCCAGTGAGCTGCCTTCGCAATCGGAGTTCTTCGTGATATCTAAGCATTTCACCGCTACACCACGAATTCCGCCCACCTCTACTGTACTCAAGACTGCCAGTTTCAACTGCAATTTTACGGTTGAGCCGCAAACTTTCACAACTGACTTAACAATCCACCTACGCTCCCTTTAAACCCAATAAATCCG 
    ##                                                                                                                                                                                                                                                        117 
    ## GACTACTCGGGTATCTAATCCTGTTTGATACCCACACTTTCGAGCCTCAGTGTCAGTTGCAGTCCAGTGAGCTGCCTTCGCAATCGGAGTTCTTCGTGATATCTAAGCATTTCACCGCTACACCACGAATTCCGCCCACCTCTACTGTACTCAAGACTGCCAGTTTCAACTGCAATTTTACGGTTGAGCCGCAAACTTTCACAACTGACTTAACAATCCACCTACGCTCCCTTTAAACCCAATAAATCCG 
    ##                                                                                                                                                                                                                                                        104 
    ## CCTACGGGCGGCTGCAGTGGGGAATATTGCACAATGGGGGAAACCCTGATGCAGCGACGCCGCGTGAGTGAAGAAGTATTTCGGTATGTAAAGCTCTATCAGCAGGGAAGAGAATGACGGTACCTGAGTAAGAAGCCCCGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGGGGCAAGCGTTATCCGGATTTACTGGGTGTAAAGGGAGCGTAGACGGCGAAGAAAGTCTGAAGTGAAAGCCCGC 
    ##                                                                                                                                                                                                                                                         33 
    ## TGAGGAATATTGGTCAATGGCCGGGAGGCTGAACCAGCCAAGTCGCGTGAGGGAAGACGGCCCTACGGGTTGTAAACCTCTTTTGTCGGGGAGCAAAGAGCGGGACGTGTCCCGCCTCGAGAGTACCCGAAGAAAAAGCATCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGCAGGCGGAAGATCAAGTCAGCGGTAAAATCG 
    ##                                                                                                                                                                                                                                                         43 
    ## TGAGGAATATTGGTCAATGGACGAGAGTCTGAACCAGCCAAGTAGCGTGAAGGATGACTGCCCTATGGGTTGTAAACTTCTTTTATATGGGAATAAAACGTTCCACGTGTGGGATTTTGTATGTACCATATGAATAAGGATCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGATCCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGAGCGTAGGTGGATTGTTAAGTCAGTTGTGAAAGTT 
    ##                                                                                                                                                                                                                                                        105 
    ## GACTACACGGGTATCTAATCCTGTTCGATACCCACACTTTCGTGCATGAGCGTCAGTTGAGCGCCGGTATGCTGCCTTCGCAATCGGAGTTCTGCGTGATATCTATGCATTTCACCGCTACACCACGCATTCCGCATACTTCTCGCTCACTCAAGAAAACCAGTTTCAACGGCTCGAAGAGGTTGAGCCTCTCAATTTTACCGCTGACTTGATCTTCCGCCTGCGCACCCTTTAAACCCAATAAATCCGG 
    ##                                                                                                                                                                                                                                                         92 
    ## GACTACAGGGGTATCTAAGCCTGTTCGATACCCACACTTTCGTGCATGAGCGTCAGTTGAGCGCCGGTATGCTGCCTTCGCAATCGGAGTTCTGCGTGATATCTATGCATTTCACCGCTACACCACGCATTCCGCATACTTCTCGCTCACTCAAGAAAACCAGTTTCAACGGCTCGAAGAGGTTGAGCCTCTCAATTTTACCGCTGACTTGATCTTCCGCCTGCGCACCCTTTAAACCCAATAAATCCGG 
    ##                                                                                                                                                                                                                                                         89 
    ## GACTACACGGGTATCTAATCCTGTTTGCTCCCCACGCTTTCGAGCCTCAACGTCAGTTACAGTCCAGTAAGCCGCCTTCGCCACCGGTGTTCCTCCTGATATCTACGCATTTCACCGCTACACCAGGAATTCCGCTTACCCCTCCTGCACTCAAGCCATACAGTTTCCAGAGCAGTTCGGCAGTTGAGCTGCCGCATTCCACTCCGGACTTGCATTGCCGTCTACGCTCCCTTTACACCCAGTAATTCCG 
    ##                                                                                                                                                                                                                                                         27 
    ## GACTACACGGGTATCTAAGCCTGTTTGCTCCCCACGCTTTCGCGCCTCAGCGTCAGTTGTCGTCCAGCAGGCCGCCTTCGCCACTGGTGTTCCTCCTAATATCTACGCATTTCACCGCTACACTAGGAATTCCGCCTGCCTCTCCGATACTCAAGGGCTACAGTTTCAAATGCAGTTCCGGGGTTGAGCCCCGGGATTTCACATCTGACTTGCAGCTCCGCCTACACGCCCTTTACACCCAGTAAATCCG 
    ##                                                                                                                                                                                                                                                         31 
    ## TGAGGAATATTGGTCAATGGGCGCGAGCCTGAACCAGCCAAGTCGCGTGAGGGAAGACGGTCCTATGGATTGTAAACCTCTTTTGTCGGGGAGCAAAAGACGCCACGCGTGGCGTTCCGAGAGTACCCGAAGAAAAAGCATCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGTAGGCGGCCTGCCAAGTCAGCGGTAAAATTG 
    ##                                                                                                                                                                                                                                                         32 
    ## CCTACGGGGGGCTGCAGTGAGGAATATTGGTCAATGGGCGCGAGCCTGAACCAGCCAAGTCGCGTGAGGGAAGACGGTCCTATGGATTGTAAACCTCTTTTGTCGGGGAGCAAAAGACGCCACGCGTGGCGTTCCGAGAGTACCCGAAGAAAAAGCATCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGTAGGCGGCCTGCCA 
    ##                                                                                                                                                                                                                                                         37 
    ## CCTACGGGCGGCTGCAGTGAGGAATATTGGTCAATGGGCGCGAGCCTGAACCAGCCAAGTAGCGTGCAGGACGACGGCCCTATGGGTTGTAAACTGCTTTTATACGGGGATAAAGTATGCCACGTGTGGTTTATTGCAGGTACCGTATGAATAAGGACCGGCTAATTCCGTGCCAGCAGCCGCGGTAATACGGAAGGTCCGGGCGTTATCCGGATTTATTGGGTTTAAAGGGAGCGTAGGCCGGGTTTTA 
    ##                                                                                                                                                                                                                                                        106 
    ## GACTACAGGGGTATCTAAGCCTGTTTGCTCCCCACGCTTTCGAGCCTCAGCGTCAGTTACCGTCCAGCAAGCCGCCTCCGCCACTGGTGTTCCTCCTAATATCTACGCATTTCACCGCTACACTAGGAATTCCACTTGCCTCTCCGGCACTCCAGCATGGCAGTTTCAAATGCAGTCCCGGGGTTGGGCCCCGGGTTTTCACATCTGACTTGCCACGCCGCCTGCGCTCCCTTTACACCCAGTAAATCCG 
    ##                                                                                                                                                                                                                                                         42 
    ## TGGGGAATATTGGGCAATGGGCGCAAGCCTGACCCAGCAACGCCGCGTGAAGGAAGAAGGCTTTCGGGTTGTAAACTTCTTTTCTGGGGGACGAACAAATGACGGTACCCCAGGAATAAGCCACGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGTGGCAAGCGTTATCCGGATTTATTGGGTGTAAAGGGCGTGTAGGCGGGAATGCAAGTCAGATGTGAAAACTATGGGCTCAACCCATAGC 
    ##                                                                                                                                                                                                                                                         10 
    ## TGGGGAATATTGCACAATGGGGGAAACCCTGATGCAGCGACGCCGCGTGAGTGAAGAAGTATTTCGGTATGTAAAGCTCTATCAGCAGGAACGAGACAAGACGGTACCTGACTAAGAAGCCCCGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGGGGCAAGCGTTATCCGGAATCACTGGGTGTAAAGGGAGCGTAGACGGCTGAGCAAGCCTGAAGTGAAAGGCGGGGGCCCAACCCCCGGAC 
    ##                                                                                                                                                                                                                                                          9 
    ## GACTACCAGGGTATCTAAGCCTGTTTGCTCCCCACGCTTTCGTGCCTCAACGTCAGAAGTAGCTTGGTAAGCTGCCTTCGCAATCGGTGTTCTGTATGATCTCTAAGCATTTCACCGCTACACCATACATTCCGCCTACCGCAACTACTCTCTAGCCGAACAGTATCAGAGGCAATTCCGAAGTTGAGCCTCGGGATTTCACCTCTAACTTATCCGACCGCCTACGCACCCTTTAAACCCAATAAATCCG 
    ##                                                                                                                                                                                                                                                        109 
    ## GACTACACGGGTATCTAATCCTGTTTGATCCCCGCACTTTCGTGCCTCAGCGTCAGTAGGGCGCCGGTATGCTGCCTTCGCAATCGGGGTTCTGCGTGATATCTATGCATTTCACCGCTACACCACGCATTCCGCATACTTCTCGCCCACTCGAGCCCGGCAGTTTCAACGGCTGTACGGGGTTGAGCCCCGCAATTTTACCGCTGACTTGGCAGGCCGCCTACGCACCCTTTAAACCCAATAAATCCGG 
    ##                                                                                                                                                                                                                                                         64 
    ## TGGGGAATATTGGGCAATGGGCGCAAGCCTGACCCAGCAACGCCGCGTGAAGGAAGAAGGCTTTCGGGTTGTAAACTTCTTTTGACAGGGAAGAGTAGAAGACGGTACCTGTCGAATAAGCCACGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGTGGCAAGCGTTGTCCGGATTTACTGGGTGTAAAGGGCGTGTAGCCGGGTTGACAAGTCAGATGTGAAATCCTGCGGCTTAACCGCAGAA 
    ##                                                                                                                                                                                                                                                         12 
    ## GACTACAAGGGTATCTAAGCCTGTTTGATACCCACACTTTCGAGCCTCAGTGTCAGTTGCAGTCCAGTGAGCTGCCTTCGCAATCGGAGTTCTTCGTGATATCTAAGCATTTCACCGCTACACCACGAATTCCGCCCACCTCTACTGTACTCAAGACTGCCAGTTTCAACTGCAATTTTACGGTTGAGCCGCAAACTTTCACAACTGACTTAACAATCCACCTACGCTCCCTTTAAACCCAATAAATCCG 
    ##                                                                                                                                                                                                                                                         41 
    ## GACTACTCGGGTATCTAATCCTGTTTGCTCCCCACGCTTTCGCGCCTCAGCGTCAGTTAATGTCCAGCAGGCCGCCTTCGCCACTGGTGTTCCTCCGTATATCTACGCATTTCACCGCTACACACGGAATTCCGCCTGCCTCTCCATCACTCAAGAACTACAGTTTCAAATGCAGGTTATGGGTTGAGCCCATAATTTTCACATCTGACTTGCAGTCCCGCCTACACGCCCTTTACACCCAGTAAATCCG 
    ##                                                                                                                                                                                                                                                         13 
    ## CCTACGGGGGGCTGCAGTCGGGAATATTGCGCAATGGAGGCAACTCTGACGCAGTGACGCCGCGTATAGGAAGAAGGTTTTCGGATTGTAAACTATTGTCCTCAGGGAAGAAAAAAGACTGTACCTGAGAAGAAAGCTCCGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGGAGCAAGCGTTATCCGGAATTATTGGGTGTAAAGGGTGCGTAGACGGGATAACAAGTTGGTTGTGAAACCCCT 
    ##                                                                                                                                                                                                                                                         11 
    ## GACTACCGGGGTATCTAATCCTGTTTGCTCCCCACGCCTTCGGGCCTCAACGTCAGTCACAGTCCAGCAGGCCGCCTTCGCCACCGGTGTTCCTCCTAATATCTACGCATTTCACCGCTACACTAGGAATTCCGCCTGCCCCTCCTGTACTCCAGCCATGCAGTTCCAAGAGCAGTCCCGGGGTTGGGCCCCGGGCTTTCACTCCTGGCTTGCATTGCCGTCTGCGCCCCCTTTACACCCAGTAAATCCG 
    ##                                                                                                                                                                                                                                                         13 
    ## GACTACTGGGGTATCTAAGCCTGTTCGCTCCCCACGCTTTCGTGCCTCAGCGTCAGTAATAGCTTGGAAAGCTGCCTTCGCAATCGGTGTTCTGTGTGATATCTAAGCATTTCACCGCTACACCACACATTCCGCCTTCCGCAACTACACTCTAGTCAGACAGTATCAGAGGCAGTTTCGATGTTGAGCATCGAGCTTTCACCTCTAACTTACCTGACCGCCTACGCACCCTTTAAACCCAATAAATCCG 
    ##                                                                                                                                                                                                                                                         11 
    ## GACTACCGGGGTATCTAAGCCTGTTCGATACCCACACTTTCGTGCATGAGCGTCAGTTGAGCGCCGGTATGCTGCCTTCGCAATCGGAGTTCTGCGTGATATCTATGCATTTCACCGCTACACCACGCATTCCGCATACTTCTCGCTCACTCAAGAAAACCAGTTTCAACGGCTCGAAGAGGTTGAGCCTCTCAATTTTACCGCTGACTTGATCTTCCGCCTGCGCACCCTTTAAACCCAATAAATCCGG 
    ##                                                                                                                                                                                                                                                         76 
    ## GACTACAAGGGTATCTAATCCTGTTTGATACCCACACTTTCGAGCCTCAGTGTCAGTTGCAGTCCAGTGAGCTGCCTTCGCAATCGGAGTTCTTCGTGATATCTAAGCATTTCACCGCTACACCACGAATTCCGCCCACCTCTACTGTACTCAAGACTGCCAGTTTCAACTGCAATTTTACGGTTGAGCCGCAAACTTTCACAACTGACTTAACAATCCACCTACGCTCCCTTTAAACCCAATAAATCCG 
    ##                                                                                                                                                                                                                                                         83 
    ## GACTACAAGGGTATCTAATCCTGTTTGCTCCCCACGCTTTCGTGCCTCAACGTCAGAAGTAGCTTGGTAAGCTGCCTTCGCAATCGGTGTTCTGTATGATCTCTAAGCATTTCACCGCTACACCATACATTCCGCCTACCGCAACTACTCTCTAGCCGAACAGTATCAGAGGCAATTCCGAAGTTGAGCCTCGGGATTTCACCTCTAACTTATCCGACCGCCTACGCACCCTTTAAACCCAATAAATCCG 
    ##                                                                                                                                                                                                                                                         45 
    ## GACTACCCGGGTATCTAATCCTGTTTGATACCCACACTTTCGAGCCTCAGTGTCAGTTGCAGTCCAGTGAGCTGCCTTCGCAATCGGAGTTCTTCGTGATATCTAAGCATTTCACCGCTACACCACGAATTCCGCCCACCTCTACTGTACTCAAGACTGCCAGTTTCAACTGCAATTTTACGGTTGAGCCGCAAACTTTCACAACTGACTTAACAATCCACCTACGCTCCCTTTAAACCCAATAAATCCG 
    ##                                                                                                                                                                                                                                                        118 
    ## GACTACCGGGGTATCTAATCCTGTTTGATCCCCGCACTTTCGTGCCTCAGCGTCAGTCGGGCGCCGGTATGCTGCCTTCGCAATCGGGGTTCTGCGTGATATCTATGCATTTCACCGCTACACCACGCATTCCGCATACTTCTCGCCCACTCGAGCCCGGCAGTTTCAACGGCTGTACGGGGTTGAGCCCCGCAATTTTACCGCTGACTTGGCAGGCCGCCTACGCACCCTTTAAACCCAATAAATCCGG 
    ##                                                                                                                                                                                                                                                        105 
    ## TGAGGAATATTGGTCAATGGCCGGAAGGCTGAACCAGCCAAGTCGCGTGAGGGAATAAGGCCCTACGGGTCGTAAACCTCTTTTGCCAGGGAGCAATGCCGTCCACGTGTGGACGGAAGGAGAGTACCTGGAGAAAAAGCATCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGCAGGCGGAAGATCAAGTCAGCGGTAAAATT 
    ##                                                                                                                                                                                                                                                         18 
    ## GACTACCGGGGTATCTAAGCCTGTTTGCTCCCCACGCTTTCGAGCCTCAACGTCAGTTACAGTCCAGTGAGCCGCCTTCGCCACCGGTGTTCCTCCTAATATCTACGCATTTCACCGCTACACTAGGAATTCCGCTCACCCCTCCTGCACTCTAGCTGTACAGTTTCCAAAGCAGTTCCGGGGTTGGGCCCCGGCATTTCACTTCAGACTTGCACAGCCGTCTACGCTCCCTTTACACCCAGTAAATCCG 
    ##                                                                                                                                                                                                                                                         33 
    ## TGGGGAATATTGCACAATGGGGGGAACCCTGATGCAGCGACGCCGCGTGAGTGAAGAAGTATTTCGGTATGTAAAGCTCTATCAGCAGGGAAGAAGATGACGGTACCTGACTAAGAAGCCCCGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGGGGCAAGCGTTATCCGGATTTACTGGGTGTAAAGGGAGCGTAGACGGCAGGGCAAGTCTGGAGTGAAAGCCCGGGGCCCAACCCCGGGACT 
    ##                                                                                                                                                                                                                                                         26 
    ## GACTACCAGGGTATCTAATCCCTTTCGCTCCCCTGGCCTTCGTGCCTCAGCGTCAGTTAATGTCCAGGAACCCGCCTTCGCCACGAGTGTTCCTCTCGATATCTACGCATTTCACTGCTACACCGAGAATTCCGGTTCCCCCTCCATTACTCTAGTCTCGCAGTATCATGTGCCGTCCGCGGGTTGAGCCCGCGCCTTTCACACACGACTTACGAAACAGCCTACGCACGCTTTACGCCCAGTGATTCCG 
    ##                                                                                                                                                                                                                                                         91 
    ## GACTACCCGGGTATCTAATCCTGTTTGCTCCCCACGCTTTCGTGCCTCAACGTCAGAAGTAGCTTGGTAAGCTGCCTTCGCAATCGGTGTTCTGTATGATCTCTAAGCATTTCACCGCTACACCATACATTCCGCCTACCGCAACTACTCTCTAGCCGAACAGTATCAGAGGCAATTCCGAAGTTGAGCCTCGGGATTTCACCTCTAACTTATCCGACCGCCTACGCACCCTTTAAACCCAATAAATCCG 
    ##                                                                                                                                                                                                                                                         44 
    ## GACTACTCGGGTATCTAAGCCTGTTCGATACCCGCACTTTCGTGCCTCAGCGTCAGTAGGGCGCCGGTATGCTGCCTTCGCAATCGGGGTTCTGCGTGATATCTATGCATTTCACCGCTACACCACGCATTCCGCATACTTCTCGCCCACTCAAGGCACCCAGTTCCGACGGCAGGCCGGGGTTGAGCCCCGACATTTGACCGCCGGCTTAAATGCCCGCCTACGCACCCTTTAAACCCAATAAATCCGG 
    ##                                                                                                                                                                                                                                                         27 
    ## GACTACCCGGGTATCTAATCCTGTTCGATACCCGCGCTTTCGTGCCTCAGCGTCAGTGACCGGCCGGTACGCTGCCTTCGCAATCGGGGTTCTGCGTGATATCTATGCATTTCACCGCTACACCACGCATTCCGCGTACTTCTCCGGTACTCAAGGCCCCCGGTTCCGACGGCTCGGCGGGGTTGAGCCCCGCAATTTTACCGCCGGCCTGAAGGCCCGCCTACGCACCCTTTAAACCCAATAAATCCGG 
    ##                                                                                                                                                                                                                                                         30 
    ## GACTACTCGGGTATCTAATCCCTTTCGCTCCCCTGGCCTTCGTGCCTCAGCGTCAGTTAATGTCCAGGAACCCGCCTTCGCCACGAGTGTTCCTCTCGATATCTACGCATTTCACTGCTACACCGAGAATTCCGGTTCCCCCTCCATTACTCTAGTCTCGCAGTATCATGTGCCGTCCGCGGGTTGAGCCCGCGCCTTTCACACACGACTTACGAAACAGCCTACGCACGCTTTACGCCCAGTGATTCCG 
    ##                                                                                                                                                                                                                                                        100 
    ## CCTACGGGTGGCTGCAGTGAGGAATATTGGTCAATGGCCGGAAGGCTGAACCAGCCAAGTCGCGTGAGGGAATAAGGCCCTACGGGTCGTAAACCTCTTTTGCCAGGGAGCAATGCCGTCCACGTGTGGACGGAAGGAGAGTACCTGGAGAAAAAGCATCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGTAGGCGGCCTGCC 
    ##                                                                                                                                                                                                                                                         72 
    ## GACTACTGGGGTATCTAATCCTGTTTGATCCCCGCACTTTCGTGCCTCAGCGTCAGTCGGGCGCCGGTATGCTGCCTTCGCAATCGGGGTTCTGCGTGATATCTATGCATTTCACCGCTACACCACGCATTCCGCATACTTCTCGCCCACTCGAGCCCGGCAGTTTCAACGGCTGTACGGGGTTGAGCCCCGCAATTTTACCGCTGACTTGGCAGGCCGCCTACGCACCCTTTAAACCCAATAAATCCGG 
    ##                                                                                                                                                                                                                                                         88 
    ## GACTACAGGGGTATCTAATCCTGTTTGATCCCCGCACTTTCGTGCCTCAGCGTCAGTAGGGCGCCGGTATGCTGCCTTCGCAATCGGGGTTCTGCGTGATATCTATGCATTTCACCGCTACACCACGCATTCCGCATACTTCTCGCCCACTCGAGCCCGGCAGTTTCAACGGCTGTACGGGGTTGAGCCCCGCAATTTTACCGCTGACTTGGCAGGCCGCCTACGCACCCTTTAAACCCAATAAATCCGG 
    ##                                                                                                                                                                                                                                                         99 
    ## TGAGGAATATTGGTCAATGGGCGGAAGCCTGAACCAGCCAAGTCGCGTGAGGGACTAAGGCCCTACGGGTCGTAAACCTCTTTTGCCGGGGATCAGTGCCCAGCTCGCGAGCTGGGAGGGAGCGTACCCGGAGAAAAAGCATCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGTAGGCGGGCCTTCAGGCCGGCGGTAAAATT 
    ##                                                                                                                                                                                                                                                         11 
    ## TGAGGAATATTGGTCAATGGGCGCGAGCCTGAACCAGCCAAGTCGCGTGAGGGATGACGGCCCTACGGGTTGTAAACCTCTTTTGTCGGGGAGCAAATTCCGCCACGTGTGGCGGAGTCGAGAGTACCCGAAGAAAAAGCATCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGTAGGCGGACTGTCAAGTCAGCGGTAAAATA 
    ##                                                                                                                                                                                                                                                         10 
    ## TGAGGAATATTGGTCAATGGACGCAAGTCTGAACCAGCCATGCCGCGTGCAGGAAGAAGGTGCTATGCATTGTAAACTGCTTTTATACGGGGGTAATTACAGATACGTGTATCTGGATGAAAGTACCGTATGAATAAGGATCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGATCCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGTAGGCGGTCAGGTAAGTTAGAGGTGAAAGCT 
    ##                                                                                                                                                                                                                                                          4 
    ## TGGGGAATATTGGGCAATGGGGGAAACCCTGACCCAGCAACGCCGCGTGAAGGAAGAAGGCCTTCGGGTTGTAAACTTCTTTTACCAGGGACGAAGAACGTGACGGTACCTGGAGAAAAAGCCACGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGTGGCAAGCGTTGTCCGGATTTACTGGGTGTAAAGGGCGTGTAGGCGGAGCTGCAAGTCAGATGTGAAATCCCGGGGCTCAACCCCGGA 
    ##                                                                                                                                                                                                                                                          3 
    ## GACTACTGGGGTATCTAATCCTGTTTGCTCCCCACGCTTTCGAGCCTCAACGTCAGTTACAGTCCAGTAAGCCGCCTTCGCCACTGGTGTTCCTCCTGATATCTACGCATTTCACCGCTACACCAGGAATTCCGCTTACCTCTCCTGCACTCCAGCCATTCAGTTTCCAATGCACTCCCACGGTTGAGCCCTGGGTTTTCACATCGGACTTGCATTGCCGTCTACGCTCCCTTTACACCCAGTAAATCCG 
    ##                                                                                                                                                                                                                                                         23 
    ## CCTACGGGTGGCTGCAGTGAGGAATATTGGTCAATGGACGAGAGTCTGAACCAGCCAAGTAGCGTGAAGGATGACTGCCCTATGGGTTGTAAACTTCTTTTATATGGGAATAAAATGTTCCACGTGTGGGATTTTGTATGTACCATATGAATAAGGATCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGATCCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGAGCGTAGGTGGATTGTTA 
    ##                                                                                                                                                                                                                                                         89 
    ## TGGGGGATATTGCACAATGGGGGAAACCCTGATGCAGCGACGCCGCGTGAGTGAAGAAGTATTTCGGTATGTAAAGCTCTATCGGCAGGGAAGAAAAAAATGACAGTACCTGACTAAGAAGCCCCGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGGGGCAAGCGTTATCCGGATTTACTGGGTGTAAAGGGAGCGTAGGCGGCCATGCAAGTCAGAAGTGAAAACCCGGGGCTCAACCCTGGG 
    ##                                                                                                                                                                                                                                                         20 
    ## GACTACCGGGGTATCTAATCCTGTTCGATACCCACGCTTTCGTGCCTCAGCGTCAGACAGAAGCCGGTAGGCTGCCTTCGCAATCGGAGTTCTGCGCAATATCTAAGCATTTCACCGCTACACTGCGCATTCCGCCTACTTCTCTTCCTCTCTAGCACGGCAGTTTCAACGGCTCGCCGGGGTTGAGCCCCGGGCTTTTACCGCTGACTTACCATACCGCCTACGCACCCTTTAAACCCAATAAATCCGG 
    ##                                                                                                                                                                                                                                                         18 
    ## GACTACTGGGGTATCTAAGCCTGTTTGATACCCACACTTTCGAGCCTCAGTGTCAGTTGCAGTCCAGTGAGCTGCCTTCGCAATCGGAGTTCTTCGTGATATCTAAGCATTTCACCGCTACACCACGAATTCCGCCCACCTCTACTGTACTCAAGACTGCCAGTTTCAACTGCAATTTTACGGTTGAGCCGCAAACTTTCACAACTGACTTAACAATCCACCTACGCTCCCTTTAAACCCAATAAATCCG 
    ##                                                                                                                                                                                                                                                         78 
    ## GACTACAAGGGTATCTAAGCCTGTTCGATACCCGCGCCTTCGAGCTTCAGCGTCAGTAGCGCTGCCGTATGCTGCCTTCGCAATCGGGGTTCTTCGTGATATCTAAGCATTTCACCGCTACACCACGAATTCCGCATACGTTCCGCGCACTCAAGGACTCCAGTTCGCGCCGCAGTGTCAAGGTTGAGCCCCGACATTTCACGGCACGCTTAAAACCCGGCCTACGCTCCCTTTAAACCCAATAAATCCG 
    ##                                                                                                                                                                                                                                                         74 
    ## GACTACAAGGGTATCTAATCCTGTTCGATACCCGCGCCTTCGAGCTTCAGCGTCAGTAGCGCTGCCGTATGCTGCCTTCGCAATCGGGGTTCTTCGTGATATCTAAGCATTTCACCGCTACACCACGAATTCCGCATACGTTCCGCGCACTCAAGGACTCCAGTTCGCGCCGCAGTGTCAAGGTTGAGCCCCGACATTTCACGGCACGCTTAAAACCCGGCCTACGCTCCCTTTAAACCCAATAAATCCG 
    ##                                                                                                                                                                                                                                                         33 
    ## TGAGGAATATTGGTCAATGGCCGAGAGGCTGAACCAGCCAAGTCGCGTGAGGGAAGACAGCCCTACGGGTTGTAAACCTCTTTTGTCGGGGAGCAAAGAGCGGGACGCGTCCCGCCTCGAGAGTACCCGAAGAAAAAGCATCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGCAGGCGGAAGATCAAGTCAGCGGTAAAATCG 
    ##                                                                                                                                                                                                                                                         35 
    ## GACTACCAGGGTATCTAATCCTGTTCGATACCCGCGCCTTCGAGCTTCAGCGTCAGTAGCGCTGCCGTATGCTGCCTTCGCAATCGGGGTTCTTCGTGATATCTAAGCATTTCACCGCTACACCACGAATTCCGCATACGTTCCGCGCACTCAAGGACTCCAGTTCGCGCCGCAGTGTCAAGGTTGAGCCCCGACATTTCACGGCACGCTTAAAACCCGGCCTACGCTCCCTTTAAACCCAATAAATCCG 
    ##                                                                                                                                                                                                                                                         75 
    ## GACTACCAGGGTATCTAAGCCTGTTTGATCCCCGCACTTTCGTGCCTCAGCGTCAGTAGGGCGCCGGTATGCTGCCTTCGCAATCGGGGTTCTGCGTGATATCTATGCATTTCACCGCTACACCACGCATTCCGCATACTTCTCGCCCACTCGAGCCCGGCAGTTTCAACGGCTGTACGGGGTTGAGCCCCGCAATTTTACCGCTGACTTGGCAGGCCGCCTACGCACCCTTTAAACCCAATAAATCCGG 
    ##                                                                                                                                                                                                                                                         84 
    ## GACTACCGGGGTATCTAAGCCTGTTTGATACCCACACTTTCGAGCCTCAGTGTCAGTTGCAGTCCAGTGAGCTGCCTTCGCAATCGGAGTTCTTCGTGATATCTAAGCATTTCACCGCTACACCACGAATTCCGCCCACCTCTACTGTACTCAAGACTGCCAGTTTCAACTGCAATTTTACGGTTGAGCCGCAAACTTTCACAACTGACTTAACAATCCACCTACGCTCCCTTTAAACCCAATAAATCCG 
    ##                                                                                                                                                                                                                                                         69 
    ## GACTACAGGGGTATCTAAGCCTGTTTGCTCCCCACGCTTTCGAGCCTCAACGTCAGTTACTGTCCAGTAAGCCGCCTTCGCCACTGGTGTTCCTCCTAATATCTACGCATTTCACCGCTACACTAGGAATTCCGCTTACCTCTCCAGCACTCTAGCCGGGCAGTTTCCAAAGCAGTCCCGCAGTTGGGCCGCGGGCTTTCACTTCAGACTTGCTCTGCCGTCTACGCTCCCTTTACACCCAGTAAATCCG 
    ##                                                                                                                                                                                                                                                         34 
    ## GACTACTAGGGTATCTAATCCCTTTCGCTCCCCTGGCCTTCGTGCCTCAGCGTCAGTTAATGTCCAGGAACCCGCCTTCGCCACGAGTGTTCCTCTCGATATCTACGCATTTCACTGCTACACCGAGAATTCCGGTTCCCCCTCCATTACTCTAGTCTCGCAGTATCATGTGCCGTCCGCGGGTTGAGCCCGCGCCTTTCACACACGACTTACGAAACAGCCTACGCACGCTTTACGCCCAGTGATTCCG 
    ##                                                                                                                                                                                                                                                        109 
    ## GACTACAAGGGTATCTAAGCCTGTTTGCTCCCCACGCTTTCGCGCCTCAGCGTCAGTTAATGTCCAGCAGGCCGCCTTCGCCACTGGTGTTCCTCCGTATATCTACGCATTTCACCGCTACACACGGAATTCCGCCTGCCTCTCCATCACTCAAGAAGAACAGTTTCAAACGCAGTTCCAGGGTTGAGCCCTGGAATTTCACGTCTGACTTGCTCTCCCGCCTACACGCCCTTTACACCCAGTAAATCCG 
    ##                                                                                                                                                                                                                                                         36 
    ## GACTACCGGGGTATCTAATCCTGTTTGATCCCCGCACTTTCGTGCCTCAGCGTCAGTAGGGCGCCGGTATGCTGCCTTCGCAATCGGGGTTCTGCGTGATATCTATGCATTTCACCGCTACACCACGCATTCCGCATACTTCTCGCCCACTCGAGCCCGGCAGTTTCAACGGCTGTACGGGGTTGAGCCCCGCAATTTTACCGCTGACTTGGCAGGCCGCCTACGCACCCTTTAAACCCAATAAATCCGG 
    ##                                                                                                                                                                                                                                                        120 
    ## GACTACAGGGGTATCTAATCCTGTTCGATACCCGCGCCTTCGAGCTTCAGCGTCAGTAGCGCTGCCGTATGCTGCCTTCGCAATCGGGGTTCTTCGTGATATCTAAGCATTTCACCGCTACACCACGAATTCCGCATACGTTCCGCGCACTCAAGGACTCCAGTTCGCGCCGCAGTGTCAAGGTTGAGCCCCGACATTTCACGGCACGCTTAAAACCCGGCCTACGCTCCCTTTAAACCCAATAAATCCG 
    ##                                                                                                                                                                                                                                                         58 
    ## CCTACGGGGGGCTGCAGTGAGGAATATTGGTCAATGGACGAGAGTCTGAACCAGCCAAGTAGCGTGAAGGATGACTGCCCTATGGGTTGTAAACTTCTTTTATATGGGAATAAAACGTTCCACGTGTGGGATTTTGTATGTACCATATGAATAAGGATCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGATCCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGAGCGTAGGTGGATTGTTA 
    ##                                                                                                                                                                                                                                                         78 
    ## CCTACGGGGGGCTGCAGTGAGGAATATTGGTCAATGGACGCAAGTCTGAACCAGCCATGCCGCGTGCAGGAAGAAGGTGCTATGCATTGTAAACTGCTTTTATACGGGGGTAATTACAGATACGTGTATCTGGATGAAAGTACCGTATGAATAAGGATCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGATCCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGTAGGCGGTCAGGTA 
    ##                                                                                                                                                                                                                                                          8 
    ## CCTACGGGGGGCTGCAGTGGGGGATATTGCACAATGGGGGAAACCCTGATGCAGCGACGCCGCGTGAGTGAAGAAGTATTTCGGTACGTAAAGCTCTATCAGCAGGGAAGAATAATGACGGTACCTGACTAAGAAGCCCCGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGGGGCAAGCGTTATCCGGATTTACTGGGTGTAAAGGGAGCGTAGACGGCTGTGCAAGTCTGAAGTGAAATGCCG 
    ##                                                                                                                                                                                                                                                         40 
    ## GACTACCGGGGTATCTAATCCTGTTTGCTCCCCACGCTTTCGAGCCTCAACGTCAGTTTCAGTCCAGCAGGCCGCCTTCGCCGCCGGTGTTCCTCCTGATATCTACGCATTTCACCGCTACACCAGGAATTCCGCCTGCCCCTCCTGTACTCAAGTTACACAGTTTCCAAAGCAGTCCCGGGGTTGGGCCCCGGGCTTTCACTCCAGACTTATACAACCGTCTACGCTCCCTTTACACCCAGTAAATCCG 
    ##                                                                                                                                                                                                                                                         18 
    ## GACTACCGGGGTATCTAATCCTGTTTGCTCCCCACGCTTTCGCGCCTCAGCGTCAGTTGCTGTCCAGTTGACCGCCTTCGCCACTGGTGTTCCTCCTAATATCTACGCATTTCACCGCTACACTAGGAATTCCGTCAACCTCTCCAGTACTCAAGAACGACAGTTTCAAATGCAGTTCAGGGGTTGAGCCCCTGGATTTCACATCTGACTTGCCATCCCGCCTACGCGCCCTTTACACCCAGTAAATCCG 
    ##                                                                                                                                                                                                                                                          5 
    ## CCTACGGGCGGCTGCAGTGAGGAATATTGGTCAATGGCCGGAAGGCTGAACCAGCCAAGTCGCGTGAGGGAATAAGGCCCTACGGGTCGTAAACCTCTTTTGTCAGGGAGCAAGGCCGCCCACGTGTGGACGGAAGGAGAGTACCTGGAGAAAAAGCATCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGTAGGCGGCCTGCC 
    ##                                                                                                                                                                                                                                                         62 
    ## CCTACGGGGGGCTGCAGTGAGGAATATTGGTCAATGGCCGGGAGGCTGAACCAGCCAAGTCGCGTGAGGGAAGACGGCCCTACGGGTTGTAAACCTCTTTTGTCGGGGAGCAAAGAGCGGGACGTGTCCCGCCTCGAGAGTACCCGAAGAAAAAGCATCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGCAGGCGGAAGATCA 
    ##                                                                                                                                                                                                                                                         69 
    ## GACTACAGGGGTATCTAATCCTGTTTGCTCCCCACGCTTTCGAGCCTCAGCGTCAGTCATCGTCCAGCAAGCCGCCTTCGCCACTGGTGTTCCTCCTAATATCTACGCATTTCACCGCTACACTAGGAATTCCGCTTGCCTCTCCGACACTCCAGTCTGACAGTTTCCAATGCAGTCCCGGGGTTGAGCCCCGGCCTTTCACACCAGACTTGCCAGACCGCCTGCGCTCCCTTTACACCCAGTAAATCCG 
    ##                                                                                                                                                                                                                                                         31 
    ## CCTACGGGGGGCTGCAGTGGGGAATATTGCACAATGGGGGAAACCCTGATGCAGCGACGCCGCGTGAAGGAAGAAGTATTTCGGTATGTAAACTTCTATCAGCAGGGAAGAAAATGACGGTACCTGACTAAGAAGCCCCGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGGGGCAAGCGTTATCCGGATTTACTGGGTGTAAAGGGAGCGTAGGCGGTGCAGCAAGTCAGAAGTGAAAGCCCGG 
    ##                                                                                                                                                                                                                                                         14 
    ## GACTACAGGGGTATCTAATCCTGTTTGCTCCCCACGCTTTCGAGCCTCAACGTCAGTTACAGTCCAGTAAGCCGCCTTCGCCACTGGTGTTCCTCCTAATATCTACGCATTTCACCGCTACACTAGGAATTCCACTTACCTCTCCTGCACTCCAGCCGTACAGTTCCCAAAGCAGTCCGGGGGTTGGGCCCCCGCCTTTCACTTCAGGCTTGCTCAGCCGTCTACGCTCCCTTTACACCCAGTGATTCCG 
    ##                                                                                                                                                                                                                                                         47 
    ## CCTACGGGGGGCTGCAGTGAGGAATATTGGTCAATGGGCGGAAGCCTGAACCAGCCAAGTCGCGTGAGGGATTAAGGCCCTGCGGGTCGTAAACCTCTTTTGCCGGGGAGCAGTGGTCCGGACGCGTCCGGACCGGGAGAGTACCCGGAGAAAAAGCATCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGTAGGCGGGCATTT 
    ##                                                                                                                                                                                                                                                         13 
    ## GACTACTCGGGTATCTAAGCCTGTTCGATACCCACACTTTCGTGCATGAGCGTCAGTTGAGCGCCGGTATGCTGCCTTCGCAATCGGAGTTCTGCGTGATATCTATGCATTTCACCGCTACACCACGCATTCCGCATACTTCTCGCTCACTCAAGAAAACCAGTTTCAACGGCTCGAAGAGGTTGAGCCTCTCAATTTTACCGCTGACTTGATCTTCCGCCTGCGCACCCTTTAAACCCAATAAATCCGG 
    ##                                                                                                                                                                                                                                                         62 
    ## GACTACACGGGTATCTAATCCCTTTCGCTCCCCTGGCCTTCGTGCCTCAGCGTCAGTTAATGTCCAGGAACCCGCCTTCGCCACGAGTGTTCCTCTCGATATCTACGCATTTCACTGCTACACCGAGAATTCCGGTTCCCCCTCCATTACTCTAGTCTCGCAGTATCATGTGCCGTCCGCGGGTTGAGCCCGCGCCTTTCACACACGACTTACGAAACAGCCTACGCACGCTTTACGCCCAGTGATTCCG 
    ##                                                                                                                                                                                                                                                         73 
    ## CCTACGGGGGGCTGCAGTGGGGAATATTGCACAATGGGGGGAACCCTGATGCAGCGACGCCGCGTGAGTGAAGAAGTATTTCGGTATGTAAAGCTCTATCAGCAGGGAAGAAGATGACGGTACCTGACTAAGAAGCCCCGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGGGGCAAGCGTTATCCGGATTTACTGGGTGTAAAGGGAGCGTAGACGGCAGGGCAAGTCTGGAGTGAAAGCCCGG 
    ##                                                                                                                                                                                                                                                         17 
    ## GACTACTGGGGTATCTAAGCCCTTTCGCTCCCCTGGCCTTCGTGCCTCAGCGTCAGTTAATGTCCAGGAACCCGCCTTCGCCACGAGTGTTCCTCTCGATATCTACGCATTTCACTGCTACACCGAGAATTCCGGTTCCCCCTCCATTACTCTAGTCTCGCAGTATCATGTGCCGTCCGCGGGTTGAGCCCGCGCCTTTCACACACGACTTACGAAACAGCCTACGCACGCTTTACGCCCAGTGATTCCG 
    ##                                                                                                                                                                                                                                                         53 
    ## GACTACTAGGGTATCTAAGCCTGTTCGATACCCACACTTTCGTGCATGAGCGTCAGTTGAGCGCCGGTATGCTGCCTTCGCAATCGGAGTTCTGCGTGATATCTATGCATTTCACCGCTACACCACGCATTCCGCATACTTCTCGCTCACTCAAGAAAACCAGTTTCAACGGCTCGAAGAGGTTGAGCCTCTCAATTTTACCGCTGACTTGATCTTCCGCCTGCGCACCCTTTAAACCCAATAAATCCGG 
    ##                                                                                                                                                                                                                                                         59 
    ## GACTACCGGGGTATCTAATCCTGTTTGCTCCCCACGCTTTCGTGCCTCAACGTCAGAAGTAGCTTGGTAAGCTGCCTTCGCAATCGGTGTTCTGTATGATCTCTAAGCATTTCACCGCTACACCATACATTCCGCCTACCGCAACTACTCTCTAGCCGAACAGTATCAGAGGCAATTCCGAAGTTGAGCCTCGGGATTTCACCTCTAACTTATCCGACCGCCTACGCACCCTTTAAACCCAATAAATCCG 
    ##                                                                                                                                                                                                                                                         76 
    ## CCTACGGGTGGCTGCAGTGGGGAATATTGGGCAATGGGCGCAAGCCTGACCCAGCAACGCCGCGTGAAGGAAGAAGGCTTTCGGGTTGTAAACTTCTTTTCTGGGGGACGAACAAATGACGGTACCCCAGGAATAAGCCACGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGTGGCAAGCGTTATCCGGATTTATTGGGTGTAAAGGGCGTGTAGGCGGGAATGCAAGTCAGATGTGAAAACTA 
    ##                                                                                                                                                                                                                                                         37 
    ## CCTACGGGGGGCTGCAGTGAGGAATATTGGTCAATGGCCGGAAGGCTGAACCAGCCAAGTCGCGTGAGGGAATAAGGCCCTACGGGTCGTAAACCTCTTTTGCCAGGGAGCAATGCCGTCCACGTGTGGACGGAAGGAGAGTACCTGGAGAAAAAGCATCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGCAGGCGGAAGATC 
    ##                                                                                                                                                                                                                                                         29 
    ## CCTACGGGTGGCTGCAGTCGAGAATCATTCACAATGGGGGAAACCCTGATGGTGCGACGCCGCGTGGGGGAATGAAGGTCTTCGGATTGTAAACCCCTGTCATGTGGGAGCAAATTAAAAAGATAGTACCACAAGAGGAAGAGACGGCTAACTCTGTGCCAGCAGCCGCGGTAATACAGAGGTCTCAAGCGTTGTTCGGAATCACTGGGCGTAAAGCGTGCGTAGGCTGTTTCGTAAGTCGTGTGTGAAA 
    ##                                                                                                                                                                                                                                                         50 
    ## GACTACTCGGGTATCTAATCCTGTTTGATCCCCGCACTTTCGTGCCTCAGCGTCAGTAGGGCGCCGGTATGCTGCCTTCGCAATCGGGGTTCTGCGTGATATCTATGCATTTCACCGCTACACCACGCATTCCGCATACTTCTCGCCCACTCGAGCCCGGCAGTTTCAACGGCTGTACGGGGTTGAGCCCCGCAATTTTACCGCTGACTTGGCAGGCCGCCTACGCACCCTTTAAACCCAATAAATCCGG 
    ##                                                                                                                                                                                                                                                         64 
    ## GACTACTAGGGTATCTAATCCTGTTTGATCCCCGCACTTTCGTGCCTCAGCGTCAGTAGGGCGCCGGTATGCTGCCTTCGCAATCGGGGTTCTGCGTGATATCTATGCATTTCACCGCTACACCACGCATTCCGCATACTTCTCGCCCACTCGAGCCCGGCAGTTTCAACGGCTGTACGGGGTTGAGCCCCGCAATTTTACCGCTGACTTGGCAGGCCGCCTACGCACCCTTTAAACCCAATAAATCCGG 
    ##                                                                                                                                                                                                                                                         48 
    ## TGGGGGATATTGGACAATGGGGGAAACCCTGATCCAGCGACGCCGCGTGAGTGAAGAAGTATTTCGGTATGTAAAGCTCTATCAGCAGGGAAGAAAGAAATGACGGTACCTGAGTAAGAAGCCCCGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGGGGCAAGCGTTATCCGGAATTACTGGGTGTAAAGGGAGCGTAGACGGTAATGCAAGTCTGGAGTGAAAGGCGGGGGCCCAACCCCCGG 
    ##                                                                                                                                                                                                                                                         22 
    ## GACTACAGGGGTATCTAATCCTGTTCGATACCCACGCTTTCGTGCATGAGCGTCAGTCGGGCGCCGGTACGCTGCCTTCGCAATCGGAGTTCTGCGTGATATCTATGCATTTCACCGCTACACCACGCATTCCGCGTACTTCTCACCCACTCAAGAAAACCAGTTTCAACGGCTGAAAGAGGTTGAGCCTCTCGATTTTACCGCTGACTTGATCTTCCGCCTGCGCACCCTTTAAACCCAATAAATCCGG 
    ##                                                                                                                                                                                                                                                         83 
    ## CCTACGGGGGGCTGCAGTGAGGAATATTGGTCAATGGCCGAGAGGCTGAACCAGCCAAGTCGCGTGAGGGAAGACAGCCCTACGGGTTGTAAACCTCTTTTGTCGGGGAGCAAAGAGCGGGACGCGTCCCGCCTCGAGAGTACCCGAAGAAAAAGCATCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGCAGGCGGAAGATCA 
    ##                                                                                                                                                                                                                                                         34 
    ## GACTACAGGGGTATCTAATCCTGTTTGCTCCCCACGCTTTCGTGCCTCAACGTCAGAAGTAGCTTGGTAAGCTGCCTTCGCAATCGGTGTTCTGTATGATCTCTAAGCATTTCACCGCTACACCATACATTCCGCCTACCGCAACTACTCTCTAGCCGAACAGTATCAGAGGCAATTCCGAAGTTGAGCCTCGGGATTTCACCTCTAACTTATCCGACCGCCTACGCACCCTTTAAACCCAATAAATCCG 
    ##                                                                                                                                                                                                                                                         90 
    ## CCTACGGGCGGCTGCAGTGAGGAATATTGGTCAATGGGCGCGAGCCTGAACCAGCCAAGTCGCGTGAGGGATGACGGCCCTACGGGTTGTAAACCTCTTTTGTCGGGGAGCAAATTCCGCCACGTGTGGCGGAGTCGAGAGTACCCGAAGAAAAAGCATCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGTAGGCGGACTGTC 
    ##                                                                                                                                                                                                                                                         15 
    ## GACTACTCGGGTATCTAAGCCTGTTTGATACCCACACTTTCGAGCCTCAGTGTCAGTTGCAGTCCAGTGAGCTGCCTTCGCAATCGGAGTTCTTCGTGATATCTAAGCATTTCACCGCTACACCACGAATTCCGCCCACCTCTACTGTACTCAAGACTGCCAGTTTCAACTGCAATTTTACGGTTGAGCCGCAAACTTTCACAACTGACTTAACAATCCACCTACGCTCCCTTTAAACCCAATAAATCCG 
    ##                                                                                                                                                                                                                                                         42 
    ## CCTACGGGCGGCTGCAGTCGGGAATATTGCGCAATGGAGGAAACTCTGACGCAGTGACGCCGCGTATAGGAAGAAGGTTTTCGGATTGTAAACTATTGTCGTTAGGGAAGATAAAAGACAGTACCTAAGGAGGAAGCTCCGGCTAACTATGTGCCAGCAGCCGCGGTAATACATAGGGAGCAAGCGTTATCCGGAATTATTGGGTGTAAAGGGTGCGTAGACGGAGGAACAAGTTAGTTGTGAAATCCCT 
    ##                                                                                                                                                                                                                                                         37 
    ## GACTACAAGGGTATCTAAGCCTGTTTGCTCCCCACACTTTCGAGCCTCAGCGTCAGTAAAAGCCCAGCAGGCCGCCTTCGCCACTGGTGTTCCTCCTAATATCTACGCATTTCACCGCTACACTAGGAATTCCGCCTGCCTCTACTTCACTCAAGAACTGCAGTTTTGAGTGCGGCTACCGGTTGAGCCGGTAGATTTGACACCCAACTTGCAGTCCCGCCTACGCTCCCTTTACACCCAGTAATTCCGG 
    ##                                                                                                                                                                                                                                                          4 
    ## TGGGGGATATTGCACAATGGGGGAAACCCTGATGCAGCGATGCCGCGTGAATGAAGACGGCCTTCGGGTTGTAAAGTTCTGTCGCAGGGGACGAAAATGACGGTACCCTGCAAGAAAGCTCCGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGGAGCAAGCGTTGTCCGGAATTACTGGGTGTAAAGGGAGCGTAGGCGGGAGGATAAGTTGAATGTGAAATCTATGGGCTCAACCCATAGCTG 
    ##                                                                                                                                                                                                                                                          2 
    ## CCTACGGGGGGCTGCAGTGGGGGATATTGCACAATGGGGGAAACCCTGATGCAGCGACGCCGCGTGAGTGAAGAAGCATTTCGGTGTGTAAAGCTCTATCAGCAGGGAAGAAGAAATGACGGTACCTGACTAAGAAGCCCCGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGGGGCAAGCGTTATCCGGATTTACTGGGTGTAAAGGGAGCGTAGACGGCAATGCAAGTCCGATGTGAAAACCC 
    ##                                                                                                                                                                                                                                                         16 
    ## GACTACAGGGGTATCTAAGCCTGTTTGATACCCACACTTTCGAGCCTCAGTGTCAGTTGCAGTCCAGTGAGCTGCCTTCGCAATCGGAGTTCTTCGTGATATCTAAGCATTTCACCGCTACACCACGAATTCCGCCCACCTCTACTGTACTCAAGACTGCCAGTTTCAACTGCAATTTTACGGTTGAGCCGCAAACTTTCACAACTGACTTAACAATCCACCTACGCTCCCTTTAAACCCAATAAATCCG 
    ##                                                                                                                                                                                                                                                         53 
    ## GACTACCGGGGTATCTAAGCCTGTTCGATACCCGCGCCTTCGAGCTTCAGCGTCAGTAGCGCTGCCGTATGCTGCCTTCGCAATCGGGGTTCTTCGTGATATCTAAGCATTTCACCGCTACACCACGAATTCCGCATACGTTCCGCGCACTCAAGGACTCCAGTTCGCGCCGCAGTGTCAAGGTTGAGCCCCGACATTTCACGGCACGCTTAAAACCCGGCCTACGCTCCCTTTAAACCCAATAAATCCG 
    ##                                                                                                                                                                                                                                                         41 
    ## GACTACCCGGGTATCTAAGCCTGTTCGATACCCACACTTTCGTGCATGAGCGTCAGTTGAGCGCCGGTATGCTGCCTTCGCAATCGGAGTTCTGCGTGATATCTATGCATTTCACCGCTACACCACGCATTCCGCATACTTCTCGCTCACTCAAGAAAACCAGTTTCAACGGCTCGAAGAGGTTGAGCCTCTCAATTTTACCGCTGACTTGATCTTCCGCCTGCGCACCCTTTAAACCCAATAAATCCGG 
    ##                                                                                                                                                                                                                                                         50 
    ## GACTACTCGGGTATCTAAGCCTGTTTGATCCCCGCACTTTCGTGCCTCAGCGTCAGTCGGGCGCCGGTATGCTGCCTTCGCAATCGGGGTTCTGCGTGATATCTATGCATTTCACCGCTACACCACGCATTCCGCATACTTCTCGCCCACTCGAGCCCGGCAGTTTCAACGGCTGTACGGGGTTGAGCCCCGCAATTTTACCGCTGACTTGGCAGGCCGCCTACGCACCCTTTAAACCCAATAAATCCGG 
    ##                                                                                                                                                                                                                                                         25 
    ## GACTACCCGGGTATCTAAGCCTGTTTGCTCCCCACGCTTTCGCGCCTCAGCGTCAGTTACTGTCCAGCAATCCGCCTTCGCCACTGGTGTTCCTCCGTATATCTACGCATTTCACCGCTACACACGGAATTCCGATTGCCTCTCCAGCACTCAAGAAATACAGTTTCAAATGCAGGCTATGGGTTGAGCCCATAGTTTTCACATCTGACTTGCATTCCCGCCTACACGCCCTTTACACCCAATAAATCCG 
    ##                                                                                                                                                                                                                                                         25 
    ## GACTACTGGGGTATCTAAGCCTGTTTGCTACCCACGCTTTCGAGCCTCAGCGTCAGTTATAATCCAGTAAGCCGCCTTCGCCACTGGTGTTCCTCCTAATATCTACGCATTCCACCGCTACACTAGGAATTCCGCTTACCCCTCTTACACTCAAGTCCGTCAGTCTTGAAAGCAGTTCCGGGGTTGAGCCCCGGGATTTCACTTCCAACTTAACAGACCGCCTACGCTCCCTTTACGCCCAGTCATTCCG 
    ##                                                                                                                                                                                                                                                         10 
    ## GACTACAAGGGTATCTAAGCCTGTTTGATCCCCGCACTTTCGTGCCTCAGCGTCAGTCGGGCGCCGGTATGCTGCCTTCGCAATCGGGGTTCTGCGTGATATCTATGCATTTCACCGCTACACCACGCATTCCGCATACTTCTCGCCCACTCGAGCCCGGCAGTTTCAACGGCTGTACGGGGTTGAGCCCCGCAATTTTACCGCTGACTTGGCAGGCCGCCTACGCACCCTTTAAACCCAATAAATCCGG 
    ##                                                                                                                                                                                                                                                         26 
    ## GACTACCAGGGTATCTAAGCCTGTTCGATACCCACACTTTCGTGCATGAGCGTCAGTTGAGCGCCGGTATGCTGCCTTCGCAATCGGAGTTCTGCGTGATATCTATGCATTTCACCGCTACACCACGCATTCCGCATACTTCTCGCTCACTCAAGAAAACCAGTTTCAACGGCTCGAAGAGGTTGAGCCTCTCAATTTTACCGCTGACTTGATCTTCCGCCTGCGCACCCTTTAAACCCAATAAATCCGG 
    ##                                                                                                                                                                                                                                                         51 
    ## TAAGGGATATTGGTCAATGGGGGAAACCCTGAACCAGCAACGCCGCGTGAGGGAAGACGGTTTTCGGATTGTAAACCTCTGTCCTCTGTGAAGATGATGACGGTAGCAGAGGAGGAAGCTCCGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGGAGCGAGCGTTGTCCGGATTTACTGGGTGTAAAGGGTGCGTAGGCGGCTTGGCAAGTCAGTAGTGAAATCCATGGGCTTAACCCATGAACT 
    ##                                                                                                                                                                                                                                                          6 
    ## TGAGGAATATTGGTCAATGGCCGTAAGGCTGAACCAGCCAAGTCGCGTGAGGGAAGACGGCCCTACGGGTTGTAAACCTCTTTTGTCAGGGAGCAAGGGGCAGGTCGTGACCTGCTGTGAGAGTACCTGAAGAAAAAGCATCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGCAGGCGGAGGCACAAGTCAGCGGTAAAATCG 
    ##                                                                                                                                                                                                                                                          8 
    ## CCTACGGGCGGCTGCAGTGGGGAATATTGGGCAATGGAGGAAACTCTTACCCAGCAATGCCGCGTGAATGAAGAAGGTCTTCGGATTGTAAAGTTCTTTAATTGGGGACGAAAAAAATGACGGTACCCAAGGAATAAGCCCCGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGGGGCAAGCGTTGTTCGGAATGACTGGGCGTAAAGGGTGCGTAGGCGGTGTAGCAAGTTAGAAGTGAAATAC 
    ##                                                                                                                                                                                                                                                         10 
    ## GACTACACGGGTATCTAATCCTGTTTGCTCCCCACGCTTTCGCGCCTCAGCGTCAGTTGCTGTCCAGCAGACCGCCTTCGCCACTGGTGTTCCTCCTAATATCTACGCATTTCACCGCTACACTAGGAATTCCGTCTGCCTCTCCAGTACTCAAGAACTACAGTTTCAAATGCAGGCCACAGGTTGAGCCCGTGGTTTTCACATCTGACTTGCAGTCCCGCCTACACGCCCTTTACACCCAGTAAATCCG 
    ##                                                                                                                                                                                                                                                         30 
    ## TGAGGAATATTGGTCAATGGGCGCGAGCCTGAACCAGCCAAGTCGCGTGAGGGAAGACGGTCCTATGGATTGTAAACCTCTTTTGTCGGGGAGCAAAAGACGCCACGCGTGGCGTTCCGAGAGTACCCGAAGAAAAAGCATCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGATCCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGAGCGTAGGTGGATTGTTAAGTCAGTTGTGAAAGTT 
    ##                                                                                                                                                                                                                                                         13 
    ## TGGGGAATATTGGGCAATGGGCGCAAGCCTGACCCAGCAACGCCGCGTGAAGGAAGAAGGCTTTCGGGTTGTAAACTTCTTTTCTCAGGGACGAAGCAAGTGACGGTACCTGAGGAATAAGCCACGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGTGGCAAGCGTTATCCGGATTTACTGGGTGTAAAGGGCGTGTAGGCGGGACTGCAAGTCAGATGTGAAAACCACGGGCTCAACCTGTGG 
    ##                                                                                                                                                                                                                                                         15 
    ## GACTACACGGGTATCTAAGCCTGTTCGATACCCACACTTTCGTGCATGAGCGTCAGTTGAGCGCCGGTATGCTGCCTTCGCAATCGGAGTTCTGCGTGATATCTATGCATTTCACCGCTACACCACGCATTCCGCATACTTCTCGCTCACTCAAGAAAACCAGTTTCAACGGCTCGAAGAGGTTGAGCCTCTCAATTTTACCGCTGACTTGATCTTCCGCCTGCGCACCCTTTAAACCCAATAAATCCGG 
    ##                                                                                                                                                                                                                                                         49 
    ## GACTACTCGGGTATCTAAGCCTGTTTGCTCCCCACGCTTTCGTGCCTCAACGTCAGAAGTAGCTTGGTAAGCTGCCTTCGCAATCGGTGTTCTGTATGATCTCTAAGCATTTCACCGCTACACCATACATTCCGCCTACCGCAACTACTCTCTAGCCGAACAGTATCAGAGGCAATTCCGAAGTTGAGCCTCGGGATTTCACCTCTAACTTATCCGACCGCCTACGCACCCTTTAAACCCAATAAATCCG 
    ##                                                                                                                                                                                                                                                         21 
    ## TGAGGAATATTGGTCAATGGACGGAAGTCTGAACCAGCCATGCCGCGTGCAGGAAGACGGCTCTATGAGTTGTAAACTGCTTTTGTACAAGGGTAAACCTGAATACGTGTATTCAGCTGAAAGTACTGTACGAATAAGGATCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGATCCAAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGTAGGCGGATTGATAAGTTAGAGGTGAAATGT 
    ##                                                                                                                                                                                                                                                          9 
    ## GACTACAAGGGTATCTAAGCCTGTTCGATACCCACACTTTCGTGCATGAGCGTCAGTTGAGCGCCGGTATGCTGCCTTCGCAATCGGAGTTCTGCGTGATATCTATGCATTTCACCGCTACACCACGCATTCCGCATACTTCTCGCTCACTCAAGAAAACCAGTTTCAACGGCTCGAAGAGGTTGAGCCTCTCAATTTTACCGCTGACTTGATCTTCCGCCTGCGCACCCTTTAAACCCAATAAATCCGG 
    ##                                                                                                                                                                                                                                                         46 
    ## GACTACTCGGGTATCTAATCCTGTTCGATACCCACACTTTCGTGCATGAGCGTCAGTTGAGCGCCGGTATGCTGCCTTCGCAATCGGAGTTCTGCGTGATATCTATGCATTTCACCGCTACACCACGCATTCCGCATACTTCTCGCCCACTCGAGCCCGGCAGTTTCAACGGCTGTACGGGGTTGAGCCCCGCAATTTTACCGCTGACTTGGCAGGCCGCCTACGCACCCTTTAAACCCAATAAATCCGG 
    ##                                                                                                                                                                                                                                                         46 
    ## GACTACTCGGGTATCTAATCCTGTTTGATCCCCGCACTTTCGTGCCTCAGCGTCAGTCGGGCGCCGGTATGCTGCCTTCGCAATCGGGGTTCTGCGTGATATCTATGCATTTCACCGCTACACCACGCATTCCGCATACTTCTCGCCCACTCGAGCCCGGCAGTTTCAACGGCTGTACGGGGTTGAGCCCCGCAATTTTACCGCTGACTTGGCAGGCCGCCTACGCACCCTTTAAACCCAATAAATCCGG 
    ##                                                                                                                                                                                                                                                         43 
    ## GACTACTAGGGTATCTAAGCCTGTTTGATACCCACACTTTCGAGCCTCAGTGTCAGTTGCAGTCCAGTGAGCTGCCTTCGCAATCGGAGTTCTTCGTGATATCTAAGCATTTCACCGCTACACCACGAATTCCGCCCACCTCTACTGTACTCAAGACTGCCAGTTTCAACTGCAATTTTACGGTTGAGCCGCAAACTTTCACAACTGACTTAACAATCCACCTACGCTCCCTTTAAACCCAATAAATCCG 
    ##                                                                                                                                                                                                                                                         42 
    ## GACTACCAGGGTATCTAATCCTGTTTGCTCCCCACGCTTTCGTGCCTCAACGTCAGAAGTAGCTTGGTAAGCTGCCTTCGCAATCGGTGTTCTGTATGATCTCTAAGCATTTCACCGCTACACCATACATTCCGCCTACCGCAACTACTCTCTAGCCGAACAGTATCAGAGGCAATTCCGAAGTTGAGCCTCGGGATTTCACCTCTAACTTATCCGACCGCCTACGCACCCTTTAAACCCAATAAATCCG 
    ##                                                                                                                                                                                                                                                         39 
    ## GACTACTCGGGTATCTAATCCTGTTTGCTCCCCACGCTTTCGTGCCTCAACGTCAGAAGTAGCTTGGTAAGCTGCCTTCGCAATCGGTGTTCTGTATGATCTCTAAGCATTTCACCGCTACACCATACATTCCGCCTACCGCAACTACTCTCTAGCCGAACAGTATCAGAGGCAATTCCGAAGTTGAGCCTCGGGATTTCACCTCTAACTTATCCGACCGCCTACGCACCCTTTAAACCCAATAAATCCG 
    ##                                                                                                                                                                                                                                                         47 
    ## CCTACGGGCGGCTGCAGTGAGGAATATTGGTCAATGGGAGAGATCCTGAACCAGCCAAGCCGCGTGAGGGAAGACGGCCCTATGGGTTGTAAACCTCTTTTGTCGGAGAACAAAACCCGGGACGAGTCCCGGACTGCGTGTATCCGAAGAAAAAGCATCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGTAGGCGGTCCGTTA 
    ##                                                                                                                                                                                                                                                         27 
    ## GACTACTGGGGTATCTAATCCTGTTTGCTCCCCACGCTTTCGAGCCTCAACGTCAGTTACAGTCCAGTAAGCCGCCTTCGCCACTGGTGTTCCTCCTAATATCTACGCATTTCACCGCTACACTAGGAATTCCGCTTACCTCTCCTGCACTCAAGTTACACAGTTTCCAGAGCAGTCCGGGGGTTGGGCCCCCGCCTTTCACTCCAGACTTGCATCACCGTCTACGCTCCCTTTACACCCAGTAATTCCG 
    ##                                                                                                                                                                                                                                                         53 
    ## CCTACGGGCGGCTGCAGTGAGGAATATTGGTCAATGGACGAGAGTCTGAACCAGCCAAGTAGCGTGAAGGATGACTGCCCTATGGGTTGTAAACTTCTTTTATATGGGAATAAAACGTTCCACGTGTGGGATTTTGTATGTACCATATGAATAAGGATCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGATCCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGAGCGTAGGTGGATTGTTA 
    ##                                                                                                                                                                                                                                                         63 
    ## GACTACACGGGTATCTAATCCTGTTTGATCCCCGCACTTTCGTGCCTCAGCGTCAGTCGGGCGCCGGTATGCTGCCTTCGCAATCGGGGTTCTGCGTGATATCTATGCATTTCACCGCTACACCACGCATTCCGCATACTTCTCGCCCACTCGAGCCCGGCAGTTTCAACGGCTGTACGGGGTTGAGCCCCGCAATTTTACCGCTGACTTGGCAGGCCGCCTACGCACCCTTTAAACCCAATAAATCCGG 
    ##                                                                                                                                                                                                                                                         60 
    ## CCTACGGGGGGCTGCAGTGAGGAATATTGGTCAATGGACGCAAGTCTGAACCAGCCATGCCGCGTGCAGGAAGACGGCTCTATGAGTTGTAAACTGCTTTTGTACGAGGGTAAACGCAGATACGCGTATCTGTCTGAAAGTATCGTACGAATAAGGATCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGATTCAAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGTAGGCGGTCGGATA 
    ##                                                                                                                                                                                                                                                         59 
    ## GACTACCGGGGTATCTAATCCTGTTTGCTCCCCACGCTTTCGCGCCTCAGCGTCAGTTAATGTCCAGCAGGCCGCCTTCGCCACTGGTGTTCCTCCCAATATCTACGCATTTCACCGCTACACTGGGAATTCCGCCTGCCTCTCCATCACTCAAGACACGCAGTTCCAAAAGCAGTTTGGGGGTTAAGCCCCCAGATTTCACTTCTGGCTTACATGCCCGCCTACACGCCCTTTACACCCAGTAAATCCG 
    ##                                                                                                                                                                                                                                                         18 
    ## CCTACGGGCGGCTGCAGTGGGGAATATTGCACAATGGGGGGAACCCTGATGCAGCGACGCCGCGTGAGTGAAGAAGTATTTCGGTACGTAAAGCTCTATCAGCAGGGAAGAAAAAAGGACGGTACCTGAGTAAGAAGCCCCGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGGGGCAAGCGTTATCCGGAATTACTGGGTGTAAAGGGAGCGTAGACGGCAATGCAAGTCCGGAGTGGAATGCG 
    ##                                                                                                                                                                                                                                                         12 
    ## CCTACGGGCGGCTGCAGTGAGGAATATTGGTCAATGGCCGAGAGGCTGAACCAGCCAAGTCGCGTGAGGGAAGACGGCCCTACGGGTTGTAAACCTCTTTTGTCGGGGAGCAAAGAGCGGGACGCGTCCCGCCTCGAGAGTACCCGAAGAAAAAGCATCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGCAGGCGGAAGATCA 
    ##                                                                                                                                                                                                                                                         53 
    ## CCTACGGGTGGCTGCAGTGAGGAATATTGGTCAATGGACGCAAGTCTGAACCAGCCATGCCGCGTGCAGGAAGACGGCTCTATGAGTTGTAAACTGCTTTTGTACGAGGGTAAACGCAGATACGTGTATCTGTCTGAAAGTATCGTACGAATAAGGATCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGATTCAAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGTAGGCGGTCGGATA 
    ##                                                                                                                                                                                                                                                         54 
    ## GACTACCAGGGTATCTAATCCTGTTTGCTCCCCACGCTTTCGAGCCTCAGCGTCATTTACAGTCCAGTAAGCCGCCTTCGCCACTGGTGTTCCTCCTAATATCTACGCATTTCACCGCTACACTAGGAATTCCACTTACCTCTCCTGCAATCTAGTTCAACAGTTTCAAAAGCAGTCCCGGAGTTGAGCCCCGGGCTTTCACTTCTGACTTGCTGCACCGCCTACGCTCCCTTTACACCCAGTAAATCCG 
    ##                                                                                                                                                                                                                                                         17 
    ## TGGGGAATATTGCACAATGGGGGAAACCCTGATGCAGCGACGCCGCGTGAGCGATGAAGTATTTCGGTATGTAAAGCTCTATCAGCAGGGAAGAAAATGACGGTACCTGACTAAGAAGCCCCGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGGGGCAAGCGTTATCCGGATTTACTGGGTGTAAAGGGAGCGTAGACGGTGTGGCAAGCCAGATGTGAAAGCCCGGGGCTCAACCCCGGGACT 
    ##                                                                                                                                                                                                                                                         11 
    ## TGAGGAATATTGGTCAATGGGAGAGATCCTGAACCAGCCAAGCCGCGTGAGGGAAGACGGCCCTATGGGTTGTAAACCTCTTTTGTCGGAGAACAAAACCCGGGACGAGTCCCGGACTGCGTGTATCCGAAGAAAAAGCATCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGTAGGCGGTCCGTTAAGTCAGCGGTAAAATTG 
    ##                                                                                                                                                                                                                                                         45 
    ## TGGGGAATATTGCACAATGGGGGAAACCCTGATGCAGCAACGCCGCGTGAGTGAAGAAGTATTTCGGTATGTAAAGCTCTATCAGCAGGGAAGAAGATGACGGTACCTGAGTAAGAAGCCCCGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGGGGCAAGCGTTATCCGGATTTACTGGGTGTAAAGGGAGCGTAGACGGTGAGACAAGTCTGAAGTGAAAGCCCGGGGCTCAACCCCGGGACT 
    ##                                                                                                                                                                                                                                                         12 
    ## CCTACGGGGGGCTGCAGTGGGGAATATTGCGCAATGGGCGAAAGCCTGACGCAGCGACGCCGCGTGAGGGATGAAGGTCTTCGGATCGTAAACCTCTGTCAGCAGGGAAGAACGGTCACTGTGCTAATCAGCAGTGAATTGACGGTACCTGCAAAGGAAGCACCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGGTGCGAGCGTTAATCGGAATTACTGGGCGTAAAGCGCTCGTAGGCGGTA 
    ##                                                                                                                                                                                                                                                          4 
    ## GACTACCCGGGTATCTAATCCTGTTTGATCCCCGCACTTTCGTGCCTCAGCGTCAGTAGGGCGCCGGTATGCTGCCTTCGCAATCGGGGTTCTGCGTGATATCTATGCATTTCACCGCTACACCACGCATTCCGCATACTTCTCGCCCACTCGAGCCCGGCAGTTTCAACGGCTGTACGGGGTTGAGCCCCGCAATTTTACCGCTGACTTGGCAGGCCGCCTACGCACCCTTTAAACCCAATAAATCCGG 
    ##                                                                                                                                                                                                                                                         42 
    ## GACTACAAGGGTATCTAATCCTGTTTGATCCCCGCACTTTCGTGCCTCAGCGTCAGTCGGGCGCCGGTATGCTGCCTTCGCAATCGGGGTTCTGCGTGATATCTATGCATTTCACCGCTACACCACGCATTCCGCATACTTCTCGCCCACTCGAGCCCGGCAGTTTCAACGGCTGTACGGGGTTGAGCCCCGCAATTTTACCGCTGACTTGGCAGGCCGCCTACGCACCCTTTAAACCCAATAAATCCGG 
    ##                                                                                                                                                                                                                                                         35 
    ## GACTACTAGGGTATCTAATCCTGTTCGATACCCACGCTTTCGTGCATGAGCGTCAGTCGGGCGCCGGTACGCTGCCTTCGCAATCGGAGTTCTGCGTGATATCTATGCATTTCACCGCTACACCACGCATTCCGCGTACTTCTCACCCACTCAAGAAAACCAGTTTCAACGGCTGAAAGAGGTTGAGCCTCTCGATTTTACCGCTGACTTGATCTTCCGCCTGCGCACCCTTTAAACCCAATAAATCCGG 
    ##                                                                                                                                                                                                                                                         51 
    ## GACTACCGGGGTATCTAAGCCCTTTCGCTCCCCTGGCCTTCGTGCCTCAGCGTCAGTTAATGTCCAGGAACCCGCCTTCGCCACGAGTGTTCCTCTCGATATCTACGCATTTCACTGCTACACCGAGAATTCCGGTTCCCCCTCCATTACTCTAGTCTCGCAGTATCATGTGCCGTCCGCGGGTTGAGCCCGCGCCTTTCACACACGACTTACGAAACAGCCTACGCACGCTTTACGCCCAGTGATTCCG 
    ##                                                                                                                                                                                                                                                         49 
    ## GACTACAGGGGTATCTAATCCTGTTTGCTCCCCACACTTTCGAGCCTCAGCGTCAGTTAAAGCCCAGTTGGCCGCCTTCGCCACCGGTGTTCCTCCGAATATCTACGCATTTCACCGCTACACTCGGAATTCCGCCAACCTCTACTTCACTCAAGAAAGCCAGTTTCAACTGCAGTCTACAGGTTAAGCCCGTAGTTTTCACAGCTGACTTGGCTTCCCGCCTACGCTCCCTTTACACCCAGTAATTCCG 
    ##                                                                                                                                                                                                                                                         96 
    ## TGGGGAATATTGCACAATGGGCGAAAGCCTGATGCAGCAACGCCGCGTGAAGGAAGAAGGTCTTCGGATCGTAAACTTCTGTCCTTGGGGAAGATAATGACGGTACCCTTGGAGGAAGCCCCGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGGGGCAAGCGTTATCCGGAATTATTGGGCGTAAAGAGTGCGTAGGTGGTTACTTAAGCGCGGGGTTTAAGGCAATGGCTCAACCATTGTTCG 
    ##                                                                                                                                                                                                                                                          3 
    ## CCTACGGGGGGCTGCAGTGGGGAATATTGCACAATGGGCGAAAGCCTGATGCAGCAACGCCGCGTGAAGGAAGACGGTTTTCGGATTGTAAACTTCTATCAATAGGGAAGAAAGAAATGACGGTACCTAAATAAGAAGCCCCGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGGGGCAAGCGTTATCCGGAATTACTGGGTGTAAAGGGAGAGTAGGCGGCAAGGTAAGCGATATGTGAAAGCC 
    ##                                                                                                                                                                                                                                                         10 
    ## GACTACACGGGTATCTAATCCTGTTTGCTCCCCACGCTTTCGTGCCTCAACGTCAGAAGTAGCTTGGTAAGCTGCCTTCGCAATCGGTGTTCTGTATGATCTCTAAGCATTTCACCGCTACACCATACATTCCGCCTACCGCAACTACTCTCTAGCCGAACAGTATCAGAGGCAATTCCGAAGTTGAGCCTCGGGATTTCACCTCTAACTTATCCGACCGCCTACGCACCCTTTAAACCCAATAAATCCG 
    ##                                                                                                                                                                                                                                                         54 
    ## GACTACCGGGGTATCTAATCCTGTTTGCTCCCCACGCTTTCGCGCCTCACCGTCAGTTGCCGTCCAGTCATCCGCCTTCGCCACTGGTGTTCTTCCTTATATCTACGCATTTCACCGCTACACAAGGAATTCCGATGACCTCTCCGGTACTCAAGAAAAACAGTTTCAAATGCAGTTCCGCGGTTGAGCCGCGGGATTTCACATCTGACTTGTCTTCCCGGCTGCACGCCCTTTACACCCAGTAAATCCG 
    ##                                                                                                                                                                                                                                                         33 
    ## GACTACTGGGGTATCTAAGCCTGTTTGATCCCCGCACTTTCGTGCCTCAGCGTCAGTCGGGCGCCGGTATGCTGCCTTCGCAATCGGGGTTCTGCGTGATATCTATGCATTTCACCGCTACACCACGCATTCCGCATACTTCTCGCCCACTCGAGCCCGGCAGTTTCAACGGCTGTACGGGGTTGAGCCCCGCAATTTTACCGCTGACTTGGCAGGCCGCCTACGCACCCTTTAAACCCAATAAATCCGG 
    ##                                                                                                                                                                                                                                                         41 
    ## TGGGGAATATTGCACAATGGGGGAAACCCTGATGCAGCAACGCCGCGTGAGTGAAGAAGTATTTCGGTATGTAAAGCTCTATCAGCAGGGAAGAAGAAAGACGGTACCTGACTAAGAAGCCCCGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGGGGCAAGCGTTATCCGGATTTACTGGGTGTAAAGGGAGCGTAGACGGCTGTGCAAGTCTGAAGTGAAAGCCCGTGGCTCAACCGCGGAAC 
    ##                                                                                                                                                                                                                                                         19 
    ## TGAGGAATATTGGTCAATGGACGCAAGTCTGAACCAGCCATGCCGCGTGCAGGAAGACGGCTCTATGAGTTGTAAACTGCTTTTGTACGAGGGTAAACGCAGATACGCGTATCTGTCTGAAAGTATCGTACGAATAAGGATCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGATTCAAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGTAGGCGGTCGGATAAGTTAGAGGTGAAATCC 
    ##                                                                                                                                                                                                                                                         52 
    ## GACTACACGGGTATCTAATCCTGTTCGATACCCACGCTTTCGTGCATGAGCGTCAGTTGGGCGCCGGTACGCTGCCTTCGCAATCGGAGTTCTGCGTGATATCTATGCATTTCACCGCTACACCACGCATTCCGCGTACTTCTCACCCACTCAAGAAAACCAGTTTCAACGGCTGAAAGAGGTTGAGCCTCTCGATTTTACCGCTGACTTGATCTTCCGCCTGCGCACCCTTTAAACCCAATAAATCCGG 
    ##                                                                                                                                                                                                                                                         56 
    ## GACTACAGGGGTATCTAAGCCCTTTCGCTCCCCTGGCCTTCGTGCCTCAGCGTCAGTTAATGTCCAGGAACCCGCCTTCGCCACGAGTGTTCCTCTCGATATCTACGCATTTCACTGCTACACCGAGAATTCCGGTTCCCCCTCCATTACTCTAGTCTCGCAGTATCATGTGCCGTCCGCGGGTTGAGCCCGCGCCTTTCACACACGACTTACGAAACAGCCTACGCACGCTTTACGCCCAGTGATTCCG 
    ##                                                                                                                                                                                                                                                         27 
    ## GACTACCGGGGTATCTAATCCTGTTCGATACCCGCGCCTTCGAGCTTCAGCGTCAGTAGCGCTGCCGTATGCTGCCTTCGCAATCGGGGTTCTTCGTGATATCTAAGCATTTCACCGCTACACCACGAATTCCGCATACGTTCCGCGCACTCAAGGACTCCAGTTCGCGCCGCAGTGTCAAGGTTGAGCCCCGACATTTCACGGCACGCTTAAAACCCGGCCTACGCTCCCTTTAAACCCAATAAATCCG 
    ##                                                                                                                                                                                                                                                         39 
    ## GACTACTCGGGTATCTAAGCCTGTTCGATACCCGCGCCTTCGAGCTTCAGCGTCAGTAGCGCTGCCGTATGCTGCCTTCGCAATCGGGGTTCTTCGTGATATCTAAGCATTTCACCGCTACACCACGAATTCCGCATACGTTCCGCGCACTCAAGGACTCCAGTTCGCGCCGCAGTGTCAAGGTTGAGCCCCGACATTTCACGGCACGCTTAAAACCCGGCCTACGCTCCCTTTAAACCCAATAAATCCG 
    ##                                                                                                                                                                                                                                                         18 
    ## TGAGGAATATTGGTCAATGGGCGGAAGCCTGAACCAGCCAAGTCGCGTGAGGGAAGACGGTCCTATGGATTGTAAACCTCTTTTGCCGGGGAGCAATGCCACCTTTGCGAAGGTGGAGGGAGAGTACCCGGAGAAAAAGCATCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGTAGGCGGGCTGTTAAGTCAGCGGTCAAATG 
    ##                                                                                                                                                                                                                                                          6 
    ## GACTACACGGGTATCTAATCCTGTTTGCTCCCCACGCTTTCGAGCCTCAACGTCAGTTTCCGTCCAGTAAGCCGCCTTCGCCACTGGTGTTCCTCCTAATATCTACGCATTTCACCGCTACACTAGGAATTCCGCTTACCTCTCCGGTACTCTAGTTACATAGTTTCCAATGCAGTCCCGGGGTTGAGCCCCGGGCTTTCACATCAGACTTACATTACCGTCTACGCTCCCTTTACACCCAGTAAATCCG 
    ##                                                                                                                                                                                                                                                         42 
    ## TGAGGAATATTGGTCAATGGCCGAGAGGCTGAACCAGCCAAGTCGCGTGAGGGAAGACGGCCCTACGGGTTGTAAACCTCTTTTGTCGGGGAGCAAAGAGCGGGACGCGTCCCGCCTCGAGAGTACCCGAAGAAAAAGCATCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGCAGGCGGAAGATCAAGTCAGCGGTAAAATCG 
    ##                                                                                                                                                                                                                                                         44 
    ## TGGGGAATATTGCACAATGGGGGAAACCCTGATGCAGCGACGCCGCGTGGGTGACGAAGTATTTCGGTATGTAAAGCCCTATCAGCAGGGAAGAAGATGACAGTACCTGACTAAGAAGCCCCGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGGGGCAAGCGTTATCCGGATTTACTGGGTGTAAAGGGAGCGTAGGTGGCATGGCAAGCCAGAAGTGAAAACCCGGGGCTTAACCCCGCGGAT 
    ##                                                                                                                                                                                                                                                          3 
    ## GACTACCAGGGTATCTAATCCTGTTTGATCCCCGCACTTTCGTGCCTCAGCGTCAGTCGGGCGCCGGTATGCTGCCTTCGCAATCGGGGTTCTGCGTGATATCTATGCATTTCACCGCTACACCACGCATTCCGCATACTTCTCGCCCACTCGAGCCCGGCAGTTTCAACGGCTGTACGGGGTTGAGCCCCGCAATTTTACCGCTGACTTGGCAGGCCGCCTACGCACCCTTTAAACCCAATAAATCCGG 
    ##                                                                                                                                                                                                                                                         31 
    ## GACTACACGGGTATCTAAGCCTGTTTGATACCCACACTTTCGAGCCTCAGTGTCAGTTGCAGTCCAGTGAGCTGCCTTCGCAATCGGAGTTCTTCGTGATATCTAAGCATTTCACCGCTACACCACGAATTCCGCCCACCTCTACTGTACTCAAGACTGCCAGTTTCAACTGCAATTTTACGGTTGAGCCGCAAACTTTCACAACTGACTTAACAATCCACCTACGCTCCCTTTAAACCCAATAAATCCG 
    ##                                                                                                                                                                                                                                                         34 
    ## CCTACGGGCGGCTGCAGTCGGGAATATTGCGCAATGGAGGAAACTCTGACGCAGTGACGCCGCGTGCAGGAAGAAGGTTTTCGGATTGTAAACTGCTTTAGACAGGGAAGAAACAAGACAGTACCTGTAGAATAAGCTCCGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGGAGCAAGCGTTATCCGGATTTATTGGGTGTAAAGGGTGCGTAGACGGGAAGGCAAGTTGGTTGTGAAATCCCT 
    ##                                                                                                                                                                                                                                                          5 
    ## CCTACGGGCGGCTGCAGTGAGGAATATTGGTCAATGGACGGAAGTCTGAACCAGCCATGCCGCGTGCAGGAAGACGGCTCTATGAGTTGTAAACTGCTTTTGTACAAGGGTAAACCTGAATACGTGTATTCAGCTGAAAGTACTGTACGAATAAGGATCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGATCCAAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGTAGGCGGATTGATA 
    ##                                                                                                                                                                                                                                                         13 
    ## CCTACGGGCGGCTGCAGTGGGGAATATTGCACAATGGGCGAAAGCCTGATGCAGCAACGCCGCGTGAAGGAAGACGGTTTTCGGATTGTAAACTTCTGTTCTTAGTGAAGAATAATGACGGTAACTAAGGAGCAAGCCACGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGTGGCAAGCGTTGTCCGGAATTACTGGGTGTAAAGGGAGCGTAGGCGGGAAGCCAAGTCAGCTGTGAAAACTAC 
    ##                                                                                                                                                                                                                                                         30 
    ## GACTACTCGGGTATCTAATCCTGTTCGATACCCGCGCCTTCGAGCTTCAGCGTCAGTAGCGCTGCCGTATGCTGCCTTCGCAATCGGGGTTCTTCGTGATATCTAAGCATTTCACCGCTACACCACGAATTCCGCATACGTTCCGCGCACTCAAGGACTCCAGTTCGCGCCGCAGTGTCAAGGTTGAGCCCCGACATTTCACGGCACGCTTAAAACCCGGCCTACGCTCCCTTTAAACCCAATAAATCCG 
    ##                                                                                                                                                                                                                                                         31 
    ## TGGGGAATATTGCACAATGGGGGAAACCCTGATGCAGCGACGCCGCGTGAGTGAAGAAGCGTTTCGGCGCGTAAAGCTCTGTCAGCGGGGAAGAAACACGGTTCGCGAGAGCCGCGGACGGTACCCGACCAAGAAGCCCCGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGGGGCAAGCGTTATCCGGATTCACTGGGTGTAAAGGGAGCGTAGACGGCGAGACAAGTCTGGAGTGAAAGGCCG 
    ##                                                                                                                                                                                                                                                          2 
    ## GACTACAAGGGTATCTAATCCCTTTCGCTCCCCACGCTTTCGAGCCTCAACGTCAGTTACTGTCCAGTAAGCCGCCTTCGCCACTGGTGTTCTTCCTAATATCTACGCATTTCACCGCTACACTAGGAATTCCGCTTACCTCTCCAGCACTCCAGCTTAACAGTTTCCAAAGCAGTCCCAGGGTTGAGCCCTGGGCTTTCACTTCAGACTTGCTATGCCGTCTACGCTCCCTTTACACCCAGTAAATCCG 
    ##                                                                                                                                                                                                                                                         11 
    ## GACTACAAGGGTATCTAATCCTGTTTGATCCCCGCACTTTCGTGCCTCAGCGTCAGTAGGGCGCCGGTATGCTGCCTTCGCAATCGGGGTTCTGCGTGATATCTATGCATTTCACCGCTACACCACGCATTCCGCATACTTCTCGCCCACTCGAGCCCGGCAGTTTCAACGGCTGTACGGGGTTGAGCCCCGCAATTTTACCGCTGACTTGGCAGGCCGCCTACGCACCCTTTAAACCCAATAAATCCGG 
    ##                                                                                                                                                                                                                                                         35 
    ## TGAGGAATATTGGTCAATGGCCGGAAGGCTGAACCAGCCAAGTCGCGTGAGGGAATAAGGCCCTATGGGTCGTAAACCTCTTTTGTCAGGGAGCAAAGGCGTCCACGAGTGGACGAAAGGAGAGTACCTGAAGAAAAAGCATCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGTAGGCGGAGTGTCAAGTCAGCGGTAAAATT 
    ##                                                                                                                                                                                                                                                         17 
    ## GACTACTAGGGTATCTAATCCTGTTTGCTCCCCACGCTTTCGTGCCTCAACGTCAGAAGTAGCTTGGTAAGCTGCCTTCGCAATCGGTGTTCTGTATGATCTCTAAGCATTTCACCGCTACACCATACATTCCGCCTACCGCAACTACTCTCTAGCCGAACAGTATCAGAGGCAATTCCGAAGTTGAGCCTCGGGATTTCACCTCTAACTTATCCGACCGCCTACGCACCCTTTAAACCCAATAAATCCG 
    ##                                                                                                                                                                                                                                                         45 
    ## GACTACACGGGTATCTAATCCTGTTCGATACCCGCGCCTTCGAGCTTCAGCGTCAGTAGCGCTGCCGTATGCTGCCTTCGCAATCGGGGTTCTTCGTGATATCTAAGCATTTCACCGCTACACCACGAATTCCGCATACGTTCCGCGCACTCAAGGACTCCAGTTCGCGCCGCAGTGTCAAGGTTGAGCCCCGACATTTCACGGCACGCTTAAAACCCGGCCTACGCTCCCTTTAAACCCAATAAATCCG 
    ##                                                                                                                                                                                                                                                         23 
    ## TGAGGAATATTGGTCAATGGCCGGAAGGCTGAACCAGCCAAGTCGCGTGAGGGAATAAGGCCCTACGGGTCGTAAACCTCTTTTGTCAGGGAGCAAGGCCGCCCACGTGTGGACGGAAGGAGAGTACCTGGAGAAAAAGCATCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGCAGGCGGAAGATCAAGTCAGCGGTAAAATT 
    ##                                                                                                                                                                                                                                                         10 
    ## CCTACGGGGGGCTGCAGTGGGGAATATTGGGCAATGGGGGAAACCCTGACCCAGCAACGCCGCGTGAAGGAAGAAGGCCTTCGGGTTGTAAACTTCTTTTACCAGGGACGAAGAACGTGACGGTACCTGGAGAAAAAGCCACGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGTGGCAAGCGTTGTCCGGATTTACTGGGTGTAAAGGGCGTGTAGGCGGAGCTGCAAGTCAGATGTGAAATCC 
    ##                                                                                                                                                                                                                                                          7 
    ## CCTACGGGCGGCTGCAGTGGGGAATATTGCACAATGGGGGAAACCCTGATGCAGCAACGCCGCGTGAGTGAAGAAGTATTTCGGTATGTAAAGCTCTATCAGCAGGGAAGAAGATGACGGTACCTGATTAAGAAGCCCCGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGGGGCAAGCGTTATCCGGATTTACTGGGTGTAAAGGGAGCGCAGACGGCAGTGCAAGTCTGGAGTGAAAGCCCGG 
    ##                                                                                                                                                                                                                                                          8 
    ## GACTACTGGGGTATCTAAGCCTGTTTGCTCCCCACGCTTTCGTGCCTCAGTGTCAGTTACAGTCCAGCAACTCGCCTTCGCCACCGGTGTTCTTCCTAATATCTACGCATTCCACCGCTACACTAGGAATTCCAGTTGCCCCTCCTGCACTCAAGTCCAACAGTTTTAGTAGTAGTGCCGGAGTTGAGCCCCGGAGTTACGCTACTAACTTGCTGAACCACCTACGCACCCTTTACGCCCAGTCATTCCG 
    ##                                                                                                                                                                                                                                                         39 
    ## CCTACGGGTGGCTGCAGTGAGGAATATTGGTCAATGGCCGGAAGGCTGAACCAGCCAAGTCGCGTGAGGGAATAAGGCCCTACGGGTCGTAAACCTCTTTTGTCAGGGAGCAAGGCCGCCCACGTGTGGACGGAAGGAGAGTACCTGGAGAAAAAGCATCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGTAGGCGGCCTGCC 
    ##                                                                                                                                                                                                                                                         30

``` r
# Inspect the mergers data.frame for the first sample
head(mergers[[1]])
```

    ##                                                                                                                                                                                                                                                                                                                                                                                                                                                  sequence
    ## 1                   TGAGGAATATTGGTCAATGGGCGCGAGCCTGAACCAGCCAAGTCGCGTGAGGGAAGACGGTCCTATGGATTGTAAACCTCTTTTGTCGGGGAGCAAAAGACGCCACGCGTGGCGTTCCGAGAGTACCCGAAGAAAAAGCATCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGCAGGCGGAAGATCAAGTCAGCGGTAAAATTGAGAGGCTCAACCTCTTCGAGCCGTTGAAACTGGTTTTCTTGAGTGAGCGAGAAGTATGCGGAATGCGTGGTGTAGCGGTGAAATGCATAGATATCACGCAGAACTCCGATTGCGAAGGCAGCATACCGGCGCTCAACTGACGCTCATGCACGAAAGTGTGGGTATCGAACA
    ## 2  CCTACGGGGGGCTGCAGTGAGGAATATTGGTCAATGGGCGCGAGCCTGAACCAGCCAAGTCGCGTGAGGGAAGACGGTCCTATGGATTGTAAACCTCTTTTGTCGGGGAGCAAAAGACGCCACGCGTGGCGTTCCGAGAGTACCCGAAGAAAAAGCATCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGCAGGCGGAAGATCAAGTCAGCGGTAAAATTGAGAGGCTCAACCTCTTCGAGCCGTTGAAACTGGTTTTCTTGAGTGAGCGAGAAGTATGCGGAATGCGTGGTGTAGCGGTGAAATGCATAGATATCACGCAGAACTCCGATTGCGAAGGCAGCATACCGGCGCTCAACTGACGCTCATGCACGAAAGTGTGGGTATCGAACA
    ## 3                  TGAGGAATATTGGTCAATGGACGAGAGTCTGAACCAGCCAAGTAGCGTGAAGGATGACTGCCCTATGGGTTGTAAACTTCTTTTATATGGGAATAAAATGTTCCACGTGTGGGATTTTGTATGTACCATATGAATAAGGATCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGATCCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGAGCGTAGGTGGATTGTTAAGTCAGTTGTGAAAGTTTGCGGCTCAACCGTAAAATTGCAGTTGAAACTGGCAGTCTTGAGTACAGTAGAGGTGGGCGGAATTCGTGGTGTAGCGGTGAAATGCTTAGATATCACGAAGAACTCCGATTGCGAAGGCAGCTCACTGGACTGCAACTGACACTGAGGCTCGAAAGTGTGGGTATCAAACA
    ## 4                  TGAGGAATATTGGTCAATGGGCGCGAGCCTGAACCAGCCAAGTAGCGTGCAGGACGACGGCCCTATGGGTTGTAAACTGCTTTTATACGGGGATAAAGTATGCCACGTGTGGTTTATTGCAGGTACCGTATGAATAAGGACCGGCTAATTCCGTGCCAGCAGCCGCGGTAATACGGAAGGTCCGGGCGTTATCCGGATTTATTGGGTTTAAAGGGAGCGTAGGCCGGGTTTTAAGCGTGCCGTGAAATGTCGGGGCTCAACCTTGACACTGCGGCGCGAACTGGAGTCCTTGAGTGCGCGGAACGTATGCGGAATTCGTGGTGTAGCGGTGAAATGCTTAGATATCACGAAGAACCCCGATTGCGAAGGCAGCATACGGCAGCGCTACTGACGCTGAAGCTCGAAGGCGCGGGTATCGAACA
    ## 5                                     TCGGGAATATTGCGCAATGGAGGAAACTCTGACGCAGTGACGCCGCGTATAGGAAGAAGGTTTTCGGATTGTAAACTATTGTCGTTAGGGAAGATAAAAGACAGTACCTAAGGAGGAAGCTCCGGCTAACTATGTGCCAGCAGCCGCGGTAATACATAGGGAGCAAGCGTTATCCGGAATTATTGGGTGTAAAGGGTGCGTAGACGGAGGAACAAGTTAGTTGTGAAATCCCTCGGCTTAACTGAGGAACTGCAACTAAAACTATTCCCCTTGAGTGTCGGAGAGGAAAGTGGAATTCCTAGTGTAGCGGTGAAATGCGTAGATATTAGGAGGAACACCAGTGGCGAAGGCGACTTTCTGGACGATAACTGACGTTGAGGCACGAAAGTGTGGGGAGCAAACA
    ## 6 CCTACGGGGGGCTGCAGTGAGGAATATTGGTCAATGGACGAGAGTCTGAACCAGCCAAGTAGCGTGAAGGATGACTGCCCTATGGGTTGTAAACTTCTTTTATATGGGAATAAAATGTTCCACGTGTGGGATTTTGTATGTACCATATGAATAAGGATCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGATCCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGAGCGTAGGTGGATTGTTAAGTCAGTTGTGAAAGTTTGCGGCTCAACCGTAAAATTGCAGTTGAAACTGGCAGTCTTGAGTACAGTAGAGGTGGGCGGAATTCGTGGTGTAGCGGTGAAATGCTTAGATATCACGAAGAACTCCGATTGCGAAGGCAGCTCACTGGACTGCAACTGACACTGAGGCTCGAAAGTGTGGGTATCAAACA
    ##   abundance forward reverse nmatch nmismatch nindel prefer accept
    ## 1      1604       1       1     49         0      0      2   TRUE
    ## 2      1148       2       1     32         0      0      2   TRUE
    ## 3       761       3       2     48         0      0      2   TRUE
    ## 4       669       4       3     48         0      0      2   TRUE
    ## 5       628       5       4     67         0      0      2   TRUE
    ## 6       560       7       2     31         0      0      2   TRUE

``` r
# Remove intermediate files to save memory
rm(derepF); rm(derepR)
```

------------------------------------------------------------------------

## Step 6: Sequence Table Construction

``` r
seqtab <- makeSequenceTable(mergers)

# Number of sequence variants
dim(seqtab)
```

    ## [1]   2 839

``` r
# Size distribution of sequence variants
table(nchar(getSequences(seqtab)))
```

    ## 
    ## 401 402 403 404 405 406 408 419 420 421 422 423 424 425 426 427 437 438 439 440 
    ##   2  26  12   6  14   1   1  27  19  72  77   2   1   3   2   3   7 112 105  33 
    ## 441 442 443 444 446 458 
    ## 103   9  14   6 112  70

``` r
# Save sequence table
saveRDS(seqtab, "seqtab.rds")
```

------------------------------------------------------------------------

## Step 7: Chimera Removal

Search for and remove chimeric sequences from the sequence table.

``` r
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
```

    ## Identified 504 bimeras out of 839 input sequences.

``` r
# Check the impact of chimera removal
dim(seqtab.nochim)
```

    ## [1]   2 335

``` r
sum(seqtab.nochim)/sum(seqtab)
```

    ## [1] 0.6426458

``` r
# Size distribution after chimera removal
table(nchar(getSequences(seqtab.nochim)))
```

    ## 
    ## 401 402 403 404 405 406 408 419 420 421 422 423 426 427 437 438 439 440 441 442 
    ##   2  22  12   6  14   1   1  19  12  23  38   1   2   3   1  15  29  26  33   9 
    ## 443 444 446 458 
    ##   8   6  13  39

------------------------------------------------------------------------

## Step 8: Sequence Filtering

Filter sequence variants by length, retaining only those of biologically
relevant sizes. Here 400 to 431 bp sequences are being kept, but you
will want to change based on the size distribution you get after chimera
removal. Save the ASV table

``` r
# Filter sequence variants between 400 and 431 bp
seqtab2.nochim <- seqtab.nochim[, nchar(colnames(seqtab.nochim)) %in% seq(400, 431)]
dim(seqtab2.nochim)
```

    ## [1]   2 156

``` r
table(nchar(getSequences(seqtab2.nochim)))
```

    ## 
    ## 401 402 403 404 405 406 408 419 420 421 422 423 426 427 
    ##   2  22  12   6  14   1   1  19  12  23  38   1   2   3

``` r
# Save the sequence variant table
seqnum <- paste0("ASV", seq(ncol(seqtab2.nochim))) 
uniqueSeqs <- as.list(colnames(seqtab2.nochim))   
seqtab2.nochim.transposed <- t(seqtab2.nochim)  
rownames(seqtab2.nochim.transposed) <- as.character(seqnum) 
write.csv(seqtab2.nochim.transposed, file="ASV_abundance.csv") 
```

------------------------------------------------------------------------

## Step 9: Write ASV Sequences to FASTA

Export ASV sequences to a FASTA file for downstream analysis.

``` r
write.fasta(uniqueSeqs, seqnum, "Asvs.fasta")
```

------------------------------------------------------------------------

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
```

    ##                                                                                                                                                                                                                                                                                                                                                                                                                                        Kingdom   
    ## TGAGGAATATTGGTCAATGGGCGCGAGCCTGAACCAGCCAAGTCGCGTGAGGGAAGACGGTCCTATGGATTGTAAACCTCTTTTGTCGGGGAGCAAAAGACGCCACGCGTGGCGTTCCGAGAGTACCCGAAGAAAAAGCATCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGCAGGCGGAAGATCAAGTCAGCGGTAAAATTGAGAGGCTCAACCTCTTCGAGCCGTTGAAACTGGTTTTCTTGAGTGAGCGAGAAGTATGCGGAATGCGTGGTGTAGCGGTGAAATGCATAGATATCACGCAGAACTCCGATTGCGAAGGCAGCATACCGGCGCTCAACTGACGCTCATGCACGAAAGTGTGGGTATCGAACA  "Bacteria"
    ## TGAGGAATATTGGTCAATGGACGAGAGTCTGAACCAGCCAAGTAGCGTGAAGGATGACTGCCCTATGGGTTGTAAACTTCTTTTATATGGGAATAAAATGTTCCACGTGTGGGATTTTGTATGTACCATATGAATAAGGATCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGATCCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGAGCGTAGGTGGATTGTTAAGTCAGTTGTGAAAGTTTGCGGCTCAACCGTAAAATTGCAGTTGAAACTGGCAGTCTTGAGTACAGTAGAGGTGGGCGGAATTCGTGGTGTAGCGGTGAAATGCTTAGATATCACGAAGAACTCCGATTGCGAAGGCAGCTCACTGGACTGCAACTGACACTGAGGCTCGAAAGTGTGGGTATCAAACA "Bacteria"
    ## TGAGGAATATTGGTCAATGGGCGCGAGCCTGAACCAGCCAAGTAGCGTGCAGGACGACGGCCCTATGGGTTGTAAACTGCTTTTATACGGGGATAAAGTATGCCACGTGTGGTTTATTGCAGGTACCGTATGAATAAGGACCGGCTAATTCCGTGCCAGCAGCCGCGGTAATACGGAAGGTCCGGGCGTTATCCGGATTTATTGGGTTTAAAGGGAGCGTAGGCCGGGTTTTAAGCGTGCCGTGAAATGTCGGGGCTCAACCTTGACACTGCGGCGCGAACTGGAGTCCTTGAGTGCGCGGAACGTATGCGGAATTCGTGGTGTAGCGGTGAAATGCTTAGATATCACGAAGAACCCCGATTGCGAAGGCAGCATACGGCAGCGCTACTGACGCTGAAGCTCGAAGGCGCGGGTATCGAACA "Bacteria"
    ## TCGAGAATCATTCACAATGGGGGAAACCCTGATGGTGCGACGCCGCGTGGGGGAATGAAGGTCTTCGGATTGTAAACCCCTGTCATGTGGGAGCAAATTAAAAAGATAGTACCACAAGAGGAAGAGACGGCTAACTCTGTGCCAGCAGCCGCGGTAATACAGAGGTCTCAAGCGTTGTTCGGAATCACTGGGCGTAAAGCGTGCGTAGGCTGTTTCGTAAGTCGTGTGTGAAAGGCGCGGGCTCAACCCGCGGACGGCACATGATACTGCGAGACTAGAGTAATGGAGGGGGAACCGGAATTCTCGGTGTAGCAGTGAAATGCGTAGATATCGAGAGGAACACTCGTGGCGAAGGCGGGTTCCTGGACATTAACTGACGCTGAGGCACGAAGGCCAGGGGAGCGAAAG               "Bacteria"
    ## TCGGGAATATTGCGCAATGGAGGAAACTCTGACGCAGTGACGCCGCGTATAGGAAGAAGGTTTTCGGATTGTAAACTATTGTCGTTAGGGAAGATAAAAGACAGTACCTAAGGAGGAAGCTCCGGCTAACTATGTGCCAGCAGCCGCGGTAATACATAGGGAGCAAGCGTTATCCGGAATTATTGGGTGTAAAGGGTGCGTAGACGGAGGAACAAGTTAGTTGTGAAATCCCTCGGCTTAACTGAGGAACTGCAACTAAAACTATTCCCCTTGAGTGTCGGAGAGGAAAGTGGAATTCCTAGTGTAGCGGTGAAATGCGTAGATATTAGGAGGAACACCAGTGGCGAAGGCGACTTTCTGGACGATAACTGACGTTGAGGCACGAAAGTGTGGGGAGCAAACA                    "Bacteria"
    ## TGAGGAATATTGGTCAATGGCCGGAAGGCTGAACCAGCCAAGTCGCGTGAGGGAATAAGGCCCTACGGGTCGTAAACCTCTTTTGCCAGGGAGCAATGCCGTCCACGTGTGGACGGAAGGAGAGTACCTGGAGAAAAAGCATCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGTAGGCGGCCTGCCAAGTCAGCGGTAAAATTGCGGGGCTCAACCCCGTACAGCCGTTGAAACTGCCGGGCTCGAGTGGGCGAGAAGTATGCGGAATGCGTGGTGTAGCGGTGAAATGCATAGATATCACGCAGAACCCCGATTGCGAAGGCAGCATACCGGCGCCCGACTGACGCTGAGGCACGAAAGTGCGGGGATCAAACA "Bacteria"
    ##                                                                                                                                                                                                                                                                                                                                                                                                                                        Phylum             
    ## TGAGGAATATTGGTCAATGGGCGCGAGCCTGAACCAGCCAAGTCGCGTGAGGGAAGACGGTCCTATGGATTGTAAACCTCTTTTGTCGGGGAGCAAAAGACGCCACGCGTGGCGTTCCGAGAGTACCCGAAGAAAAAGCATCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGCAGGCGGAAGATCAAGTCAGCGGTAAAATTGAGAGGCTCAACCTCTTCGAGCCGTTGAAACTGGTTTTCTTGAGTGAGCGAGAAGTATGCGGAATGCGTGGTGTAGCGGTGAAATGCATAGATATCACGCAGAACTCCGATTGCGAAGGCAGCATACCGGCGCTCAACTGACGCTCATGCACGAAAGTGTGGGTATCGAACA  "Bacteroidota"     
    ## TGAGGAATATTGGTCAATGGACGAGAGTCTGAACCAGCCAAGTAGCGTGAAGGATGACTGCCCTATGGGTTGTAAACTTCTTTTATATGGGAATAAAATGTTCCACGTGTGGGATTTTGTATGTACCATATGAATAAGGATCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGATCCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGAGCGTAGGTGGATTGTTAAGTCAGTTGTGAAAGTTTGCGGCTCAACCGTAAAATTGCAGTTGAAACTGGCAGTCTTGAGTACAGTAGAGGTGGGCGGAATTCGTGGTGTAGCGGTGAAATGCTTAGATATCACGAAGAACTCCGATTGCGAAGGCAGCTCACTGGACTGCAACTGACACTGAGGCTCGAAAGTGTGGGTATCAAACA "Bacteroidota"     
    ## TGAGGAATATTGGTCAATGGGCGCGAGCCTGAACCAGCCAAGTAGCGTGCAGGACGACGGCCCTATGGGTTGTAAACTGCTTTTATACGGGGATAAAGTATGCCACGTGTGGTTTATTGCAGGTACCGTATGAATAAGGACCGGCTAATTCCGTGCCAGCAGCCGCGGTAATACGGAAGGTCCGGGCGTTATCCGGATTTATTGGGTTTAAAGGGAGCGTAGGCCGGGTTTTAAGCGTGCCGTGAAATGTCGGGGCTCAACCTTGACACTGCGGCGCGAACTGGAGTCCTTGAGTGCGCGGAACGTATGCGGAATTCGTGGTGTAGCGGTGAAATGCTTAGATATCACGAAGAACCCCGATTGCGAAGGCAGCATACGGCAGCGCTACTGACGCTGAAGCTCGAAGGCGCGGGTATCGAACA "Bacteroidota"     
    ## TCGAGAATCATTCACAATGGGGGAAACCCTGATGGTGCGACGCCGCGTGGGGGAATGAAGGTCTTCGGATTGTAAACCCCTGTCATGTGGGAGCAAATTAAAAAGATAGTACCACAAGAGGAAGAGACGGCTAACTCTGTGCCAGCAGCCGCGGTAATACAGAGGTCTCAAGCGTTGTTCGGAATCACTGGGCGTAAAGCGTGCGTAGGCTGTTTCGTAAGTCGTGTGTGAAAGGCGCGGGCTCAACCCGCGGACGGCACATGATACTGCGAGACTAGAGTAATGGAGGGGGAACCGGAATTCTCGGTGTAGCAGTGAAATGCGTAGATATCGAGAGGAACACTCGTGGCGAAGGCGGGTTCCTGGACATTAACTGACGCTGAGGCACGAAGGCCAGGGGAGCGAAAG               "Verrucomicrobiota"
    ## TCGGGAATATTGCGCAATGGAGGAAACTCTGACGCAGTGACGCCGCGTATAGGAAGAAGGTTTTCGGATTGTAAACTATTGTCGTTAGGGAAGATAAAAGACAGTACCTAAGGAGGAAGCTCCGGCTAACTATGTGCCAGCAGCCGCGGTAATACATAGGGAGCAAGCGTTATCCGGAATTATTGGGTGTAAAGGGTGCGTAGACGGAGGAACAAGTTAGTTGTGAAATCCCTCGGCTTAACTGAGGAACTGCAACTAAAACTATTCCCCTTGAGTGTCGGAGAGGAAAGTGGAATTCCTAGTGTAGCGGTGAAATGCGTAGATATTAGGAGGAACACCAGTGGCGAAGGCGACTTTCTGGACGATAACTGACGTTGAGGCACGAAAGTGTGGGGAGCAAACA                    "Firmicutes"       
    ## TGAGGAATATTGGTCAATGGCCGGAAGGCTGAACCAGCCAAGTCGCGTGAGGGAATAAGGCCCTACGGGTCGTAAACCTCTTTTGCCAGGGAGCAATGCCGTCCACGTGTGGACGGAAGGAGAGTACCTGGAGAAAAAGCATCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGTAGGCGGCCTGCCAAGTCAGCGGTAAAATTGCGGGGCTCAACCCCGTACAGCCGTTGAAACTGCCGGGCTCGAGTGGGCGAGAAGTATGCGGAATGCGTGGTGTAGCGGTGAAATGCATAGATATCACGCAGAACCCCGATTGCGAAGGCAGCATACCGGCGCCCGACTGACGCTGAGGCACGAAAGTGCGGGGATCAAACA "Bacteroidota"     
    ##                                                                                                                                                                                                                                                                                                                                                                                                                                        Class             
    ## TGAGGAATATTGGTCAATGGGCGCGAGCCTGAACCAGCCAAGTCGCGTGAGGGAAGACGGTCCTATGGATTGTAAACCTCTTTTGTCGGGGAGCAAAAGACGCCACGCGTGGCGTTCCGAGAGTACCCGAAGAAAAAGCATCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGCAGGCGGAAGATCAAGTCAGCGGTAAAATTGAGAGGCTCAACCTCTTCGAGCCGTTGAAACTGGTTTTCTTGAGTGAGCGAGAAGTATGCGGAATGCGTGGTGTAGCGGTGAAATGCATAGATATCACGCAGAACTCCGATTGCGAAGGCAGCATACCGGCGCTCAACTGACGCTCATGCACGAAAGTGTGGGTATCGAACA  "Bacteroidia"     
    ## TGAGGAATATTGGTCAATGGACGAGAGTCTGAACCAGCCAAGTAGCGTGAAGGATGACTGCCCTATGGGTTGTAAACTTCTTTTATATGGGAATAAAATGTTCCACGTGTGGGATTTTGTATGTACCATATGAATAAGGATCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGATCCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGAGCGTAGGTGGATTGTTAAGTCAGTTGTGAAAGTTTGCGGCTCAACCGTAAAATTGCAGTTGAAACTGGCAGTCTTGAGTACAGTAGAGGTGGGCGGAATTCGTGGTGTAGCGGTGAAATGCTTAGATATCACGAAGAACTCCGATTGCGAAGGCAGCTCACTGGACTGCAACTGACACTGAGGCTCGAAAGTGTGGGTATCAAACA "Bacteroidia"     
    ## TGAGGAATATTGGTCAATGGGCGCGAGCCTGAACCAGCCAAGTAGCGTGCAGGACGACGGCCCTATGGGTTGTAAACTGCTTTTATACGGGGATAAAGTATGCCACGTGTGGTTTATTGCAGGTACCGTATGAATAAGGACCGGCTAATTCCGTGCCAGCAGCCGCGGTAATACGGAAGGTCCGGGCGTTATCCGGATTTATTGGGTTTAAAGGGAGCGTAGGCCGGGTTTTAAGCGTGCCGTGAAATGTCGGGGCTCAACCTTGACACTGCGGCGCGAACTGGAGTCCTTGAGTGCGCGGAACGTATGCGGAATTCGTGGTGTAGCGGTGAAATGCTTAGATATCACGAAGAACCCCGATTGCGAAGGCAGCATACGGCAGCGCTACTGACGCTGAAGCTCGAAGGCGCGGGTATCGAACA "Bacteroidia"     
    ## TCGAGAATCATTCACAATGGGGGAAACCCTGATGGTGCGACGCCGCGTGGGGGAATGAAGGTCTTCGGATTGTAAACCCCTGTCATGTGGGAGCAAATTAAAAAGATAGTACCACAAGAGGAAGAGACGGCTAACTCTGTGCCAGCAGCCGCGGTAATACAGAGGTCTCAAGCGTTGTTCGGAATCACTGGGCGTAAAGCGTGCGTAGGCTGTTTCGTAAGTCGTGTGTGAAAGGCGCGGGCTCAACCCGCGGACGGCACATGATACTGCGAGACTAGAGTAATGGAGGGGGAACCGGAATTCTCGGTGTAGCAGTGAAATGCGTAGATATCGAGAGGAACACTCGTGGCGAAGGCGGGTTCCTGGACATTAACTGACGCTGAGGCACGAAGGCCAGGGGAGCGAAAG               "Verrucomicrobiae"
    ## TCGGGAATATTGCGCAATGGAGGAAACTCTGACGCAGTGACGCCGCGTATAGGAAGAAGGTTTTCGGATTGTAAACTATTGTCGTTAGGGAAGATAAAAGACAGTACCTAAGGAGGAAGCTCCGGCTAACTATGTGCCAGCAGCCGCGGTAATACATAGGGAGCAAGCGTTATCCGGAATTATTGGGTGTAAAGGGTGCGTAGACGGAGGAACAAGTTAGTTGTGAAATCCCTCGGCTTAACTGAGGAACTGCAACTAAAACTATTCCCCTTGAGTGTCGGAGAGGAAAGTGGAATTCCTAGTGTAGCGGTGAAATGCGTAGATATTAGGAGGAACACCAGTGGCGAAGGCGACTTTCTGGACGATAACTGACGTTGAGGCACGAAAGTGTGGGGAGCAAACA                    "Clostridia"      
    ## TGAGGAATATTGGTCAATGGCCGGAAGGCTGAACCAGCCAAGTCGCGTGAGGGAATAAGGCCCTACGGGTCGTAAACCTCTTTTGCCAGGGAGCAATGCCGTCCACGTGTGGACGGAAGGAGAGTACCTGGAGAAAAAGCATCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGTAGGCGGCCTGCCAAGTCAGCGGTAAAATTGCGGGGCTCAACCCCGTACAGCCGTTGAAACTGCCGGGCTCGAGTGGGCGAGAAGTATGCGGAATGCGTGGTGTAGCGGTGAAATGCATAGATATCACGCAGAACCCCGATTGCGAAGGCAGCATACCGGCGCCCGACTGACGCTGAGGCACGAAAGTGCGGGGATCAAACA "Bacteroidia"     
    ##                                                                                                                                                                                                                                                                                                                                                                                                                                        Order               
    ## TGAGGAATATTGGTCAATGGGCGCGAGCCTGAACCAGCCAAGTCGCGTGAGGGAAGACGGTCCTATGGATTGTAAACCTCTTTTGTCGGGGAGCAAAAGACGCCACGCGTGGCGTTCCGAGAGTACCCGAAGAAAAAGCATCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGCAGGCGGAAGATCAAGTCAGCGGTAAAATTGAGAGGCTCAACCTCTTCGAGCCGTTGAAACTGGTTTTCTTGAGTGAGCGAGAAGTATGCGGAATGCGTGGTGTAGCGGTGAAATGCATAGATATCACGCAGAACTCCGATTGCGAAGGCAGCATACCGGCGCTCAACTGACGCTCATGCACGAAAGTGTGGGTATCGAACA  "Bacteroidales"     
    ## TGAGGAATATTGGTCAATGGACGAGAGTCTGAACCAGCCAAGTAGCGTGAAGGATGACTGCCCTATGGGTTGTAAACTTCTTTTATATGGGAATAAAATGTTCCACGTGTGGGATTTTGTATGTACCATATGAATAAGGATCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGATCCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGAGCGTAGGTGGATTGTTAAGTCAGTTGTGAAAGTTTGCGGCTCAACCGTAAAATTGCAGTTGAAACTGGCAGTCTTGAGTACAGTAGAGGTGGGCGGAATTCGTGGTGTAGCGGTGAAATGCTTAGATATCACGAAGAACTCCGATTGCGAAGGCAGCTCACTGGACTGCAACTGACACTGAGGCTCGAAAGTGTGGGTATCAAACA "Bacteroidales"     
    ## TGAGGAATATTGGTCAATGGGCGCGAGCCTGAACCAGCCAAGTAGCGTGCAGGACGACGGCCCTATGGGTTGTAAACTGCTTTTATACGGGGATAAAGTATGCCACGTGTGGTTTATTGCAGGTACCGTATGAATAAGGACCGGCTAATTCCGTGCCAGCAGCCGCGGTAATACGGAAGGTCCGGGCGTTATCCGGATTTATTGGGTTTAAAGGGAGCGTAGGCCGGGTTTTAAGCGTGCCGTGAAATGTCGGGGCTCAACCTTGACACTGCGGCGCGAACTGGAGTCCTTGAGTGCGCGGAACGTATGCGGAATTCGTGGTGTAGCGGTGAAATGCTTAGATATCACGAAGAACCCCGATTGCGAAGGCAGCATACGGCAGCGCTACTGACGCTGAAGCTCGAAGGCGCGGGTATCGAACA "Bacteroidales"     
    ## TCGAGAATCATTCACAATGGGGGAAACCCTGATGGTGCGACGCCGCGTGGGGGAATGAAGGTCTTCGGATTGTAAACCCCTGTCATGTGGGAGCAAATTAAAAAGATAGTACCACAAGAGGAAGAGACGGCTAACTCTGTGCCAGCAGCCGCGGTAATACAGAGGTCTCAAGCGTTGTTCGGAATCACTGGGCGTAAAGCGTGCGTAGGCTGTTTCGTAAGTCGTGTGTGAAAGGCGCGGGCTCAACCCGCGGACGGCACATGATACTGCGAGACTAGAGTAATGGAGGGGGAACCGGAATTCTCGGTGTAGCAGTGAAATGCGTAGATATCGAGAGGAACACTCGTGGCGAAGGCGGGTTCCTGGACATTAACTGACGCTGAGGCACGAAGGCCAGGGGAGCGAAAG               "Verrucomicrobiales"
    ## TCGGGAATATTGCGCAATGGAGGAAACTCTGACGCAGTGACGCCGCGTATAGGAAGAAGGTTTTCGGATTGTAAACTATTGTCGTTAGGGAAGATAAAAGACAGTACCTAAGGAGGAAGCTCCGGCTAACTATGTGCCAGCAGCCGCGGTAATACATAGGGAGCAAGCGTTATCCGGAATTATTGGGTGTAAAGGGTGCGTAGACGGAGGAACAAGTTAGTTGTGAAATCCCTCGGCTTAACTGAGGAACTGCAACTAAAACTATTCCCCTTGAGTGTCGGAGAGGAAAGTGGAATTCCTAGTGTAGCGGTGAAATGCGTAGATATTAGGAGGAACACCAGTGGCGAAGGCGACTTTCTGGACGATAACTGACGTTGAGGCACGAAAGTGTGGGGAGCAAACA                    "Clostridia UCG-014"
    ## TGAGGAATATTGGTCAATGGCCGGAAGGCTGAACCAGCCAAGTCGCGTGAGGGAATAAGGCCCTACGGGTCGTAAACCTCTTTTGCCAGGGAGCAATGCCGTCCACGTGTGGACGGAAGGAGAGTACCTGGAGAAAAAGCATCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGTAGGCGGCCTGCCAAGTCAGCGGTAAAATTGCGGGGCTCAACCCCGTACAGCCGTTGAAACTGCCGGGCTCGAGTGGGCGAGAAGTATGCGGAATGCGTGGTGTAGCGGTGAAATGCATAGATATCACGCAGAACCCCGATTGCGAAGGCAGCATACCGGCGCCCGACTGACGCTGAGGCACGAAAGTGCGGGGATCAAACA "Bacteroidales"     
    ##                                                                                                                                                                                                                                                                                                                                                                                                                                        Family           
    ## TGAGGAATATTGGTCAATGGGCGCGAGCCTGAACCAGCCAAGTCGCGTGAGGGAAGACGGTCCTATGGATTGTAAACCTCTTTTGTCGGGGAGCAAAAGACGCCACGCGTGGCGTTCCGAGAGTACCCGAAGAAAAAGCATCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGCAGGCGGAAGATCAAGTCAGCGGTAAAATTGAGAGGCTCAACCTCTTCGAGCCGTTGAAACTGGTTTTCTTGAGTGAGCGAGAAGTATGCGGAATGCGTGGTGTAGCGGTGAAATGCATAGATATCACGCAGAACTCCGATTGCGAAGGCAGCATACCGGCGCTCAACTGACGCTCATGCACGAAAGTGTGGGTATCGAACA  "Muribaculaceae" 
    ## TGAGGAATATTGGTCAATGGACGAGAGTCTGAACCAGCCAAGTAGCGTGAAGGATGACTGCCCTATGGGTTGTAAACTTCTTTTATATGGGAATAAAATGTTCCACGTGTGGGATTTTGTATGTACCATATGAATAAGGATCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGATCCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGAGCGTAGGTGGATTGTTAAGTCAGTTGTGAAAGTTTGCGGCTCAACCGTAAAATTGCAGTTGAAACTGGCAGTCTTGAGTACAGTAGAGGTGGGCGGAATTCGTGGTGTAGCGGTGAAATGCTTAGATATCACGAAGAACTCCGATTGCGAAGGCAGCTCACTGGACTGCAACTGACACTGAGGCTCGAAAGTGTGGGTATCAAACA "Bacteroidaceae" 
    ## TGAGGAATATTGGTCAATGGGCGCGAGCCTGAACCAGCCAAGTAGCGTGCAGGACGACGGCCCTATGGGTTGTAAACTGCTTTTATACGGGGATAAAGTATGCCACGTGTGGTTTATTGCAGGTACCGTATGAATAAGGACCGGCTAATTCCGTGCCAGCAGCCGCGGTAATACGGAAGGTCCGGGCGTTATCCGGATTTATTGGGTTTAAAGGGAGCGTAGGCCGGGTTTTAAGCGTGCCGTGAAATGTCGGGGCTCAACCTTGACACTGCGGCGCGAACTGGAGTCCTTGAGTGCGCGGAACGTATGCGGAATTCGTGGTGTAGCGGTGAAATGCTTAGATATCACGAAGAACCCCGATTGCGAAGGCAGCATACGGCAGCGCTACTGACGCTGAAGCTCGAAGGCGCGGGTATCGAACA "Prevotellaceae" 
    ## TCGAGAATCATTCACAATGGGGGAAACCCTGATGGTGCGACGCCGCGTGGGGGAATGAAGGTCTTCGGATTGTAAACCCCTGTCATGTGGGAGCAAATTAAAAAGATAGTACCACAAGAGGAAGAGACGGCTAACTCTGTGCCAGCAGCCGCGGTAATACAGAGGTCTCAAGCGTTGTTCGGAATCACTGGGCGTAAAGCGTGCGTAGGCTGTTTCGTAAGTCGTGTGTGAAAGGCGCGGGCTCAACCCGCGGACGGCACATGATACTGCGAGACTAGAGTAATGGAGGGGGAACCGGAATTCTCGGTGTAGCAGTGAAATGCGTAGATATCGAGAGGAACACTCGTGGCGAAGGCGGGTTCCTGGACATTAACTGACGCTGAGGCACGAAGGCCAGGGGAGCGAAAG               "Akkermansiaceae"
    ## TCGGGAATATTGCGCAATGGAGGAAACTCTGACGCAGTGACGCCGCGTATAGGAAGAAGGTTTTCGGATTGTAAACTATTGTCGTTAGGGAAGATAAAAGACAGTACCTAAGGAGGAAGCTCCGGCTAACTATGTGCCAGCAGCCGCGGTAATACATAGGGAGCAAGCGTTATCCGGAATTATTGGGTGTAAAGGGTGCGTAGACGGAGGAACAAGTTAGTTGTGAAATCCCTCGGCTTAACTGAGGAACTGCAACTAAAACTATTCCCCTTGAGTGTCGGAGAGGAAAGTGGAATTCCTAGTGTAGCGGTGAAATGCGTAGATATTAGGAGGAACACCAGTGGCGAAGGCGACTTTCTGGACGATAACTGACGTTGAGGCACGAAAGTGTGGGGAGCAAACA                    NA               
    ## TGAGGAATATTGGTCAATGGCCGGAAGGCTGAACCAGCCAAGTCGCGTGAGGGAATAAGGCCCTACGGGTCGTAAACCTCTTTTGCCAGGGAGCAATGCCGTCCACGTGTGGACGGAAGGAGAGTACCTGGAGAAAAAGCATCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGTAGGCGGCCTGCCAAGTCAGCGGTAAAATTGCGGGGCTCAACCCCGTACAGCCGTTGAAACTGCCGGGCTCGAGTGGGCGAGAAGTATGCGGAATGCGTGGTGTAGCGGTGAAATGCATAGATATCACGCAGAACCCCGATTGCGAAGGCAGCATACCGGCGCCCGACTGACGCTGAGGCACGAAAGTGCGGGGATCAAACA "Muribaculaceae" 
    ##                                                                                                                                                                                                                                                                                                                                                                                                                                        Genus                   
    ## TGAGGAATATTGGTCAATGGGCGCGAGCCTGAACCAGCCAAGTCGCGTGAGGGAAGACGGTCCTATGGATTGTAAACCTCTTTTGTCGGGGAGCAAAAGACGCCACGCGTGGCGTTCCGAGAGTACCCGAAGAAAAAGCATCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGCAGGCGGAAGATCAAGTCAGCGGTAAAATTGAGAGGCTCAACCTCTTCGAGCCGTTGAAACTGGTTTTCTTGAGTGAGCGAGAAGTATGCGGAATGCGTGGTGTAGCGGTGAAATGCATAGATATCACGCAGAACTCCGATTGCGAAGGCAGCATACCGGCGCTCAACTGACGCTCATGCACGAAAGTGTGGGTATCGAACA  NA                      
    ## TGAGGAATATTGGTCAATGGACGAGAGTCTGAACCAGCCAAGTAGCGTGAAGGATGACTGCCCTATGGGTTGTAAACTTCTTTTATATGGGAATAAAATGTTCCACGTGTGGGATTTTGTATGTACCATATGAATAAGGATCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGATCCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGAGCGTAGGTGGATTGTTAAGTCAGTTGTGAAAGTTTGCGGCTCAACCGTAAAATTGCAGTTGAAACTGGCAGTCTTGAGTACAGTAGAGGTGGGCGGAATTCGTGGTGTAGCGGTGAAATGCTTAGATATCACGAAGAACTCCGATTGCGAAGGCAGCTCACTGGACTGCAACTGACACTGAGGCTCGAAAGTGTGGGTATCAAACA "Bacteroides"           
    ## TGAGGAATATTGGTCAATGGGCGCGAGCCTGAACCAGCCAAGTAGCGTGCAGGACGACGGCCCTATGGGTTGTAAACTGCTTTTATACGGGGATAAAGTATGCCACGTGTGGTTTATTGCAGGTACCGTATGAATAAGGACCGGCTAATTCCGTGCCAGCAGCCGCGGTAATACGGAAGGTCCGGGCGTTATCCGGATTTATTGGGTTTAAAGGGAGCGTAGGCCGGGTTTTAAGCGTGCCGTGAAATGTCGGGGCTCAACCTTGACACTGCGGCGCGAACTGGAGTCCTTGAGTGCGCGGAACGTATGCGGAATTCGTGGTGTAGCGGTGAAATGCTTAGATATCACGAAGAACCCCGATTGCGAAGGCAGCATACGGCAGCGCTACTGACGCTGAAGCTCGAAGGCGCGGGTATCGAACA "Prevotellaceae UCG-001"
    ## TCGAGAATCATTCACAATGGGGGAAACCCTGATGGTGCGACGCCGCGTGGGGGAATGAAGGTCTTCGGATTGTAAACCCCTGTCATGTGGGAGCAAATTAAAAAGATAGTACCACAAGAGGAAGAGACGGCTAACTCTGTGCCAGCAGCCGCGGTAATACAGAGGTCTCAAGCGTTGTTCGGAATCACTGGGCGTAAAGCGTGCGTAGGCTGTTTCGTAAGTCGTGTGTGAAAGGCGCGGGCTCAACCCGCGGACGGCACATGATACTGCGAGACTAGAGTAATGGAGGGGGAACCGGAATTCTCGGTGTAGCAGTGAAATGCGTAGATATCGAGAGGAACACTCGTGGCGAAGGCGGGTTCCTGGACATTAACTGACGCTGAGGCACGAAGGCCAGGGGAGCGAAAG               "Akkermansia"           
    ## TCGGGAATATTGCGCAATGGAGGAAACTCTGACGCAGTGACGCCGCGTATAGGAAGAAGGTTTTCGGATTGTAAACTATTGTCGTTAGGGAAGATAAAAGACAGTACCTAAGGAGGAAGCTCCGGCTAACTATGTGCCAGCAGCCGCGGTAATACATAGGGAGCAAGCGTTATCCGGAATTATTGGGTGTAAAGGGTGCGTAGACGGAGGAACAAGTTAGTTGTGAAATCCCTCGGCTTAACTGAGGAACTGCAACTAAAACTATTCCCCTTGAGTGTCGGAGAGGAAAGTGGAATTCCTAGTGTAGCGGTGAAATGCGTAGATATTAGGAGGAACACCAGTGGCGAAGGCGACTTTCTGGACGATAACTGACGTTGAGGCACGAAAGTGTGGGGAGCAAACA                    NA                      
    ## TGAGGAATATTGGTCAATGGCCGGAAGGCTGAACCAGCCAAGTCGCGTGAGGGAATAAGGCCCTACGGGTCGTAAACCTCTTTTGCCAGGGAGCAATGCCGTCCACGTGTGGACGGAAGGAGAGTACCTGGAGAAAAAGCATCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGTAGGCGGCCTGCCAAGTCAGCGGTAAAATTGCGGGGCTCAACCCCGTACAGCCGTTGAAACTGCCGGGCTCGAGTGGGCGAGAAGTATGCGGAATGCGTGGTGTAGCGGTGAAATGCATAGATATCACGCAGAACCCCGATTGCGAAGGCAGCATACCGGCGCCCGACTGACGCTGAGGCACGAAAGTGCGGGGATCAAACA NA

``` r
# Create a copy of taxonomy and update missing classifications
tax1 <- as.matrix(tax)
rownames(tax1) <- as.character(seqnum)
y <- which(is.na(tax1) == TRUE) 
tax1[y] <- "Unclassified"       
write.csv(tax1, file="ASV_tax.csv") 
```

------------------------------------------------------------------------

## Step 11: Save Workspace

Save the entire R environment to an .Rdata file for reproducibility.

``` r
save.image(file="dada2.Rdata")
```
