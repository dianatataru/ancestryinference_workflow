# ancestryinference_workflow
Breaking up the steps of ancestryinfer to run highthroughput on large natural complex hybrid dataset

## 1. Aligning all samples to three reference genomes

Use ```bwa mem``` to map reads from hybrid to all parental references independently.

```
       bwa mem -M -R $RG1 $genome1 $read1 $read2 > $sam1
       bwa mem -M -R $RG1 $genome2 $read1 $read2 > $sam2
       bwa mem -M -R $RG1 $genome3 $read1 $read2 > $sam3

#filter for mapping quality>29 and proper pairing with samtools

#remove duplicates with picard

```
## 2. Joint filtering of genomes

Identify reads that do not map uniquely to parental genome and exclude them using ```ngsutils```. 

## 3. Read Counting or Variant Calling

Variant calling with bcftools, then reads matching each parental allele at ancestry informative sites are counted from a samtools mpileup file for each hybrid individual. 


## 4. Thinning to one AIM per read across individuals

Counts for each parental allele at ancestry informative sites are subsampled to thin to one ancestry informative site per read if multiple sites occur within one read. This thinning is performed jointly across individuals such that the same site is retained for all individuals in the dataset.

## 5. AncestryHMM

Requires ```armadillo```.
