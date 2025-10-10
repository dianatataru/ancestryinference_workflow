# ancestryinference_workflow
Breaking up the steps of ancestryinfer to run highthroughput on large natural complex hybrid dataset

## 1. Aligning all samples to three reference genomes

Use ```bwa mem``` to map reads from hybrid to all parental references independently.

## 2. Read Counting or Variant Calling

Identify reads that do not map uniquely to parental genome and exclude them using ```ngsutils```. Then reads matching each parental allele at ancestry informative sites are counted from a samtools mpileup file for each hybrid individual. This is what the pipeline does and it works well for low coverage data but my data is high coverage (23x)

## 3. Thinning to one AIM per read across individuals

Counts for each parental allele at ancestry informative sites are subsampled to thin to one ancestry informative site per read if multiple sites occur within one read. This thinning is performed jointly across individuals such that the same site is retained for all individuals in the dataset.

## 4. AncestryHMM

Requires ```armadillo```.
