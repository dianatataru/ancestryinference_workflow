# ancestryinference_workflow
Alterations to the Schumer Lab's threeway ancestryinfer program (https://github.com/Schumerlab/ancestryinfer) to be run with three Mimulus species (M. guttatus, M. lacinatus, M. nasutus). Requires downloading the ancestryinfer program to run.

Written by: Diana Tataru
Created: June 11, 2025

## 0. Create AIMS: Run genotype_filter.sh

This script is adapted for three species from the following manuscript: Fluctuating reproductive isolation and stable ancestry structure in a fine-scaled mosaic of #hybridizing Mimulus monkeyflowers" by Matthew Farnitano, Keith Karoly, and Andrea Sweigart, and github: https://github.com/mfarnitano/CAC_popgen/#blob/main/reference_panels/genotype_filter.sh  

Located in ```/lustre/project/kferris/Diana_Tataru/lac_nas_gut/AIMS```. This script creates a file with ancestry informative sites (AIMs) from 5 allopatric populations each of three species: *M. guttatus, M.laciniatus,* and *M. nasutus*. Input file for it is a joint genotyped vcf ```lacnasgut_jointgeno.vcf.gz```, created using the Ferris Lab GATK variant calling pipeline, using TOLv5 of *M. guttatus* as the reference. These are the populations used:

|    Species    |  Pop |    SRA     | Longitude  | Latitude  |
|    -------    | ---  | ---------  | ---------  | --------  |
|  M. guttatus  |  TOL | Phyt 551   | -120.63315 | 37.969917 |
|  M. guttatus  |  YVO | Unpub.     | -119.74643 | 37.723367 |
|  M. guttatus  | MAR  | SRX030542  | -123.29445 | 43.4786   |
|  M. guttatus  | LMC  | SRX030680  | -123.083917| 38.863983 |
|  M. guttatus  | IM   | SRR398937  | -122.508783|	45.57571 |
|  M. guttatus  | AHQT | SRX142379  | -110.813   | 44.431    |
| M. laciniatus | OPN  | SRR23709136| -119.4852	 | 37.8107   |
| M. laciniatus	| WLF  | Unpub.	    | -119.59385 | 37.841533 |
| M. laciniatus	| TRT  | SRX19570592| -119.70535 |	37.7165  |
| M. laciniatus	| PER  | SRX19570591| -119.3687	 | 37.055767 |
| M. laciniatus	| HUL  | Unpub.	    | -119.150168| 37.2334366|
| M. nasutus	| SF   | SRR29155563| -121.0225	 | 45.264444 |
| M. nasutus	| DPRN | SRR1259273	| -120.344	 | 37.828    |
| M. nasutus	| KOOT | SRR1259272	| -115.983	 | 48.104    |
| M. nasutus	| NHN  | SRX525051	| -124.16	 | 49.273    |
| M. nasutus	| CACN | SRR1259271	| -121.3667	 | 45.71076  |

To visualize distribution of AIMs genomewide, run ```visualizeAIMs.R ```. This also outputs number of AIMS per chromosome:

| Chr01 | Chr02 | Chr03 | Chr04 | Chr05 | Chr06 | Chr07 | Chr08 | Chr09 | Chr10 | Chr11 | Chr12 | Chr13 | Chr14 |
| ----- | ----- | ----- | ----- | ----- | ----- | ----- | ----- | ----- | ----- | ----- | ----- | ----- | ----- |
|  9260 | 15019 | 11090 | 17851 | 13041 | 16073 | 8709  | 20077 | 10296 | 15111 | 11289 | 14642 | 15133 | 25823 | 

Additional files needed to run this:
Located in ```/project/dtataru/ancestryinfer```. This has multiple input files that need to be created and changed.   
  a.  output ```AIMs_panel15_final.AIMs_counts.txt``` & ```AIMs_panel15_final.AIMs.txt``` 
  b.  Reference genomes for each of three species, I'm using the following, listed in ```/project/dtataru/hybrids/ancestryinfer/reference_genomes/```:  
      1.  MguttatusTOL_551_v5.0.fa  
      2.  Mnasutusvar_SF_822_v2.0.fa  
      3.  WLF47.fasta (consensus genome made by me with high coverage unpub. sequencing data  

Sites fixed for each species:
lac(subset): 52,000
nas:51,659
gut:24,590
total:128,249

First time around ended up with 698623 remaining_sites.txt
Ended up wtih 203414. Need to revisit this filitering step because I'm worries that those many sites that weren't fixed were mostly guttatus and then when I subset I ended up with very few laciniatus sites. Let me investigate in the aims file briefly:

Number of guttatus Alt Counts: 730625 (Ref=907553, total=1,638,178)
awk '{ sum += $6 } END { print sum }' AIMs_panel15_final.AIMs.txt

Number of laciniatus Alt Counts: 1053592 (Ref=642072, total=1,695,664)
awk '{ sum += $8 } END { print sum }' AIMs_panel15_final.AIMs.txt

Number of nasutus Alt Counts: 1576401 (Ref- 254823, total= 1,831,224)
awk '{ sum += $11 } END { print sum }' AIMs_panel15_final.AIMs.txt

Okay, so something is off here, still, and I think it has to do with the way sites were subset or some way that ref/alt is being called. I think I want to try to run this with the full set up AIMS that I got, with 1267995 total sites.

```
cd /project/dtataru/lac_nas_gut/AIMS/VCFs

awk 'BEGIN{OFS="\t"} {gsub(/Chr_/,"Chr-",$1); print $1, $2, $3, $4}' \
    panel15.TOL551.SNPs.fspi.rm.AIMs_counts.filtered.txt \
    > /project/dtataru/ancestryinfer/AIMs_panel15_final.AIMs.full.txt

awk 'BEGIN{OFS="\t"} {gsub(/Chr_/,"Chr-",$1); print $1, $2, $5, $6, $7, $8}' \
    panel15.TOL551.SNPs.fspi.rm.AIMs_counts.filtered.txt \
    > /project/dtataru/ancestryinfer/AIMs_panel15_final.AIMs_counts.full.txt
```
	  
## 1. Aligning all samples to three reference genomes

Use ```bwa mem``` to map reads from each hybrid individual to all parental references independently. Input files are fastq files for all samples, four lanes per sample and two reads per lane (Total 1,232 arrays).  In the TMPDIR, this creates a folder for each sample, with all reads aligned for each lane run. Output is three .sam files for each sample, corresponding to each species' reference genome. In the next step, these separate runs will be merged. This script is ```map_array_DT.sh ```.

```
#!/bin/bash
#SBATCH --output=/project/dtataru/ancestryinfer/logs/map_%A_%a.out
#SBATCH --error=/project/dtataru/ancestryinfer/logs/map_%A_%a.err
#SBATCH --time=1-00:00:00
#SBATCH -p single
#SBATCH -N 1            #: Number of Nodes
#SBATCH -n 1            #: Number of Tasks per Node
#SBATCH -A loni_ferrislac
#SBATCH --cpus-per-task=12
#SBATCH --array=1-4   # Job array when n is number of unique samples, % is array throttling (do 32 at a time, 8 samples)

### LOAD MODULES ###
#For this step, bwa and needed
eval "$(conda shell.bash hook)"
conda activate /home/dtataru/.conda/envs/ancestryinfer

echo "Start Job"
echo "SLURM_ARRAY_TASK_ID = ${SLURM_ARRAY_TASK_ID}"

### ASSIGN VARIABLES ARRAY SIZE > 1000 ###
mapfile -t R1_FILE_LIST < <(ls -1 /project/dtataru/hybrids/1_hybrid1data/*R1_001.fastq.gz | sort)
mapfile -t R2_FILE_LIST < <(ls -1 /project/dtataru/hybrids/1_hybrid1data/*R2_001.fastq.gz | sort)

R1_FILE_LIST_SIZE=${#R1_FILE_LIST[@]}
R2_FILE_LIST_SIZE=${#R2_FILE_LIST[@]}

# check for matching read files
for v in R1_FILE_LIST_SIZE R2_FILE_LIST_SIZE; do echo "Using $v=${!v}"; done
if [ $R1_FILE_LIST_SIZE -ne $R2_FILE_LIST_SIZE ]; then
        echo "ERROR: Encountered unmatched read files: R1_FILE_LIST_SIZE=$R1_FILE_LIST_SIZE and R2_FILE_LIST_SIZE=$R2_FILE_LIST._SIZE Exiting."
        exit 1
fi

NUMBER_FILE_PAIRS_ALREADY_DONE=1000
#NUMBER_FILE_PAIRS_ALREADY_DONE=0
NUMBER_FILE_PAIRS_TODO=$(( $R1_FILE_LIST_SIZE-$NUMBER_FILE_PAIRS_ALREADY_DONE ))
MAX_FILE_PAIRS_PER_JOB=$(( ($NUMBER_FILE_PAIRS_TODO+1000-1)/1000 ))
for v in NUMBER_FILE_PAIRS_ALREADY_DONE NUMBER_FILE_PAIRS_TODO MAX_FILE_PAIRS_PER_JOB; do
        echo "Using $v=${!v}"
done
START_FILE_INDEX=$(( $NUMBER_FILE_PAIRS_ALREADY_DONE+1+($SLURM_ARRAY_TASK_ID-1)*$MAX_FILE_PAIRS_PER_JOB))
END_FILE_INDEX=$(( $START_FILE_INDEX+$MAX_FILE_PAIRS_PER_JOB-1 ))
for v in START_FILE_INDEX END_FILE_INDEX; do echo "Using $v=${!v}"; done
for FILE_INDEX in $(seq $START_FILE_INDEX $END_FILE_INDEX); do
        echo "FILE_INDEX=$FILE_INDEX"
        if [ $FILE_INDEX -gt $R1_FILE_LIST_SIZE ]; then
                continue
        fi

R1=${R1_FILE_LIST[$((FILE_INDEX - 1))]}
R2=${R2_FILE_LIST[$((FILE_INDEX - 1))]}

SAMPLE=$(echo $R1 | cut -d "/" -f 6 | cut -d "_" -f 1-3)
HEADER=$(echo $R1 | cut -d "/" -f 6 | cut -d "_" -f 1)

genome1="/project/dtataru/hybrids/ancestryinfer/reference_genomes/MguttatusTOL_551_v5.0.fa"
genome2="/project/dtataru/hybrids/ancestryinfer/reference_genomes/Mnasutusvar_SF_822_v2.0.fa"
genome3="/project/dtataru/hybrids/ancestryinfer/reference_genomes/WLF47.fasta"

echo "R1=$R1"
echo "R2=$R2"
echo "SAMPLE=$SAMPLE"
echo "HEADER=$HEADER"
echo "genome1=$genome1"
echo "genome2=$genome2"
echo "genome3=$genome3"
echo "tag=$tag"

### SET TMPDIR ###
TMPDIR="/work/dtataru/TMPDIR/${HEADER}"
mkdir -p "$TMPDIR"
cd "$TMPDIR" 

### MAPPING ###
echo "Mapping ${SAMPLE} to three parental genomes"

# Read group string
RG="@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tPL:illumina\tLB:hyblib1\tPU:LSIslowmode"

# Run bwa mem for each parental genome
bwa mem -M -R "$RG" "$genome1" "$R1" "$R2" > "${SAMPLE}.par1.sam"
bwa mem -M -R "$RG" "$genome2" "$R1" "$R2" > "${SAMPLE}.par2.sam"
bwa mem -M -R "$RG" "$genome3" "$R1" "$R2" > "${SAMPLE}.par3.sam"

echo "Mapping complete for ${SAMPLE}"

done

```

It takes ~12 hours to run. So if I'm running 32 samples at a time (max allowed by HPC, also don't queue more than 100 at a time), that would be ~30 samples a day, and it would take ~10 days to do the alignment.

## 2. Joint filtering of genomes

Identify reads that do not map uniquely to parental genome and exclude them. Input is three .sam files for each sample (one per species reference), output is three sorted.pass.unique.bam for each sample. All of this takes place in each samples directory in TMPDIR, total array # is 308.

This is the script ```run_samtools_to_jointfiltering_DT.sh```:

```
#!/bin/bash
#SBATCH --job-name=samtools_filter
#SBATCH --output=/project/dtataru/ancestryinfer/logs/samtools_filter__%A_%a.out
#SBATCH --error=/project/dtataru/ancestryinfer/logs/samtools_filter__%A_%a.err
#SBATCH --time=3-00:00:00
#SBATCH -p single
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task=12
#SBATCH -A loni_ferrislac
#SBATCH --array=2  # For some reason this starts at array 2

### LOAD MODULES ###
#module load samtools/1.19
#module load bcftools/1.18
eval "$(conda shell.bash hook)"
conda activate /home/dtataru/.conda/envs/ancestryinfer

### ASSIGN VARIABLES ###
P=$(find /work/dtataru/TMPDIR/ -type d | sort | awk -v line=${SLURM_ARRAY_TASK_ID} 'line==NR')
SAMPLE=$(echo $P | cut -d "/" -f 5 | cut -d "_" -f 1)

genome1="/project/dtataru/hybrids/ancestryinfer/reference_genomes/MguttatusTOL_551_v5.0.fa"
genome2="/project/dtataru/hybrids/ancestryinfer/reference_genomes/Mnasutusvar_SF_822_v2.0.fa"
genome3="/project/dtataru/hybrids/ancestryinfer/reference_genomes/WLF47.fasta"

TMPDIR="/work/dtataru/TMPDIR/${SAMPLE}"
cd "$TMPDIR"

echo "Working in TMPDIR: $TMPDIR"
echo "Sample: $SAMPLE"

### PARAMETERS ###
QUALITY=29
MAX_ALIGN=0 #was set to 2000000 before, but trying out keeping all reads
RATE=0.005
THREADS=20

### MERGE LANES ###
echo "merging lanes"

for P in 1 2 3; do
    SAM_FILES=("${TMPDIR}/${SAMPLE}"*_L00?.par${P}.sam)
	 BAM_FILES=()

   for SAM in "${SAM_FILES[@]}"; do
        BAM="${SAM%.sam}.bam"
        echo "Converting $SAM to $BAM"
        samtools fixmate -O bam "$SAM" "$BAM"
        BAM_FILES+=("$BAM")
    done

	MERGED_BAM="${TMPDIR}/${SAMPLE}.merged.par${P}.bam"
    samtools merge -f -@ ${THREADS} "$MERGED_BAM" "${BAM_FILES[@]}"
done

### SAMTOOLS SORT AND FILTER ###
echo "Starting Samtools"

for P in 1 2 3; do
    SORTED="${TMPDIR}/${SAMPLE}.par${P}.sorted.bam"
    UNIQUE="${TMPDIR}/${SAMPLE}.par${P}.sorted.unique.bam"

    samtools sort -@ 12 -o $SORTED $MERGED_BAM
    samtools index $SORTED
    samtools view -b -q $QUALITY $SORTED > $UNIQUE
    samtools index $UNIQUE
done

echo "Samtools Done"

### JOINT FILTERING ###

echo "Start Joint Filtering"

samtools view -F 4 ${SAMPLE}.par1.sorted.unique.bam | cut -f1 > p1_pass
samtools view -F 4 ${SAMPLE}.par2.sorted.unique.bam | cut -f1 > p2_pass
samtools view -F 4 ${SAMPLE}.par3.sorted.unique.bam | cut -f1 > p3_pass

comm -12 <(sort p1_pass) <(sort p2_pass) | \
comm -12 - <(sort p3_pass) > pass_all

if [[ $MAX_ALIGN -gt 0 ]]; then
    echo "Subsampling to $MAX_ALIGN alignments"
    shuf -n $MAX_ALIGN pass_all -o pass_all
fi

for P in 1 2 3; do
    samtools view -N pass_all -b ${SAMPLE}.par${P}.sorted.unique.bam > ${SAMPLE}.par${P}.sorted.pass.unique.bam
	samtools index ${SAMPLE}.par${P}.sorted.pass.unique.bam
done

echo "Job Done"

```

### 3. Variant Calling 

Joint variant calling across all samples with bcftools, then genotype calls matching parent 1 (coordinate space) alleles at ancestry informative sites are counted from a joint samtools mpileup file. Note, this script calls ```vcf_to_counts_non-colinear_DTv3.pl```. Input files are each of the ${SAMPLE}.par${P}.sorted.pass.unique.bam files, which are merged, and output files are hybrid1 vcf counts and .bed files. Initially, I looped this through all parent files but I think it is only needed for par1. I have kept the format which can easily incorporate all parents if needed (just add for P in 1 2 3). The script will fail if there are empty or corrupt bams. Somehow, everytime I have run Step 2 empty hidden files have been created in the directory, so before running remove those:

```
rm /work/dtataru/TMPDIR/.par1.sorted.pass.unique.bam
rm /work/dtataru/TMPDIR/.par2.sorted.pass.unique.bam
rm /work/dtataru/TMPDIR/.par3.sorted.pass.unique.bam
```
Also, there were three "Bad bams", GBG10, HHH46, and SHG11, which I move to the BAD_BAMS folder to be able to merge the other bams. Then run the script called ```varcall_readcount_DT.sh```:

```
#!/bin/bash
#SBATCH --job-name=varcall
#SBATCH --output=/project/dtataru/ancestryinfer/logs/varcall_%j.out
#SBATCH --error=/project/dtataru/ancestryinfer/logs/varcall__%j.err
#SBATCH --time=3-00:00:00
#SBATCH -p single
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task=12
#SBATCH -A loni_ferrislac

### LOAD MODULES ###
module load samtools/1.19
module load bcftools/1.18
eval "$(conda shell.bash hook)"
conda activate /home/dtataru/.conda/envs/ancestryinfer

### ASSIGN VARIABLES ###
genome1="/project/dtataru/hybrids/ancestryinfer/reference_genomes/MguttatusTOL_551_v5.0.fa"
genome2="/project/dtataru/hybrids/ancestryinfer/reference_genomes/Mnasutusvar_SF_822_v2.0.fa"
genome3="/project/dtataru/hybrids/ancestryinfer/reference_genomes/WLF47.fasta"

P=$(find /work/dtataru/TMPDIR/ -type d | sort | awk -v line=${SLURM_ARRAY_TASK_ID} 'line==NR')
SAMPLE=$(echo $P | cut -d "/" -f 6 | cut -d "." -f 1)

PATH_SCRIPTS="/project/dtataru/ancestryinfer"
AIMS="/project/dtataru/ancestryinfer/AIMs_panel15_final.AIMs.txt"
AIM_COUNTS="/project/dtataru/ancestryinfer/AIMs_panel15_final.AIMs_counts.txt"
WORKDIR="/work/dtataru/HYBRIDS/HMM_INPUT"
BAMDIR="/work/dtataru/TMPDIR/"
THREADS=20

### CHECK BAM FILES FOR CORRUPTION ###
#script will fail if there are corrupt bams
#echo "Checking BAM files for corruption or emptiness..."
#CORRUPT_COUNT=0

#for f in /work/dtataru/TMPDIR/*/*.sorted.pass.unique.bam; do
#    if ! samtools quickcheck -v "$f"; then
#        echo "Corrupt or empty file detected: $f"
#        CORRUPT_COUNT=$((CORRUPT_COUNT+1))
#    fi
#done

#if [ $CORRUPT_COUNT -gt 0 ]; then
#    echo "ERROR: Found $CORRUPT_COUNT corrupt or empty BAM files. Exiting job."
#    exit 1
#else
#    echo "All BAM files passed samtools quickcheck."
#fi

### MERGE ALL BAMS FOR VARIANT CALLING ###
#echo "Merge BAM files"
#cd "$BAMDIR"

#for P in 1; do
#	BAM_FILES=($(find "$BAMDIR" -type f -name "*.par${P}.sorted.pass.unique.bam" | sort))
#	MERGED="${WORKDIR}/hybrids1merged.par${P}.pass.unique.bam"
#	SORTED="${WORKDIR}/hybrids1merged.par${P}.sorted.pass.unique.bam"

#   	samtools merge -r -c -p -@ ${THREADS} "$MERGED" "${BAM_FILES[@]}"
#	samtools sort -@ 12 -o "$SORTED" "$MERGED"
#	samtools index "$SORTED"
#done
#echo "BAM files merged"

```
In this script, the I changed the mpileup step from this to what is current, because default max-depth is 250 and that is way too low for 308 samples at 23x. It timed out after 3 days of trying to make it, though. Now I am going to run a second script for variant calling that runs an array by chromosome called ```varcall_simple.sh```:

```
#!/bin/bash
#SBATCH --job-name=varcall
#SBATCH --output=/project/dtataru/ancestryinfer/logs/varcall_%j.out
#SBATCH --error=/project/dtataru/ancestryinfer/logs/varcall_%j.err
#SBATCH --time=3-00:00:00
#SBATCH -p single
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task=24
#SBATCH -A loni_ferrislac
#SBATCH --array=1-14       
#SBATCH --mem=32G

module load bcftools
module load samtools

### ASSIGN VARIABLES ###
PATH_SCRIPTS="/project/dtataru/ancestryinfer"
genome1="/project/dtataru/hybrids/ancestryinfer/reference_genomes/MguttatusTOL_551_v5.0.fa"
WORKDIR="/work/dtataru/HYBRIDS/HMM_INPUT"
BAM_FILE="${WORKDIR}/hybrids1merged.par1.sorted.pass.unique.bam"
FOCAL_CHROM_LIST="/project/dtataru/BWB/ancestryinfer/focal_chrom_list.txt"
CHR=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$FOCAL_CHROM_LIST")
OUTVCF="hybrids1.par1.maxdepth6000.${CHR}.vcf"
PATH_SCRIPTS="/project/dtataru/ancestryinfer"

### VARIANT CALLING ###
echo "start variant calling"

cd "$WORKDIR"

#exit if vcf already exists
if [[ -s "$OUTVCF" ]]; then
    echo "Output exists, skipping: $OUTVCF"
    exit 0
fi

bcftools mpileup -Ou -d 6000 -r "$CHR" -f "$genome1" "$BAM_FILE" | bcftools call -m -Ov -o "$OUTVCF"

### GENERATE HMM INPUT FILES ###
echo "start generating hmm input"

### AIMS TO COUNTS ###
INFILE_AIMS="${OUTVCF}.aims"
COUNTS="${OUTVCF}_counts"
AIMS="/project/dtataru/ancestryinfer/AIMs_panel15_final.AIMs.txt"
AIM_COUNTS="/project/dtataru/ancestryinfer/AIMs_panel15_final.AIMs_counts.txt"
AIMS_MOD="${AIMS}.mod"
AIMS_BED="${AIMS}.mod.bed"
COUNTS_BED="${COUNTS}.bed"

awk '!/^#/ {print $1"_"$2"\t"$0}' "$OUTVCF" > "$VCF_MOD"
perl "${PATH_SCRIPTS}/combine_FAS_scriptome_snippet.pl" "$AIMS_MOD" "$VCF_MOD" "$INFILE_AIMS"
perl "${PATH_SCRIPTS}/vcf_to_counts_non-colinear_DTv4.pl"  "$INFILE_AIMS" "$COUNTS"
perl -F'_|\t' -lane 'print join("\t", $F[0], $F[1], @F[4..$#F])' "$COUNTS" > "$COUNTS_BED"

### COUNTS TO HMM INPUT ###
perl "${PATH_SCRIPTS}/vcf_counts_to_hmmv3.pl" "$COUNTS_BED" "$AIM_COUNTS" 0.00000002 > "${COUNTS}.hmmsites1"

echo "Job Done"
```

Number of sites for each Chromosome

|Chrom|# of AIMS|
|-----|---------|
|  1  |  8948   |
|  2  |  14386  |
|  3  |  10409  |
|  4  |  17211  |
|  5  |  12362  |
|  6  |  15483  |
|  7  |  8235   |
|  8  |  18949  |
|  9  |  10019  |
|  10 |  14372  |
|  11 |  10866  |
|  12 |  13923  |
|  13 |  14492  |
|  14 |  24575  |

total 194,230 sites

### 4. AncestryHMM

Run main ancestryhmm program. Because ancestryinfer uses read counts instead of genotype calls, there is a step in between these in that wrapper that thin to one AIM per read for each individual, to account for non-independence. Because we are using genotype calls (-g), this is not necessary. Calling this ```run_hmm_DT.sh ```:

```
#!/bin/bash
#SBATCH --job-name=ancestryhmm
#SBATCH --output=/project/dtataru/ancestryinfer/logs/ancestryhmm_%j.out
#SBATCH --error=/project/dtataru/ancestryinfer/logs/ancestryhmm__%j.err
#SBATCH --time=05:00:00
#SBATCH -p workq
#SBATCH -N 1
#SBATCH --cpus-per-task=12
#SBATCH -A loni_ferrislac
#SBATCH --array=1             # one task = one sample

### LOAD MODULES ###
module load python/3.11.5-anaconda
module load boost/1.83.0/intel-2021.5.0
module load gcc/13.2.0
module load gsl/2.7.1/intel-2021.5.0
module load bcftools
eval "$(conda shell.bash hook)"
conda activate /home/dtataru/.conda/envs/ancestryinfer
export PATH=$PATH:/project/dtataru/ancestryinfer/Ancestry_HMM/src
export LD_LIBRARY_PATH=${CONDA_PREFIX}/lib:$LD_LIBRARY_PATH

### ASSIGN VARIABLES ###
PATH_SCRIPTS="/project/dtataru/ancestryinfer"
AIMS="/project/dtataru/ancestryinfer/AIMs_panel15_final.AIMs.txt"
INPUTDIR="/work/dtataru/HYBRIDS/HMM_INPUT"
WORKDIR="/work/dtataru/HYBRIDS/HMM_OUTPUT"
FOCAL_CHROM_LIST="/project/dtataru/BWB/ancestryinfer/focal_chrom_list.txt"
cd $WORKDIR

### SAMPLE INFO ###
#only create this sample list once
#echo "Creating current.samples.list file"
#bcftools query -l hybrids1.par1.maxdepth6000.Chr-01.vcf > current.samples.list

SAMPLE_ID=$SLURM_ARRAY_TASK_ID
SAMPLE_NAME=$(sed -n "${SAMPLE_ID}p" "${INPUTDIR}/current.samples.list" | awk '{print $1}')
echo "Processing SAMPLE: $SAMPLE_NAME    SAMPLE_ID=$SAMPLE_ID"

### RUN HMM ###

### LOOP OVER CHROMS ###
while read CHROM; do

    echo "  → Processing chromosome $CHROM"

    INPUT="${INPUTDIR}/hybrids1.par1.maxdepth6000.${CHROM}.vcf_counts.hmmsites1"

    ### DETERMINE SAMPLE COLUMNS ###
    A_col=$((10 + (SAMPLE_ID-1)*2))
    a_col=$((11 + (SAMPLE_ID-1)*2))

    echo "     Sample columns: A=$A_col   a=$a_col"

    OUTFILE="${WORKDIR}/${SAMPLE_NAME}.${CHROM}.counts.hmmsites1"

    ### Extract first 9 columns + sample's two genotype columns ###
    awk -v A=$A_col -v a=$a_col -v OFS="\t" \
        '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$A,$a}' \
        "$INPUT" > "$OUTFILE"

    ### Write sample-specific sample list file ###
    echo -e "${SAMPLE_NAME}\t2" > ${SAMPLE_NAME}.${CHROM}.samples

    ### RUN ANCESTRY_HMM ###
    ancestry_hmm \
        -i "$OUTFILE" \
        -s "${SAMPLE_NAME}.${CHROM}.samples" \
        -a 3 0.33 0.33 0.34 \
        -p 0 -10000 0.33 \
        -p 1 -1000 0.33 \
        -p 2 -500 0.34 \
        -e 0.05 -g

    echo "Done with $SAMPLE_NAME on $CHROM"

done < "$FOCAL_CHROM_LIST"

echo "ALL COMPLETE for sample $SAMPLE_NAME"

```

## 5. PARSE TSVs FOR PLOTTING
Parsing ancestry from the output tsvs for downstream analysis and visualization. Output of this can be analyzed in something similar to  Banarjee et al. 2023's```Tlalica_three-way_hybrids/local_ancestry_calling/local_ancestry_plots.R```. This is called '''5_ProcessHMM_Output.sh'''

```
#!/bin/bash
#SBATCH --job-name=ProcessHMM_Output
#SBATCH --output=/project/dtataru/ancestryinfer/logs/ProcessHMM_Output_%j.out
#SBATCH --error=/project/dtataru/ancestryinfer/logs/ProcessHMM_Output_%j.err
#SBATCH --time=03:00:00
#SBATCH -p single
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -A loni_ferrislac                      

# Load Modules
module load R/3.2.4

# Set variables
PATH_SCRIPTS="/project/dtataru/ancestryinfer"
INPUT_DIR="/work/dtataru/HYBRIDS/HMM_OUT"
WORKING_DIR="work/dtataru/HYBRIDS/HMM_POSTPROCESS"
POSTERIOR_THRESH=0.8
INTERVALS_RSCRIPT="${PATH_SCRIPTS}/identify_intervals_ancestryinfer_DTv2.R"
FOCAL_CHROM_LIST="/project/dtataru/BWB/ancestryinfer/focal_chrom_list.txt"      
TAG=$(paste -sd "_" "$FOCAL_CHROM_LIST")

cd $WORKING_DIR

### Create Input Files ###

current_list="current.posterior.samples.list_${TAG}"
read_list="current.samples.read.list_${TAG}"

# loop through all .posterior files
for p in *.posterior; do
    sample="${p%.posterior}"
    echo "$sample" >> $current_list

    for hmmsites in ${INPUT_DIR}/${sample}.Chr-*.counts.hmmsites1; do
        if [[ -f "$hmmsites" ]]; then
            echo "$hmmsites" >> "$read_list"
        else
            echo "WARNING: no hmmsites1 files found for sample $sample" >&2
        fi
    done
done

#remove duplicates
sort -u "$current_list" -o "$current_list"
sort -u "$read_list" -o "$read_list"

### CONVERT TO TSV ###

perl ${PATH_SCRIPTS}/convert_rchmm_to_ancestry_tsv_3way_v2_DT.pl current.posterior.samples.list current.samples.read.list 1 ${FOCAL_CHROM_LIST}

echo "Transposing ancestry probabilities"
perl ${PATH_SCRIPTS}/transpose_tsv.pl ancestry-probs-par1_transposed_Chr-01_Chr-02_Chr-03_Chr-04_Chr-05_Chr-06_Chr-07_Chr-08_Chr-09_Chr-10_Chr-11_Chr-12_Chr-13_Chr-14.tsv
perl ${PATH_SCRIPTS}/transpose_tsv.pl ancestry-probs-par2_transposed_Chr-01_Chr-02_Chr-03_Chr-04_Chr-05_Chr-06_Chr-07_Chr-08_Chr-09_Chr-10_Chr-11_Chr-12_Chr-13_Chr-14.tsv
perl ${PATH_SCRIPTS}/transpose_tsv.pl ancestry-probs-par3_transposed_Chr-01_Chr-02_Chr-03_Chr-04_Chr-05_Chr-06_Chr-07_Chr-08_Chr-09_Chr-10_Chr-11_Chr-12_Chr-13_Chr-14.tsv
perl ${PATH_SCRIPTS}/transpose_tsv.pl ancestry-probs-par1par2_transposed_Chr-01_Chr-02_Chr-03_Chr-04_Chr-05_Chr-06_Chr-07_Chr-08_Chr-09_Chr-10_Chr-11_Chr-12_Chr-13_Chr-14.tsv
perl ${PATH_SCRIPTS}/transpose_tsv.pl ancestry-probs-par1par3_transposed_Chr-01_Chr-02_Chr-03_Chr-04_Chr-05_Chr-06_Chr-07_Chr-08_Chr-09_Chr-10_Chr-11_Chr-12_Chr-13_Chr-14.tsv
perl ${PATH_SCRIPTS}/transpose_tsv.pl ancestry-probs-par2par3_transposed_Chr-01_Chr-02_Chr-03_Chr-04_Chr-05_Chr-06_Chr-07_Chr-08_Chr-09_Chr-10_Chr-11_Chr-12_Chr-13_Chr-14.tsv

### parse ancestry in transposed tsvs ###
echo "parse transposed tsv start"
POSTERIOR_THRESH=0.8
perl ${PATH_SCRIPTS}/parse_3way_tsv_to_genotypes_file.pl Chr-01_Chr-02_Chr-03_Chr-04_Chr-05_Chr-06_Chr-07_Chr-08_Chr-09_Chr-10_Chr-11_Chr-12_Chr-13_Chr-14.tsv $POSTERIOR_THRESH > genotypes.txt

perl ${PATH_SCRIPTS}/parse_3way_tsv_ancestry.pl Chr-01_Chr-02_Chr-03_Chr-04_Chr-05_Chr-06_Chr-07_Chr-08_Chr-09_Chr-10_Chr-11_Chr-12_Chr-13_Chr-14.tsv 0.8 > ancestry_proportions.txt

echo "all files parsed"

### identify intervals ###
echo "run intervals R script"
Rscript ${PATH_SCRIPTS}/identify_intervals_ancestryinfer_DTv2.R genotypes.txt ${PATH_SCRIPTS}
echo "identified intervals"
```
