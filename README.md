# ancestryinference_workflow
Breaking up the steps of ancestryinfer to run highthroughput on large natural complex hybrid dataset

Check that all files transferred by running ```bash check_missing_fastq.sh /project/dtataru/hybrids/1_hybrid1data/``` in hybrids director/

## 1. Aligning all samples to three reference genomes

Use ```bwa mem``` to map reads from hybrid to all parental references independently. Input files are fastq files for all samples, four lanes per sample and two reads per lane (Total 1,232 arrays).  In the TMPDIR, this creates a folder for each sample, with all reads aligned for each lane run. Output is three .sam files for each sample, corresponding to each species' reference genome. In the next step, these separate runs will be merged. This script is ```map_array_DT.sh ```.

```
#!/bin/bash
#SBATCH --output=/project/dtataru/ancestryinfer/logs/map_%A_%a.out
#SBATCH --error=/project/dtataru/ancestryinfer/logs/map_%A_%a.err
#SBATCH --time=1-00:00:00
#SBATCH -p single
#SBATCH -N 1            #: Number of Nodes
#SBATCH -n 1            #: Number of Tasks per Node
#SBATCH -A loni_ferrislab
#SBATCH --cpus-per-task=12
#SBATCH --array=5-10   # Job array (1-n) when n is number of unique samples that came off the sequencer

### LOAD MODULES ###
#For this step, bwa and needed
eval "$(conda shell.bash hook)"
conda activate /home/dtataru/.conda/envs/ancestryinfer

echo "Start Job"
echo "SLURM_ARRAY_TASK_ID = ${SLURM_ARRAY_TASK_ID}"

### ASSIGN VARIABLES ARRAY SIZE > 1000 ###

R1_FILE_LIST=$(ls -1 /project/dtataru/hybrids/1_hybrid1data/*R1_001.fastq.gz)
R2_FILE_LIST=$(ls -1 /project/dtataru/hybrids/1_hybrid1data/*R2_001.fastq.gz)
# check for matching read files
R1_FILE_LIST_SIZE=$(wc -w <<<$R1_FILE_LIST)
R2_FILE_LIST_SIZE=$(wc -w <<<$R2_FILE_LIST)
for v in R1_FILE_LIST_SIZE R2_FILE_LIST_SIZE; do echo "Using $v=${!v}"; done
if [ $R1_FILE_LIST_SIZE -ne $R2_FILE_LIST_SIZE ]; then
	echo "ERROR: Encountered unmatched read files: R1_FILE_LIST_SIZE=$R1_FILE_LIST_SIZE and R2_FILE_LIST_SIZE=$R2_FILE_LIST._SIZE Exiting."
	exit 1
fi

#NUMBER_FILE_PAIRS_ALREADY_DONE=1000 
NUMBER_FILE_PAIRS_ALREADY_DONE=0
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

R1=$(cut -f $FILE_INDEX -d " "<<<$R1_FILE_LIST)
R2=$(cut -f $FILE_INDEX -d " "<<<$R2_FILE_LIST)

### ASSIGN VARIABLES ARRAY SIZE < 1000 ###
#R1=$(find /project/dtataru/hybrids/1_hybrid1data/ \
#    | grep R1_001.fastq.gz \
#    | sort \
#    | awk -v line=${SLURM_ARRAY_TASK_ID} 'line==NR')
#R2=$(find /project/dtataru/hybrids/1_hybrid1data/ \
#    | grep R2_001.fastq.gz \
#    | sort \
#    | awk -v line=${SLURM_ARRAY_TASK_ID} 'line==NR')

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

```

Started at 2:25 pm, estimated SUs:
sbatch: 129499.42 SUs available in loni_ferrislab
sbatch: 2016.00 SUs estimated for this job.
sbatch: lua: Submitted job 379027

It takes ~3 hours to run, so far. So if I'm running 10 samples at a time, that would be ~30 samples a day, and it would take ~10 days to do the alignment.

## 2. Joint filtering of genomes

Identify reads that do not map uniquely to parental genome and exclude them. Input is three .sam files for each sample (one per species reference), output is three sorted.pass.unique.bam for each sample. All of this takes place in each samples directory in TMPDIR, total array # is 308.

```
perl $path/run_samtools_to_hmm_threegenomes_v2.pl $current_job $genome1 $genome2 $genome3 $read_length $save_files $max_align $focal_chrom $rec_M_per_bp $path $quality\n
```

run_samtools_to_jointfiltering_DT.sh involves:

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
#SBATCH -A loni_ferrislab
#SBATCH --array=2  # For some reason this starts at array 2

### LOAD MODULES ###
module load samtools/1.19
module load bcftools/1.18
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

### 3. Variant Calling and Read Counting

Joint variant calling across all samples with bcftools, then reads matching each parental allele at ancestry informative sites are counted from a samtools mpileup file for each hybrid individual. Note, this script calls ```vcf_to_counts_non-colinear_DTv2.pl```. Input files are each of the ${SAMPLE}.par${P}.sorted.pass.unique.bam files, which are merged by parental species' genome, and output files are hybrid1 vcf counts and .bed files for each parent (3 total). This script will be called ```varcall_readcount_DT.sh```.

```
#!/bin/bash
#SBATCH --job-name=varcall
#SBATCH --output=/project/dtataru/ancestryinfer/logs/varcall_%A_%a.out
#SBATCH --error=/project/dtataru/ancestryinfer/logs/varcall__%A_%a.err
#SBATCH --time=3-00:00:00
#SBATCH -p single
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task=12
#SBATCH -A loni_ferrislab

### LOAD MODULES ###
module load samtools/1.19
module load bcftools/1.18
eval "$(conda shell.bash hook)"
conda activate /home/dtataru/.conda/envs/ancestryinfer

### ASSIGN VARIABLES ###
genome1="/project/dtataru/hybrids/ancestryinfer/reference_genomes/MguttatusTOL_551_v5.0.fa"
genome2="/project/dtataru/hybrids/ancestryinfer/reference_genomes/Mnasutusvar_SF_822_v2.0.fa"
genome3="/project/dtataru/hybrids/ancestryinfer/reference_genomes/WLF47.fasta"

PATH_SCRIPTS="/project/dtataru/ancestryinfer"
AIMS="/project/dtataru/ancestryinfer/AIMs_panel15_final.AIMs.txt"
WORKDIR="/work/dtataru/HMM_INPUT"
BAMDIR="/work/dtataru/TMPDIR"
THREADS=20

### MERGE ALL BAMS FOR VARIANT CALLING ###
cd "$BAMDIR"

for P in 1 2 3; do
    BAM_FILES=("${BAMDIR}/*.par${P}.sorted.pass.unique.bam)
	MERGED="${WORKDIR}/hybrids1merged.par${P}.pass.unique.bam"
	SORTED="${WORKDIR}/hybrids1merged.par${P}.sorted.pass.unique.bam"

   for BAM in "${BAM_FILES[@]}"; do
        samtools merge -r -c -p -@ ${THREADS} $MERGED "${BAM_FILES[@]}
    done

	samtools sort -@ 12 -o $SORTED $MERGED
	samtools index $SORTED
done

### VARIANT CALLING ###
echo "start variant calling"

cd "$WORKDIR"

for P in 1 2 3; do
    GENOME_VAR="genome${P}"
    GENOME=${!GENOME_VAR}
    BCF="hybrids1.par${P}.bcf"
    VCF="hybrids1.par${P}.vcf.gz"

    bcftools mpileup -r -f $GENOME -o $BCF $SORTED
    bcftools call -mO z -o $VCF $BCF
    gunzip $VCF
done

echo "finished variant calling"

### GENERATE HMM INPUT FILES ###
echo "start generating hmm input"

for P in 1 2 3; do
    VCF="hybrids1.par${P}.vcf"
    COUNTS="${VCF}_counts"
    perl $PATH_SCRIPTS/vcf_to_counts_non-colinear_DTv2.pl $VCF $AIMS $PATH_SCRIPTS
    cat $COUNTS | perl -p -e 's/_/\t/g' | \
        awk -v OFS='\t' '{print $1, $2-1, $2, $4, $5, $6}' > ${COUNTS}.bed
done

echo "Job Done"
```

### 4. Thinning to one AIM per read across individuals

Counts for each parental allele at ancestry informative sites are subsampled to thin to one ancestry informative site per read if multiple sites occur within one read. This thinning is performed jointly across individuals such that the same site is retained for all individuals in the dataset. This I will call ```counts_to_hmm_DT.sh```:

```
#!/bin/bash
#SBATCH --job-name=counts_to_hmm
#SBATCH --output=/project/dtataru/ancestryinfer/logs/counts_to_hmm_%A_%a.out
#SBATCH --error=/project/dtataru/ancestryinfer/logs/counts_to_hmm__%A_%a.err
#SBATCH --time=3-00:00:00
#SBATCH -p single
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task=12
#SBATCH -A loni_ferrislab

### LOAD MODULES ###
module load samtools/1.19
module load bcftools/1.18
eval "$(conda shell.bash hook)"
conda activate /home/dtataru/.conda/envs/ancestryinfer

### ASSIGN VARIABLES ###
genome1="/project/dtataru/hybrids/ancestryinfer/reference_genomes/MguttatusTOL_551_v5.0.fa"
genome2="/project/dtataru/hybrids/ancestryinfer/reference_genomes/Mnasutusvar_SF_822_v2.0.fa"
genome3="/project/dtataru/hybrids/ancestryinfer/reference_genomes/WLF47.fasta"

PATH_SCRIPTS="/project/dtataru/ancestryinfer"
AIMS="/project/dtataru/ancestryinfer/AIMs_panel15_final.AIMs.txt"
WORKDIR="/work/dtataru/TMPDIR/HMM_INPUT"

echo "Working in WORKDIR: $WORKDIR"
cd $WORKDIR

for P in 1 2 3; do
    VCF="hybrids1.par${P}.vcf"
    COUNTS="${VCF}_counts"
	perl ${PATH_SCRIPTS}/vcf_counts_to_hmmv2.pl $COUNTS $AIMS 0.00000002 $PATH_SCRIPTS
done

### MAKE INPUT LISTS FOR HMM ###
hybfile="HMM.hybrid.files.list"
parfile="HMM.parental.files.list"

echo "/work/dtataru/TMPDIR/HMM_INPUT/*.hmm.pass.formatted" > $hybfile
echo "/work/dtataru/TMPDIR/HMM_INPUT/*.hmm.parental.format" > $parfile




```

### 5. AncestryHMM

Requires ```armadillo```. 

From main ```Ancestry_HMM_parallel_v7.pl```:

```
my $final_file1="ancestry-probs-par1_transposed"."$tag".".tsv"; my $final_file2="ancestry-probs-par2_transposed"."$tag".".tsv"; my $final_file3="ancestry-probs-par3_transposed"."$tag".".tsv"; 
    my $final_file4="ancestry-probs-par1par2_transposed"."$tag".".tsv";
    my $final_file5="ancestry-probs-par1par3_transposed"."$tag".".tsv";
    my $final_file6="ancestry-probs-par2par3_transposed"."$tag".".tsv";

    print "output files appended with $tag\n";
    open HMMSCRIPT, ">hmm_batch.sh";
    print HMMSCRIPT "$job3_submit\n";
    print HMMSCRIPT "cat $hyb_string > HMM.hybrid.files.list"."$tag"."\n";
    print HMMSCRIPT "cat $par_string >HMM.parental.files.list"."$tag"."\n";

    print HMMSCRIPT "perl $path/combine_all_individuals_hmm_3way_v1.pl HMM.parental.files.list"."$tag HMM.hybrid.files.list"."$tag $minor_prior1 $minor_prior2 $parental_counts_status $initial_admix $initial_admix2 $focal_chrom $read_length $error_prior $tag\n";
#!print "perl combine_all_individuals_hmm_v5.pl HMM.parental.files.list HMM.hybrid.files.list $minor_prior $parental_counts_status $initial_admix $focal_chrom $read_length $error_prior $tag\n";

    print HMMSCRIPT "perl $path/convert_rchmm_to_ancestry_tsv_3way_v1.pl current.samples.list current.samples.read.list $save_files $focal_chrom\n";
    print HMMSCRIPT "perl $path/transpose_tsv.pl $final_file1\n";
    print HMMSCRIPT "perl $path/transpose_tsv.pl $final_file2\n";
    print HMMSCRIPT "perl $path/transpose_tsv.pl $final_file3\n";
    print HMMSCRIPT "perl $path/transpose_tsv.pl $final_file4\n";
    print HMMSCRIPT "perl $path/transpose_tsv.pl $final_file5\n";
    print HMMSCRIPT "perl $path/transpose_tsv.pl $final_file6\n";
    print HMMSCRIPT "rm split_jobs_list\n"; #cleanup split read lists 
    print HMMSCRIPT "rm $aims_truebed $aims_bed\n"; #cleanup different aims files
```

## 6. Reformatting HMM Output

```
perl $path/convert_rchmm_to_ancestry_tsv_v3.pl current.samples.list current.samples.read.list $save_files $focal_chrom\n"
perl $path/transpose_tsv.pl $final_file1\n"
