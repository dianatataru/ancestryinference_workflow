# ancestryinference_workflow
Breaking up the steps of ancestryinfer to run highthroughput on large natural complex hybrid dataset

## 1. Aligning all samples to three reference genomes

Use ```bwa mem``` to map reads from hybrid to all parental references independently. This is run using a script named ```map_array_DT.sh```. In the TMPDIR, this creates a folder for each sample, with all reads aligned for each lane run. In the next step, these separate runs will be merged.

```
#!/bin/bash
#SBATCH --output=/project/dtataru/ancestryinfer/logs/map_%A_%a.out
#SBATCH --error=/project/dtataru/ancestryinfer/logs/map_%A_%a.err
#SBATCH --time=7-00:00:00
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

## 2. Joint filtering of genomes

Identify reads that do not map uniquely to parental genome and exclude them using. 

```
perl $path/run_samtools_to_hmm_threegenomes_v2.pl $current_job $genome1 $genome2 $genome3 $read_length $save_files $max_align $focal_chrom $rec_M_per_bp $path $quality\n
```

run_samtools_to_jointfiltering_DT.sh involves:

```
#!/bin/bash
#SBATCH --job-name=samtools_filter
#SBATCH --output=/project/dtataru/ancestryinfer/logs/samtools_filter_%j.out
#SBATCH --error=/project/dtataru/ancestryinfer/logs/samtools_filter_%j.err
#SBATCH --time=3-00:00:00
#SBATCH -p single
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task=12
#SBATCH --mem=64G
#SBATCH -A loni_ferrislab
#SBATCH --array=1-4  # Adjust based on number of samples

### LOAD MODULES ###
module load samtools/1.18
module load bedtools/2.31.1
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
QUALITY=30
MAX_ALIGN=100000
RATE=0.005
AIMS="/project/dtataru/ancestryinfer/AIMs_panel15_final.AIMs.txt"
PATH_SCRIPTS="/project/dtataru/ancestryinfer"
THREADS=20

### MERGE LANES ###
echo "merging lanes"
for L in 1 2 3 4 & P in 1 2 3; do
    SAM="${TMPDIR}/${SAMPLE}*_L00${L}.par${P}.sam"
    BAM="${TMPDIR}/${SAMPLE}*_L00${L}.par${P}.bam"

    samtools fixmate -O bam $SAM $BAM
    samtools merge -f -@ ${THREADS} ${TMPDIR}/${SAMPLE}.merged.par${P}.bam ${BAM}
done

### SAMTOOLS SORT AND FILTER ###
echo "Starting Samtools"

for P in 1 2 3; do
    BAM_MERGED="${TMPDIR}/${SAMPLE}.par${P}.bam"
    SORTED="${TMPDIR}/${SAMPLE}.par${P}.sorted.bam"
    UNIQUE="${TMPDIR}/${SAMPLE}.par${P}.sorted.unique.bam"

    samtools sort -@ 12 -o $SORTED $BAM_MERGED
    samtools index $SORTED
    samtools view -b -q $QUALITY $SORTED > $UNIQUE
    samtools index $UNIQUE
done

echo "Samtools Done"

### JOINT FILTERING ###

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

```

### 3. Variant Calling and Read Counting

Variant calling with bcftools, then reads matching each parental allele at ancestry informative sites are counted from a samtools mpileup file for each hybrid individual. 

At some point I think I need to reformat aims like they did in the wrapper:


```
my $mpileup1 = "$unique1".".bcf";

    if($chrom_string ne 'allchrs'){
    system("bcftools mpileup -r $chrom_string -o $mpileup1 -f $genome1 $finalbam1

### GENERATE HMM INPUT FILES ###

my $aims_bed="$aims".".mod";
system("cat $aims | perl -p -e 's/_/\t/g' | awk -v OFS=\'\\t\' \'\$1=\$1\"\_\"\$2\' > $aims_bed");
my $aims_truebed="$aims".".mod.bed";
system("cat $aims | perl -p -e 's/_/\t/g' | awk -v OFS=\'\\t\' \'\$1=\$1\"\\t\"\$2\' > $aims_truebed");

hybfile="HMM.hybrid.files.list.${SAMPLE}${TAG}"
parfile="HMM.parental.files.list.${SAMPLE}${TAG}"
hybfile=${hybfile//\//} ; hybfile=${hybfile//../.}
parfile=${parfile//\//} ; parfile=${parfile//../.}

echo "${TMPDIR}/${SAMPLE}${TAG}.hmm.pass.formatted" > $hybfile
echo "${TMPDIR}/${SAMPLE}${TAG}.hmm.parental.format" > $parfile

```

### 4. Thinning to one AIM per read across individuals

Counts for each parental allele at ancestry informative sites are subsampled to thin to one ancestry informative site per read if multiple sites occur within one read. This thinning is performed jointly across individuals such that the same site is retained for all individuals in the dataset.

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
