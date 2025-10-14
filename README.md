# ancestryinference_workflow
Breaking up the steps of ancestryinfer to run highthroughput on large natural complex hybrid dataset

## 1. Aligning all samples to three reference genomes

Use ```bwa mem``` to map reads from hybrid to all parental references independently.

```
perl run_map_threegenomes_v1.pl $current_job $genome1 $genome2 $genome3 $read_type $tag\n
 ```

and this would be run_map_threegenomes_v1.pl:

```
       bwa mem -M -R $RG1 $genome1 $read1 $read2 > $sam1
       bwa mem -M -R $RG1 $genome2 $read1 $read2 > $sam2
       bwa mem -M -R $RG1 $genome3 $read1 $read2 > $sam3

#filter for mapping quality>29 and proper pairing with samtools

#remove duplicates with picard

```
## 2. Joint filtering of genomes

Identify reads that do not map uniquely to parental genome and exclude them using ```ngsutils```. 

```
perl $path/run_samtools_to_hmm_threegenomes_v2.pl $current_job $genome1 $genome2 $genome3 $read_length $save_files $max_align $focal_chrom $rec_M_per_bp $path $quality\n
```

Includes:

```
#For each parent

samtools fixmate
samtools sort
samtools index

#Joint filtering

my $par1_pass="$unique1"."_par1_passlist";
    my $par2_pass="$unique2"."_par2_passlist";
    my $par3_pass="$unique3"."_par3_passlist";

    $par1_pass=~ s/\//_/g;
    $par2_pass=~ s/\//_/g;
    $par3_pass=~ s/\//_/g;

    $par1_pass=~ s/_read_1.fastq.gz.par1.sam.sorted.unique.bam//g;
    $par2_pass=~ s/_read_1.fastq.gz.par2.sam.sorted.unique.bam//g;
    $par3_pass=~ s/_read_1.fastq.gz.par3.sam.sorted.unique.bam//g;

    system("samtools view -F 4 $unique1 | cut -f 1 > $par1_pass");
    system("samtools view -F 4 $unique2 | cut -f 1 > $par2_pass");
    system("samtools view -F 4 $unique3 | cut -f 1 > $par3_pass");
 
    my $pass_both="$par1_pass"."_both";
    $pass_both =~ s/_par1//g;
    my $pass_all="$par1_pass"."_all";
    $pass_all =~ s/_par1//g;

 ###intersect x 2
    my $file1 = $par1_pass;
    my $file2 = $par2_pass;
    open F2, $file2 or die $!;
    open JOINT, ">$pass_both";
    while (<F2>) { $h2{$_}++ };
    open F1, $file1 or die;
    $total=$.; $printed=0;
    while (<F1>) { $total++; if ($h2{$_}) { print JOINT $_; $h2{$_} = ""; $printed++; } }

    my $file3 = $par3_pass;
    my $file4 = $pass_both;
    open F3, $file3 or die $!;
    open JOINT2, ">$pass_all";
    while (<F3>) { $h3{$_}++ };
    open F4, $file4 or die;
    $total2=$.; $printed2=0;
    while (<F4>) { $total2++; if ($h3{$_}) { print JOINT2 $_; $h3{$_} = ""; $printed2++; } }

    ###filter
    my $finalbam1="$unique1";
    $finalbam1=~ s/sorted.unique/sorted.pass.unique/g;
    
    my $finalbam2="$unique2";
    $finalbam2=~s/sorted.unique/sorted.pass.unique/g;

    my $finalbam3="$unique3";
    $finalbam3=~s/sorted.unique/sorted.pass.unique/g;

    ###Adding max read filter here
    if($max_align > 0){
    print "limiting to $max_align alignments\n";
    my $tmp_list="$pass_all".".tmp";
    system("shuf -n $max_align $pass_all > $tmp_list");
    system("mv $tmp_list $pass_all");
    }#then subsample, otherwise leave pass list as is
    
    system("samtools view -N $pass_all -b $unique1 > $finalbam1");
    system("samtools view -N $pass_all -b $unique2 > $finalbam2");
    system("samtools view -N $pass_all -b $unique3 > $finalbam3");

    system("samtools index $finalbam1");
    system("samtools index $finalbam2");
    system("samtools index $finalbam3");

#    system("rm $par1_pass $par2_pass $par3_pass $pass_both");


```

### 3. Variant Calling and Read Counting

Variant calling with bcftools, then reads matching each parental allele at ancestry informative sites are counted from a samtools mpileup file for each hybrid individual. 

```
my $mpileup1 = "$unique1".".bcf";

    if($chrom_string ne 'allchrs'){
    system("bcftools mpileup -r $chrom_string -o $mpileup1 -f $genome1 $finalbam1
```

### 4. Thinning to one AIM per read across individuals

Counts for each parental allele at ancestry informative sites are subsampled to thin to one ancestry informative site per read if multiple sites occur within one read. This thinning is performed jointly across individuals such that the same site is retained for all individuals in the dataset.

### 5. AncestryHMM

Requires ```armadillo```.

## 6. Reformatting HMM Output

```
perl $path/convert_rchmm_to_ancestry_tsv_v3.pl current.samples.list current.samples.read.list $save_files $focal_chrom\n"
perl $path/transpose_tsv.pl $final_file1\n"
