# Process_HiC_SnakeMake
 
## Table of contents:
* [Description of individual steps in pipeline](https://github.com/SansamLab/Process_HiC_SnakeMake/edit/main/README.md#description-of-individual-steps-in-pipeline)
  * [1.  run_bwa_mem](https://github.com/SansamLab/Process_HiC_SnakeMake/blob/main/README.md#1--run_bwa_mem)
  * [2.  make_pairs_with_pairtools_parse](https://github.com/SansamLab/Process_HiC_SnakeMake/blob/main/README.md#2--make_pairs_with_pairtools_parse)
  * [3.  mark_duplicates_with_pairtools_dedup](https://github.com/SansamLab/Process_HiC_SnakeMake/blob/main/README.md#3--mark_duplicates_with_pairtools_dedup)
  * [4.  filter_pairs](https://github.com/SansamLab/Process_HiC_SnakeMake/blob/main/README.md#4--filter_pairs)
  * [5.  add_frag2Pairs](https://github.com/SansamLab/Process_HiC_SnakeMake/blob/main/README.md#5--add_frag2pairs)
  * [6.  run_cooler](https://github.com/SansamLab/Process_HiC_SnakeMake/blob/main/README.md#6--run_cooler)
* [Step-by-step instructions on running Snakemake pipeline:](https://github.com/SansamLab/Process_HiC_SnakeMake/blob/main/README.md#step-by-step-instructions-on-running-snakemake-pipeline)
  * [1.  Load slurm and miniconda](https://github.com/SansamLab/Process_HiC_SnakeMake/blob/main/README.md#1--load-slurm-and-miniconda)
  * [2.  Clone repository](https://github.com/SansamLab/Process_HiC_SnakeMake#2--clone-repository)
  * [3.  Start the conda environment](https://github.com/SansamLab/Process_HiC_SnakeMake/blob/main/README.md#3--start-the-conda-environment)
    * [3A.  FIRST TIME ONLY:  Setup conda environment](https://github.com/SansamLab/Process_HiC_SnakeMake/blob/main/README.md#3a--first-time-only--setup-conda-environment)
    * [3B.  Activate conda environment](https://github.com/SansamLab/Process_HiC_SnakeMake/blob/main/README.md#3b--activate-conda-environment)
  * [4.  Modify the job-specific configuration files.](https://github.com/SansamLab/Process_HiC_SnakeMake/blob/main/README.md#4--modify-the-job-specific-coniguration-files)
    * [4A.  Modify the config/config.yml file](https://github.com/SansamLab/Process_HiC_SnakeMake/blob/main/README.md#4a--modify-the-configconfigyml-file)
    * [4B.  Modify the config/samples.csv file](https://github.com/SansamLab/Process_HiC_SnakeMake/blob/main/README.md#4b--modify-the-configsamplescsv-file)
    * [4C.  IF SLURM RESOURCE CHANGES ARE NEEDED. Modify the config/cluster_config.yml file](https://github.com/SansamLab/Process_HiC_SnakeMake/blob/main/README.md#4c--if-slurm-resource-changes-are-needed-modify-the-configcluster_configyml-file)
  * [5.  Do a dry run](https://github.com/SansamLab/Process_HiC_SnakeMake/blob/main/README.md#4--do-a-dry-run)
  * [6.  Make a DAG diagram](https://github.com/SansamLab/Process_HiC_SnakeMake/blob/main/README.md#5--make-a-dag-diagram)
  * [7.  Run on cluster with slurm](https://github.com/SansamLab/Process_HiC_SnakeMake/blob/main/README.md#6--run-on-cluster-with-slurm)

## Description of individual steps in pipeline:
![DAG of Pipeline](dag.svg)

### 1.  run_bwa_mem
```bash
bwa mem \
 -t {params.bwaThreads} \
 -SP5M \
 {params.bwaIndex} \
 {input.fq1} \
 {input.fq2} | \
 samtools view -Shb - > {output.bam}
```
### 2.  make_pairs_with_pairtools_parse
```bash
# create a temporary directory for sorting the .pairs file
[ -d {params.tempdir} ] || mkdir {params.tempdir}

# Find ligation junctions in .sam, make sorted .pairsam file. 
samtools view -h {input.bam} | \ `# convert .bam file from bwa-mem alignment to .sam`
 pairtools parse \`               # parse ligated reads into pairs`
  -c {params.chrom_sizes} \`       ## path to file with chromosome sizes`
  --add-columns mapq | \`          ## option to add column with read mapq scores`
 pairtools sort \`                # sort pairs`
  --nproc {params.nproc} \`        ## number of processors to use for sorting`
  --memory {params.memory} \`      ## amount of memory for sorting`
  --tmpdir {params.tempdir} \`     ## temporary directory path for sorting`
  --output {output.sorted_pairs}`  ## path and name of .pairs file made`

# Delete the temporary directory used for sorting
rm -rf {params.tempdir}
```
### 3.  mark_duplicates_with_pairtools_dedup
```bash
# Find and mark PCR/optical duplicates in .pairsam file from step 2.
pairtools dedup \
 --mark-dups \`                    # duplicate pairs are marked as DD in pair_type and as a duplicate in the sam entries.`
 --output-dups - \`                # output duplicates together with deduped pairs`
 --output-unmapped - \`            # output unmapped pairs together with deduped pairs`
 --output {output.marked_pairs} \` # name of output .pairs file with duplicates marked`
 {input.sorted_pairs}`             # .pairsam file input from step 2`
 
# index the .pairsam
pairix {output.marked_pairs}
```
### 4.  filter_pairs
```bash
## Generate lossless bam from the pairsam file
pairtools split \
 --output-sam {output.lossless_bam} \`# name of .bam file produced`
 {input.marked_pairs}`                # .pairsam file input from step 3`
 
# Select UU, UR, RU reads
 ## UU = unique-unique, both alignments are unique
 ## UR or RU = unique-rescued, one alignment unique, the other rescued
pairtools select \
 '(pair_type == "UU") or (pair_type == "UR") or (pair_type == "RU")' \
 --output-rest {output.unmapped_sam} \ `# name of file with the remainder of read pairs`
 --output {params.temp_file} \`         # temporary output file with the selected read pairs`
 {input.marked_pairs}`                  # .pairsam file input from step 3`
 
# Generate .pairs file from the UU, UR, and RU reads selected above
pairtools split \
 --output-pairs {params.temp_file1} \`  # temporary .pairs output file`
 {params.temp_file}`                    # input .pairsam file generated with pairtools select above`
 
# Make a .pairs file with only pairs in chromosomes of interest
pairtools select 'True' \
 --chrom-subset {params.chrom_sizes} \` # path to a chromosomes file containing a chromosome subset of interest.`
 -o {output.dedup_pairs} \`             # ouput path and filename for .pairs file in chromosomes of interest`
 {params.temp_file1}`                   # input .pairsam file generated with pairtools split above`

# index the .pairs file
pairix {output.dedup_pairs}            
```
### 5.  add_frag2Pairs

```bash
# convert to fragment map
gunzip -ck {input.dedup_pairs} | \`   # unzip dsthe .pairs file generated in step 4`
 workflow/scripts/fragment_4dnpairs.pl \
  -a - `                              # allows replacing existing frag1/frag2 columns`\
  {params.frag2_pairs_basename} \ `   # filename for .pairs file generated`
  {params.restriction_file}`          # a restriction site file, which lists on each line, the sorted locations of the enzyme restriction sites.`

# Compress the .pairs file generated
bgzip -f {params.frag2_pairs_basename}

# Index the .pairs file generated
pairix -f {output.frag2_pairs}
```

### 6.  run_cooler

```bash
# use for all chromosomes and contigs
cp {params.chrom_sizes} {params.cooler_tempchrsize}
            
# the cload command requires the chrom size file to exist besides the chrom size bin file.
cooler cload pairix \
 -p {params.cooler_n_cores} \
 -s {params.cooler_max_split} \
 {params.cooler_tempchrsize}:{params.cooler_bin_size} \
 {input.frag2_pairs} \
 {output.cooler}
```

## Step-by-step instructions on running Snakemake pipeline:

### 1.  Load slurm and miniconda
Note. The commands to do this will be different on your machine. These commands are specific to an HPC using slurm with these modules installed.

```bash
ml slurm
ml miniconda
```
### 2.  Clone repository
```bash
git clone https://github.com/SansamLab/Process_HiC_SnakeMake.git
# rename folder with project name
mv Process_HiC_SnakeMake/ My_HiC_Project_Folder/
# change directory into root of your project folder
cd My_HiC_Project_Folder
```
### 3.  Start the conda environment
### 3A.  FIRST TIME ONLY:  Setup conda environment
```bash
# -f is the location of the environment .yml file. 
## The relative path assumes that you are in the root directory of this repository.
# -p is the path where you want to install this environment
conda env create -f workflow/envs/HiCSnakemake.yml -p /s/sansam-lab/HiC_Conda_Environment 
```

### 3B.  Activate conda environment
```bash
conda activate /s/sansam-lab/HiC_Conda_Environment
```

### 4.  Modify the job-specific coniguration files.
#### 4A.  Modify the config/config.yml file

You must enter paths to the following:
* bwa_genome:
  * location of bwa indexed genome for the alignment
* chrom_sizes
  * chromosome sizes file
* juicer_RE_file
  * restriction enzyme file generated with juicer

#### 4B.  Modify the config/samples.csv file

The samples.csv file in the config folder has paths to the test fastq files. You must replace those paths with those for your own fastq files. The first column of each row is the sample name. This name will be used for all output files.

#### 4C.  IF SLURM RESOURCE CHANGES ARE NEEDED. Modify the config/cluster_config.yml file

CPU and memory requests for each rule in the pipeline are detailed in this file. If you are using SLURM, you may need to alter this file to fit your needs/system.

### 5.  Do a dry run.
A dry run produces a text output showing exactly what commands will be executed. Look this over carefully before submitting the full job. It is normal to see warnings about changes made to the code, input, and params.
```bash
snakemake -npr
```

### 6.  Make a DAG diagram.
```bash
snakemake --dag | dot -Tpdf > dag.pdf
```

### 7.  Run on cluster with slurm.
This snakemake pipeline could be executed without slurm, but if an hpc with slurm is used, the following will start the pipeline with the parameters defined in the config/cluster_config.yml file.
```bash
sbatch --wrap="\
snakemake \
-R \
-j 999 \
--cluster-config config/cluster_config.yml \
--cluster '\
sbatch \
-A {cluster.account} \
-p {cluster.partition} \
--cpus-per-task {cluster.cpus-per-task}  \
--mem {cluster.mem} \
--output {cluster.output}'"
```

### 8.  Check results, and when finished, exit environment.
The results will be saved to the "results" folder. Look over log files generated in either the logs/ or logs/snakelogs folders (depending on whether slurm was used).
```bash
conda deactivate
```
