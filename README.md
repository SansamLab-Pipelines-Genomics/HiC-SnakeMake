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
bwa mem -t {params.bwaThreads} -SP5M {params.bwaIndex} {input.fq1} {input.fq2} | samtools view -Shb - > {output.bam}
```
### 2.  make_pairs_with_pairtools_parse
```bash
[ -d {params.tempdir} ] || mkdir {params.tempdir}
samtools view -h {input.bam} | \
 pairtools parse -c {params.chrom_sizes} --add-columns mapq | \
 pairtools sort --nproc {params.nproc} --memory {params.memory} --tmpdir {params.tempdir} --output {output.sorted_pairs}
rm -rf {params.tempdir}
```
### 3.  mark_duplicates_with_pairtools_dedup
```bash
pairtools dedup \
 --mark-dups \
 --output-dups - \
 --output-unmapped - \
 --output {output.marked_pairs} \
 {input.sorted_pairs}
pairix {output.marked_pairs}
```
### 4.  filter_pairs
```bash
## Generate lossless bam
pairtools split \
 --output-sam {output.lossless_bam} \
 {input.marked_pairs}
 
# Select UU, UR, RU reads
pairtools select \
 '(pair_type == "UU") or (pair_type == "UR") or (pair_type == "RU")' \
 --output-rest {output.unmapped_sam} \
 --output {params.temp_file} \
 {input.marked_pairs}
 
pairtools split \
 --output-pairs {params.temp_file1} \
 {params.temp_file}
 
pairtools select 'True' \
 --chrom-subset {params.chrom_sizes} \
 -o {output.dedup_pairs} \
 {params.temp_file1}
 
pairix {output.dedup_pairs}  # sanity check & indexing
```
### 5.  add_frag2Pairs
### 6.  run_cooler

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
-t {cluster.time} \
--mem {cluster.mem} \
--output {cluster.output}'"
```

### 8.  Check results, and when finished, exit environment.
The results will be saved to the "results" folder. Look over log files generated in either the logs/ or logs/snakelogs folders (depending on whether slurm was used).
```bash
conda deactivate
```
