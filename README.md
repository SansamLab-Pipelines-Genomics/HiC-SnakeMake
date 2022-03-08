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
  * [2.  Clone repository]()
  * [3.  Start the conda environment](https://github.com/SansamLab/Process_HiC_SnakeMake/blob/main/README.md#3--start-the-conda-environment)
    * [3A.  FIRST TIME ONLY:  Setup conda environment](https://github.com/SansamLab/Process_HiC_SnakeMake/blob/main/README.md#3a--first-time-only--setup-conda-environment)
    * [3B.  Activate conda environment](https://github.com/SansamLab/Process_HiC_SnakeMake/blob/main/README.md#3b--activate-conda-environment)
  * [4.  Modify the job-specific coniguration files.](https://github.com/SansamLab/Process_HiC_SnakeMake/blob/main/README.md#4--modify-the-job-specific-coniguration-files)
    * [4A.  Modify the config/config.yml file](https://github.com/SansamLab/Process_HiC_SnakeMake/blob/main/README.md#4a--modify-the-configconfigyml-file)
    * [4B.  Modify the config/samples.csv file](https://github.com/SansamLab/Process_HiC_SnakeMake/blob/main/README.md#4b--modify-the-configsamplescsv-file)
    * [4C.  IF SLURM RESOURCE CHANGES ARE NEEDED. Modify the config/cluster_config.yml file](https://github.com/SansamLab/Process_HiC_SnakeMake/blob/main/README.md#4c--if-slurm-resource-changes-are-needed-modify-the-configcluster_configyml-file)
  * [5.  Do a dry run](https://github.com/SansamLab/Process_HiC_SnakeMake/blob/main/README.md#4--do-a-dry-run)
  * [6.  Make a DAG diagram](https://github.com/SansamLab/Process_HiC_SnakeMake/blob/main/README.md#5--make-a-dag-diagram)
  * [7.  Run on cluster with slurm](https://github.com/SansamLab/Process_HiC_SnakeMake/blob/main/README.md#6--run-on-cluster-with-slurm)

## Description of individual steps in pipeline:
![DAG of Pipeline](dag.svg)

### 1.  run_bwa_mem
### 2.  make_pairs_with_pairtools_parse
### 3.  mark_duplicates_with_pairtools_dedup
### 4.  filter_pairs
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

#### 4B.  Modify the config/samples.csv file

#### 4C.  IF SLURM RESOURCE CHANGES ARE NEEDED. Modify the config/cluster_config.yml file


### 5.  Do a dry run.
```bash
snakemake -npr
```

### 6.  Make a DAG diagram.
```bash
snakemake --dag | dot -Tpdf > dag.pdf
```

### 7.  Run on cluster with slurm.
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

### 8.  When finished, exit environment.
```bash
conda deactivate
```
