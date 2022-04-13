### 1.  Setup directory
#### 1.1 Make main directory
```bash
mkdir McGarvey2022_HiC_Analysis
```
#### 1.2 Get raw .fastq data for HiCs from SRA
This pipeline lacks the functionality to do this. We recommend using the fasterq-dump to transfer the files from SRA. In this example, we create a directory in which the .fastq files are transferred to a subdirectory called "fastqs". The locations of the .fastq files defined in the "Siefert_Samples.csv" provided in this repo as "../fastqs/SRR4036047_1.fastq".
```bash
# load the module with fasterq-dump. note:  this will likely differ on your system
ml ncbi_sra/2.10.7
cd McGarvey2022_HiC_Analysis/
# create a subdirectory to which the .fastq files will be transferred
mkdir fastqs
cd fastqs
# transfer each pair of files for the following SRR entries. 
# In this example we use "sbatch --wrap" to run the fasterq-dump command on our hpc.
SRAIDS=( "SRR12008033" \
"SRR12008034" )
for t in ${SRAIDS[@]}; do
  sbatch --cpus-per-task 12 --wrap="fasterq-dump --split-files --threads 12 $t"
  done
# change back to the McGarvey2022_HiC_Analysis/ directory
cd ..
```
### 1.3 Index genome and prepare restriction enzyme file
This pipeline lacks the funtionality for this. We recommend transferring the genome sequence with wget. At this time, we are uncertain how alternative chromosome sequences will affect the quantitation of aligned sequencing reads. Therefore, we have chosen to use the primary genome assembly, which lacks alternative chromosome sequences. For zebrafish, we separate the primary sequences from alternative sequences using grep.

#### 1.3.1 Load necessary software
```bash
ml seqtk/1.2 bwa/0.7.15 python/3.7.0
```

#### 1.3.2 Copy GRCz11 fasta to local directory
```bash
cd McGarvey2022_HiC_Analysis
mkdir GRCz11
cd GRCz11
sbatch --wrap=\
"wget -e robots=off --reject 'index.html' https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/035/GCF_000002035.6_GRCz11/GCF_000002035.6_GRCz11_genomic.fna.gz"
```

#### 1.3.3. Get primary sequences (non-Alts) and rename chromosomes
```bash
sbatch --wrap=\
"zcat GCF_000002035.6_GRCz11_genomic.fna.gz | seqtk seq -l 3000000000 | grep Primary -A 1 | sed 's/>/>chrUn/g' | seqtk seq -l 100 | sed 's/^.*chromosome />chr/g' | sed 's/,.*$//g' | gzip > GCF_000002035.6_GRCz11_primary_genomic.fna.gz; rm GCF_000002035.6_GRCz11_genomic.fna.gz"
```

#### 1.3.4. Index primary genome with bwa
```bash
ml bwa
sbatch --mem 32G --wrap="bwa index GCF_000002035.6_GRCz11_primary_genomic.fna.gz"
```

#### 1.3.5. Make restriction enzyme file
```bash
# get script from repository
wget https://raw.githubusercontent.com/SansamLab/HiC-SnakeMake/main/workflow/scripts/generate_site_positions.py

# unzip fasta
sbatch --mem 16G --wrap=\
"gunzip --keep GCF_000002035.6_GRCz11_primary_genomic.fna.gz"

# run script
sbatch --mem 24G --wrap=\
"python generate_site_positions.py \
'HindIII' \
'GRCz11' \
'GCF_000002035.6_GRCz11_primary_genomic.fna'"
```

#### 1.3.6. Make chromosome sizes file from genome fasta
```bash
zcat GCF_000002035.6_GRCz11_primary_genomic.fna.gz | awk '$0 ~ ">" {if (NR > 1) {print c;} c=0;printf substr($0,2,100) "\t"; } $0 !~ ">" {c+=length($0);} END { print c; }' > GCF_000002035.6_GRCz11_primary_genomic.sizes
```

### 2.  Transfer .fastq files from SRA.
This pipeline lacks the functionality to do this. We recommend using the fasterq-dump to transfer the files from SRA. In this example, we create a directory in which the .fastq files are transferred to a subdirectory called "fastqs". The locations of the .fastq files defined in the "McGarvey2022.csv" provided in this repo as "../fastqs/SRR4036047_1.fastq".

```bash
# load the module with fasterq-dump. note:  this will likely differ on your system
ml ncbi_sra
# create directory for the entire snakemake run
cd McGarvey2022_HiC_Analysis/
# create a subdirectory to which the .fastq files will be transferred
mkdir fastqs
cd fastqs
# transfer each pair of files for the following SRR entries. 
# In this example we use "sbatch --wrap" to run the fasterq-dump command on our hpc.
SRAIDS=( "SRR12008033" \
"SRR12008034" )
for t in ${SRAIDS[@]}; do
  sbatch --cpus-per-task 12 --wrap="fasterq-dump --split-files --threads 12 $t"
  done
# change back to the McGarvey2022_HiC_Analysis/ directory
cd ..
```

### 3.  Load slurm and miniconda
Note. The commands to do this will be different on your machine. These commands are specific to an HPC with slurm and miniconda modules installed.

```bash
module purge
ml slurm/20.02
ml miniconda/4.11.0
```
### 4.  Clone repository
```bash
cd McGarvey2022_HiC_Analysis
git clone https://github.com/SansamLab/HiC-SnakeMake.git
# rename folder with project name
mv HiC-SnakeMake/ McGarvey2022_HiC_Project_Folder/
# change directory into root of your project folder
cd McGarvey2022_HiC_Project_Folder
```
### 5.  Start the snakemake conda environment
#### 5.1.  Setup the snakemake conda environment
```bash
# sbatch submits a job script to slurm
# the --wrap option wraps the quoted command string in a sh script and submits
# -f is the location of the environment .yml file. 
## The relative path assumes that you are in the root directory of this repository.
# -p is the path where you want to install this environment
sbatch --mem 16G --wrap="conda env create -f workflow/envs/SnakemakeEnv.yml -p ../SnakemakeEnv" 
```

#### 5.2.  Activate the snakemake conda environment
```bash
# minimize conflicts, close all modules on slurm
module purge
ml slurm/20.02
ml miniconda/4.11.0
# activate the Snakemake Environment
conda activate ../SnakemakeEnv/
```

### 6. Modify the job-specific configuration files.

#### 6.1. Modify the config/config.yml file
The config.yml file is preconfigured for the test data set, so the file must be changed. A config file for the McGarvey data is included in this repository. To use that file rename it "config.yml".
```bash
rm config/config.yml
mv config/McGarvey2022_config.yml config/config.yml
```

#### 6.2. Modify the config/McGarvey_Samples.csv file
The McGarvey_Samples.csv file in the config folder has relative paths to the fastq files transferred from SRA. If the path differs in your system update it in McGarvey_Samples.csv.

### 7A. Run pipeline with conda environments (Alternative 1)
#### 7A.1. Install necessary conda environments
```
sbatch --mem 32G \
--wrap="\
snakemake \
--cores all \
--use-conda \
--conda-prefix condEnvs/ \
--conda-create-envs-only \
--conda-frontend conda"
```
#### 7A.2. Run pipeline with conda environments
```bash
While within the root directory of the repository clone, enter the following command.
sbatch --constraint=westmere \
--wrap="\
snakemake \
-R \
-j 999 \
--use-conda \
--conda-prefix condEnvs/ \
--conda-frontend conda \
--latency-wait 100 \
--cluster-config config/cluster_config.yml \
--cluster '\
sbatch \
-A {cluster.account} \
-p {cluster.partition} \
--cpus-per-task {cluster.cpus-per-task}  \
--mem {cluster.mem} \
--output {cluster.output}'"
```

### 7B. Run pipeline with installed modules (Alternative 2)
#### 7B.1. Modify Snakefile with modules installed on your hpc
Each rule in the workflow/Snakefile file has modules listed. These should be changed to match the names of the modules on your hpc. For example:
![rule change example](https://github.com/SansamLab/RepliTimer/blob/main/resources/ruleChangeExample.png)

#### 7B.2. Run pipeline with modules installed on hpc
While within the root directory of the repository clone, enter the following command.
```bash
sbatch --constraint=westmere \
--wrap="\
snakemake \
-R \
-j 999 \
--use-envmodules \
--latency-wait 100 \
--cluster-config config/cluster_config.yml \
--cluster '\
sbatch \
-A {cluster.account} \
-p {cluster.partition} \
--cpus-per-task {cluster.cpus-per-task}  \
--mem {cluster.mem} \
--output {cluster.output}'"
```
### 8.  Check output in the results/ directory
Use the tree program to list all of the files with sizes.
```bash
tree -hF results
```

You should get an output like this:
```
