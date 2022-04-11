
### 1. Load necessary software
```bash
ml seqtk bwa
```

### 2. Copy GRCz11 fasta to local directory
```bash
sbatch --wrap=\
"wget -e robots=off --reject 'index.html' https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/035/GCF_000002035.6_GRCz11/GCF_000002035.6_GRCz11_genomic.fna.gz"
```

### 3. Get primary sequences (non-Alts) and rename chromosomes
sbatch --wrap=\
"zcat GCF_000002035.6_GRCz11_genomic.fna.gz | seqtk seq -l 3000000000 | grep Primary -A 1 | sed 's/>/>chrUn/g' | seqtk seq -l 100 | sed 's/^.*chromosome />chr/g' | sed 's/,.*$//g' | gzip > GCF_000002035.6_GRCz11_primary_genomic.fna.gz"

### 4. Index primary genome with bwa
sbatch --mem 32G --wrap="bwa index GCF_000002035.6_GRCz11_primary_genomic.fna.gz"
