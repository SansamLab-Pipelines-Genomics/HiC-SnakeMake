## Instructions for making the restriction enzyme file

The script below was taken from Juicer.
```bash
workflow/scripts/python generate_site_positions.py \
	'HindIII' \`							  `# restriction enzyme. Note that you may need alter the script file to add your restriction enzyme.`
	'hg38' \`								  `# your name for the genome`
	'/pathToGenome/genomeSequenceFastaFile.fa'`# path and name of genome fasta`
```
