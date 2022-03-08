The "danRer11.chrom.sizes.modified" file has only the numbered chromosomes from danRer11.  All "alt" and unknown contigs were removed with the following:

```
grep -P '(chr[[:digit:]])' /Volumes/hts_core/Shared/zebrafish/danRer11/danRer11.chrom.sizes \
| grep -v alt > resources/chromosomeSizes/danRer11.chrom.sizes.modified
```