# Gene Fusion Program
## Usage
`rna_fastq` folder in the same directory as `gfp.sm` will need to be configured with paired-end RNA-seq fastq fies `[A-Za-z0-9]+_R1.fastq.gz` and `[A-Za-z0-9]+_R2.fastq.gz`. 

Some reference files and Annovar may need some configuration as well as curation of panel of normal for filtering. 

```
$ snakemake -np -s gfp.sm -j 24
```
