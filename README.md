# Identification of stress-granule RNA by TRIBE 
This repository contains a Nextflow pipeline for processing raw (sc) TRIBE data

### General outline
Raw reads are processed closely following the [GATK best-practices workflow for RNAseq short variant discovery](https://gatk.broadinstitute.org/hc/en-us/articles/360035531192-RNAseq-short-variant-discovery-SNPs-Indels-). Raw reads go through the following steps:
1. Adapter and homopolymer sequence trimming
2. Depleting rRNA sequences
3. Aligning to the reference genome
4. Deduplication
5. Variant identification and filtering

### Basic usage
Create a nextflow.config in your data directory, setting, e.g.,:
```
params {
  // Workflow flags
  reads = "fastq/*_R{1,2}.fastq.gz"
  outdir = 'TRIBE'
  genomeAnnotations = "group_references/ensembl/95/drosophila_melanogaster/Drosophila_melanogaster.BDGP6.95.gtf"
  genomeFasta = "group_references/ensembl/95/drosophila_melanogaster/Drosophila_melanogaster.BDGP6.dna.toplevel.ERCC92.fa"
  depleteFasta = "dmel_rRNA.fa"
  knownVariants = "group_references/ensembl/95/drosophila_melanogaster/drosophila_melanogaster.vcf.gz"
  intervalList = "group_references/ensembl/95/drosophila_melanogaster/Drosophila_melanogaster.BDGP6.dna.toplevel.ERCC92.complete.interval_list"
  exonsIntervals = "reference/95/Drosophila_melanogaster.BDGP6.95_exons.interval_list"
  intronsIntervals = "reference/95/Drosophila_melanogaster.BDGP6.95_introns.interval_list"
  exonsBed = "group_references/ensembl/95/drosophila_melanogaster/exons_unique.bed.gz"
}
```

and running `nextflow run main.nf -resume`

You may process bulk data by adding `-profile bulk`

