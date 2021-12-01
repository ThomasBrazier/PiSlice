# PiSlice

Estimate Pi and other population genomics statistics (PiN, PiS, Tajima's D, GC content) for a given list of windows coordinates.

## Dependencies

* 'pysam' to read fasta files and 'sammtools faidx' is a dependency of 'pysam'
* 'cyvcf2' to read and manipulate population genetics dataset in vcf format

## Input

A txt/csv file with genomic coordinates (chromosome, start, end in basepairs) is used as a template for genomic windows in which estimate the statistics.

Depending on the statistics to compute, fasta, gff and vcf files will be processed.
Genomic datasets must be indexed (faidx/tabix) and compressed with bgzip.

## Statistics

* Pi
* PiN, PiS and PiN/PiS
* Tajima's D
* SNP density
* Gene density
* GC/GC1/GC2/GC3 content

## Output
