# PiSlice

Estimate Pi and other population genomics statistics (PiN, PiS, Tajima's D, GC content) for a given list of windows coordinates.

## Dependencies

* 'pysam' to read fasta files and 'sammtools faidx' is a dependency of 'pysam'
* 'cyvcf2' to read and manipulate population genetics dataset in vcf format

## Data & Input

The type of data depends on the statistic computed.

Data:
* Genomic data (i.e. DNA sequences) as a fasta file
* Polymorphism data (i.e. SNPs), as a vcf file
* Annotation data, as a gff file

Query:
* a csv file with coordinates of genomic windows (chromosome, start, end)

A txt/csv file with genomic coordinates (chromosome, start, end in basepairs) is used as a template for genomic windows in which estimate the statistics.

Depending on the statistics to compute, fasta, gff and vcf files will be processed.
Genomic datasets must be indexed (faidx/tabix) and compressed with bgzip.

## Statistics

Genomes:
* Gene density
* GC/GC1/GC2/GC3 content
* GC/GC1/GC2/GC3 in first exon
* Gene structure (e.g. length, number of exons)

Population genomics:
* Basic statistics: number of polymorphic sites, heterozygosity
* Pi
* PiN, PiS and PiN/PiS
* Tajima's D
* SNP density

## Output

* A data frame (csv) which contains for each genomic window queried an estimate of population genomic statistics as required by user

### Data

### Figures

## Comparison with other packages

* BioPython
* Bio++
* sgkit
* tskit
* Genepop

## Code convention

* snake_case
* scripts, functions and variables all begin by lower case
