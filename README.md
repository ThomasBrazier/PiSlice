# PiSlice

Estimate Pi and other population genomics statistics (PiN, PiS, Tajima's D, GC content) for a given list of windows coordinates.

Genomic analyses from chromosome to intron/exon scale.

The package is made for flexible use, and can be run as a Python package or a standalone CLI.

You can:
* Estimate directly statistics from appropriate input stream in STDIN (i.e. DNA sequence extracted on your own)
* Use a more contextual approach from a list of windows coordinates and input data files necessary to extract data (i.e; fasta, gff and/or vcf)
* Biologically aware computation capable to handle genome and gene architecture (e.g. introns/exons rank). PiSlice tries to recover automatically the structure of genes (i.e. splicing variants, exon rank, CDS phase and strand)

Output to STDOUT or files.


## Dependencies

* 'pysam' to read fasta files and 'sammtools faidx' is a dependency of 'pysam'
* 'cyvcf2' to read and manipulate population genetics dataset in vcf format

## Data & Input

The type of data depends on the statistic computed.

Data:
* Genomic data (i.e. DNA sequences) as a fasta file, bgzipped
* Polymorphism data (i.e. SNPs), as a vcf file
* Annotation data, as a gff3 file, ncbi annotation format for parsing attributes

Query:
* a csv file with coordinates of genomic windows (chromosome, start, end)

A txt/csv file with genomic coordinates (chromosome, start, end in basepairs) is used as a template for genomic windows in which estimate the statistics.

Depending on the statistics to compute, fasta, gff and vcf files will be processed.
Genomic datasets must be indexed (faidx/tabix) and compressed with bgzip.

## Statistics

Genomes:
* Gene density
* GC/GC1/GC2/GC3 content
* GC/GC1/GC2/GC3 in first exon (does not exist yet!)
* Gene structure (e.g. length, number of exons)
* CpG
* Codon usage bias

Population genomics:
* Basic statistics: number of polymorphic sites, heterozygosity
* Pi
* PiN, PiS and PiN/PiS (PiN and PiS complicated to estimate!)
* Tajima's D
* SNP density
* SFS
* DFE


## Output

* A data frame (csv) which contains for each genomic window queried an estimate of population genomic statistics as required by user

* Extract sequences or coordinates of any feature of interest (e.g. CDS rank 1), ability to cut the genome in pieces (gene, splicing variants, exons, CDS). See features of gffutils.

### Data

### Figures

## Comparison with other packages

* BioPython
* Bio++
* sgkit
* tskit
* Genepop
* gffread to extract CDS or exons
* seqtk

## Code convention

* snake_case
* scripts, functions and variables all begin by lower case
