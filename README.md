# PiSlice

Package in active development (BETA).

Estimate Pi and other population genomics statistics (PiN, PiS, Tajima's D, GC content) for a given list of windows coordinates.

Genomic analyses from chromosome to gene scale (e.g. explore the intron/exons structure).

The package is made for flexible use, and can be run as a Python package or a standalone CLI.

In the final version, you will be able to:
* Estimate directly statistics from appropriate input stream in STDIN (i.e. DNA sequence extracted on your own)
* Use a more contextual approach from a list of windows coordinates and input data files necessary to extract data (i.e; fasta, gff and/or vcf)
* Biologically aware computation capable to handle genome and gene architecture (e.g. introns/exons rank). PiSlice tries to recover automatically the structure of genes (i.e. splicing variants, exon rank, CDS phase and strand)

Output to STDOUT or files.


## Installation

### Dependencies

* ‘pandas‘
* ‘pysam‘ to read fasta files and ‘sammtools faidx‘ is a dependency of ‘pysam‘
* ‘cyvcf2‘ to read and manipulate population genetics dataset in vcf format
* ‘htslib‘ for bgzipped fasta


### Package

PiSlice has a setup and can be installed locally with pip. Run the following command in the PiSlice git directory.

‘‘‘
pip install .
‘‘‘

Or

‘‘‘
pip install -e .
‘‘‘

to make it editable after installation.


## Data & Input

The type of data depends on the statistic computed.

Data:
* Genomic data (i.e. DNA sequences) as a fasta file, bgzipped
* Polymorphism data (i.e. SNPs), as a vcf file, gzipped
* Annotation data, as a gff3 file, ncbi annotation format for parsing attributes

Query:
* a csv file with coordinates of genomic windows (chromosome, start, end)
* a parsed gff, genomic statistics will be computed for each feature (i.e. row)

A txt/csv file with genomic coordinates (chromosome, start, end in basepairs) is used as a template for genomic windows in which estimate the statistics.

Depending on the statistics to compute, fasta, gff and vcf files will be processed.
Genomic datasets (fasta) must be indexed (faidx/tabix) and compressed with bgzip.

### Fasta

### GFF

[![Canonical gene model according to the GFF3 specification](https://github.com/The-Sequence-Ontology/Specifications/blob/master/img/figure1.png)](https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md)

### vcf


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
* WindowScanr: sliding window analysis in R https://github.com/tavareshugo/WindowScanR

## Code convention

* snake_case
* scripts, functions and variables all begin by lower case
