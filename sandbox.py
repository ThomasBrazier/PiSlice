import pandas as pd

import input
import popstatistics as pop

# Input/output
#--------------------------------------------------------------------------------------
# Import a vcf
#--------------------------------------------------------------------------------------
vcf_file = "data/Osativa.vcf.gz"
variants = input.vcf(vcf_file, strict_gt=True)
# Explore the polymorphism dataset
# Samples
print(variants.samples)
# Chromosome names
print(variants.seqnames)
# Chromosome size
print(variants.seqlens)
# Sample individuals
#variants.sample_individual('IRGC121316@c88dcbba.0')

# Sample region
start = 1000
end = 10000
chromosome = 1
# Get a 10Mb window on chromosome 1 for all samples
snps = variants('1:1000-10000')
# Or use the method .sample_variant()
snps = variants.sample_variant(1, 1000, 10000)
# snps is a generator object
# Iterate over snps
#for v in variants:
    #print(str(v))
    # single sample of 0|1 in vcf becomes [[0, 1, True]]
    # 2 samples of 0/0 and 1|1 would be [[0, 0, False], [1, 1, True]]
    #print(v.genotypes)

# Construct a genotype dataset for further analyses
geno = input.genotype(variants, chromosome, start, end)



#--------------------------------------------------------------------------------------
# Import a fasta file
#--------------------------------------------------------------------------------------
# https://onestopdataanalysis.com/read-fasta-file-python/
import input
fasta_file = "data/Osativa_GCA001433935.fna.gz"
fasta_file = "data/Zmays_GCA_000005005.6.fa.gz"
genome = input.fasta(fasta_file)

# Explore the genome dataset
# How many chromosomes?
print(genome.nreferences)
# Get chromosome names
print(genome.references)
# Get chromosome sizes
print(genome.lengths)

# Get a sequence in a given window (start-end)
# Search on chromosome 1
# The sequence retrieved is a character string
chromosome = genome[genome.references[0]]
chromosome = genome["AP014957.1"]
type(chromosome)
len(chromosome)
# Or use the methods in pysam
chromosome = genome.fetch(genome.references[1])
type(chromosome)
len(chromosome)

# Sample a sequence
start = 1
end = 10000
chromosome = genome.references[1]
window = genome.sample_sequence(chromosome, start, end)

type(window)
len(window)
# We got a 100kb window
#chromosome = genome.sample_chromosome(chromosome)

#--------------------------------------------------------------------------------------
# Import a GFF annotation file
#--------------------------------------------------------------------------------------
import pandas as pd
import input
#fasta_file = "data/Osativa_GCA001433935.fna.gz"
fasta_file = "data/Zmays_GCA_000005005.6.fa.gz"
#fasta_file = "data/Athaliana_GCA000001735.fna.gz"
fasta = input.fasta(fasta_file)
import popstatistics as pop
#gff_file = "data/Osativa_GCA001433935.gff3.gz"
gff_file = "data/Zmays_GCA_000005005.6.gff3.gz"
#gff_file = "data/Athaliana_GCA000001735.gff3.gz"
gff = input.gff(gff_file)
### !! gff and fasta must have the same chromosome names


# Compute a single GC and GC1, GC2, GC3 (i.e. single sequence or list of sequences)
chromosome = fasta.references[0]
sequence = fasta.sample_chromosome(chromosome)
len(sequence)
pop.gc(sequence)

chromosome = fasta.references[0]
sequence = fasta.sample_chromosome(chromosome)
chromosome = "1"
start = 1
end = len(sequence)
pop.gc_cds(fasta, gff, chromosome, start, end)


# Compute GC and GC1, GC2, GC3 for multiple sequences (multiple outputs)
# Create a data frame with all chromosomes full length (one chromosome per row)
nb_chromosome = 5
windows = pd.DataFrame({
    'Chromosome': fasta.references[0:nb_chromosome],
    'Start': [1] * nb_chromosome,
    'End': list(map(lambda x: len(fasta.fetch(x)), fasta.references[0:nb_chromosome]))
})
results = pop.piSlice(fasta=fasta, gff=gff, windows=windows, statistics=["gc_cds"])

# Compute GC and GC1, GC2, GC3 for all CDS sequences (multiple outputs)
# Create a data frame with all CDS
# TODO Take care of chromosome names, bug to fix when no sequence is found for a chromosome name
gff_cds = gff[(gff['feature'] == "CDS") & (gff['seqname'].apply(lambda x: x in ["1","2","3","4","5","6","7","8","9","10"]))]
windows = pd.DataFrame({
    'Chromosome': list(gff_cds['seqname']),
    'Start': list(gff_cds['start']),
    'End': list(gff_cds['end'])
})
results = pop.piSlice(fasta=fasta, gff=gff, windows=windows, statistics=["gc"])
# Take times - TODO optimization





#--------------------------------------------------------------------------------------
# Statistics
#--------------------------------------------------------------------------------------
# GC content
pop.gc(window)
# GC exact

# GC1, GC2, GC3






# Compute Pi
# Subset SNPs in a given window
snps
# Estimate pi on polymorphism data
pop.pi(snps)

# Compute Pi in a series of non-overlapping windows