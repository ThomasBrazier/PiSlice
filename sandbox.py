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
chromosome = genome.fetch(genome.references[0])
type(chromosome)
len(chromosome)

# Sample a sequence
start = 2000
end = 2100
chromosome = 1
window = genome.sample_sequence(genome.references[0], start, end)
window = genome.sample_sequence(chromosome, start, end)

type(window)
len(window)
# We got a 100kb window
#chromosome = genome.sample_chromosome(chromosome)

#--------------------------------------------------------------------------------------
# Import a GFF annotation file
#--------------------------------------------------------------------------------------
import input
gff_file = "data/Osativa_GCA001433935.gff3.gz"
gff = input.gff(gff_file)

# from BCBio import GFF
# import gzip
# gff_file = "data/Osativa_GCA001433935.gff3"
# file = open(gff_file, 'r', encoding="utf-8")
# # Limiting to features of interest
# limit_info = dict(gff_id=["1"], gff_source=["CDS"])
# for rec in GFF.parse(file, limit_info=limit_info):
#     print(rec.seq)
# file.close()
#
#
# # Read gff as a Pandas DataFrame
# import pandas as pd
# import gzip
# gff_file = "data/Osativa_GCA001433935.gff3.gz"
# file = open(gff_file, 'r', encoding="utf-8")
# gff = pd.read_csv(gff_file, sep="\t", comment="#", low_memory=False,
#                   names=["seqname", "source", "feature", "start", "end",
#                          "score", "strand", "frame", "attribute"])
# file.close()
# gff

sequence = genome.sample_chromosome(chromosome)
features = gff
start = 2000
end = 12000
chromosome = 1
pop.gcpos(sequence, features, chromosome, start, end)



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