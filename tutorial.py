import input as io
import statistics as stat

# Input/output
# Import a vcf
vcf = "data/Osativa.vcf.gz"
variants = io.readVcf(vcf)

# Explore the polymorphism dataset
# Samples
print(variants.samples)
# Chromosome names
print(variants.seqnames)
# Chromosome size
print(variants.seqlens)
# Get a 10Mb window on chromosome 1 for all samples
snps = variants('1:1000-10000')
# snps is a generator object
print(snps)
# Iterate over snps
for v in snps:
    print(str(v))
    # single sample of 0|1 in vcf becomes [[0, 1, True]]
    # 2 samples of 0/0 and 1|1 would be [[0, 0, False], [1, 1, True]]
    #print(v.genotypes)


# Import a fasta file
# https://onestopdataanalysis.com/read-fasta-file-python/
fasta = "data/Osativa_GCA001433935.fna.gz"
genome = io.readFasta(fasta)

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
type(chromosome)
len(chromosome)
# Or use the methods in pysam
chromosome = genome.fetch(genome.references[0])
type(chromosome)
len(chromosome)
# .fetch() can get a window within a chromosome
start = 2000
end = 2100
window = genome.fetch(genome.references[0], start, end)
type(window)
len(window)
# We got a 100kb window

# Statistics
# Compute Pi
# Subset SNPs in a given window
snps
# Estimate pi on polymorphism data
stat.pi(snps)

# Compute Pi in a series of non-overlapping windows