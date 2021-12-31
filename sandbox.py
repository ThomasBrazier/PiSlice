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
fasta_file = "data/Oryza_sativa_GCA_001433935.1/GCA_001433935.1_IRGSP-1.0_genomic.fna.gz"
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
# gff_file = "data/Oryza_sativa_GCA_001433935.1/GCA_001433935.1_IRGSP-1.0_genomic.gff.gz"
# gff_file = "data/Arabidopsis_thaliana_GCA_000001735.2/GCA_000001735.2_TAIR10.1_genomic.gff.gz"
# gff = input.read_gff(gff_file)
# gff = gff.gff.parse_attributes(infer_rank=True)

# gff.gff.feature("gene")
# gff.gff.region(1,20000,"CP002684.1")
# gff.gff.parent("exon-gnl|JCVI|mRNA.AT1G01010.1-2")
# gff.gff.children("gene-AT1G01020")
# gff.gff.rank(1)

save_filename = "data/Oryza_sativa_GCA_001433935.1/GCA_001433935.1_IRGSP-1.0_genomic.csv.gz"
#save_filename = "data/Arabidopsis_thaliana_GCA_000001735.2/GCA_000001735.2_TAIR10.1_genomic.csv.gz"
#input.write_gff2csv(gff, save_filename)
gff = input.read_gff(save_filename)

# Test statistics functions
import popstatistics as pop
# Compute a single GC and GC1, GC2, GC3 (i.e. single sequence or list of sequences)
chromosome = genome.references[0]
sequence = genome.sample_chromosome(chromosome)
len(sequence)
pop.gc(sequence)

chromosome = genome.references[0]
sequence = genome.sample_chromosome(chromosome)
start = 1
end = len(sequence)
pop.gc_codon(genome, gff, chromosome, start, end)


# Compute GC and GC1, GC2, GC3 for multiple sequences (multiple outputs)
# Create a data frame with all chromosomes full length (one chromosome per row)
nb_chromosome = 5
windows = pd.DataFrame({
    'Chromosome': genome.references[0:nb_chromosome],
    'Start': [1] * nb_chromosome,
    'End': list(map(lambda x: len(genome.fetch(x)), genome.references[0:nb_chromosome]))
})
results = pop.piSlice(windows=windows, statistics=["gene_count", "gc", "gc_codon"], fasta=genome, gff=gff)
results
# Results congruent with estimates of Ressayre et al. 2015 for A. thaliana and O. sativa

# Get GC for CDS rank 1 & next ones
results_rank1 = pop.piSlice(windows=windows, statistics=["gc", "gc_codon"], fasta=genome, gff=gff.gff.rank(1))
results_rank2 = pop.piSlice(windows=windows, statistics=["gc", "gc_codon"], fasta=genome, gff=gff.gff.rank(2))
results_rank3 = pop.piSlice(windows=windows, statistics=["gc", "gc_codon"], fasta=genome, gff=gff.gff.rank(3))



# Compute GC and GC1, GC2, GC3 for all CDS sequences (multiple outputs)
# Create a data frame with all CDS
# TODO Take care of chromosome names, bug to fix when no sequence is found for a chromosome name
gff_cds = gff[(gff['feature'] == "CDS") & (gff['seqname'].apply(lambda x: x in ["1","2","3","4","5","6","7","8","9","10"]))]
windows = pd.DataFrame({
    'Chromosome': list(gff_cds['seqname']),
    'Start': list(gff_cds['start']),
    'End': list(gff_cds['end'])
})
results = pop.piSlice(windows=windows, statistics=["gc", "gc_codon"], fasta=genome, gff=gff)
# Take times - TODO optimization




#--------------------------------------------------------------------------------------
# Try gffutils
#--------------------------------------------------------------------------------------
import gffutils
gff_file = "data/Oryza_sativa_GCA_001433935.1/GCA_001433935.1_IRGSP-1.0_genomic.gff.gz"
db_name = "data/Oryza_sativa_GCA_001433935.1/gff.db"
db = gffutils.create_db(gff_file, db_name, force=True, keep_order=True,
                        merge_strategy='merge', sort_attribute_values=True)
# Creating DB takes time and resources
# But do it once and then import existing DB
db = gffutils.FeatureDB(db_name, keep_order=True)
dir(db)

for i in db.featuretypes():
    print("Feature: %s: %d" % (i, db.count_features_of_type(i)))

for g in db.features_of_type('gene'):
    print(g)
gene = db['gene-OSNPB_120641600']
print(gene)

for i in db.all_features():
    print(i)

# import input
# gff_file = "data/Oryza_sativa_GCA_001433935.1/GCA_001433935.1_IRGSP-1.0_genomic.gff.gz"
# gff_db = input.gff(gff_file, merge_strategy='merge', sort_attribute_values=True)
# db_name = "data/Oryza_sativa_GCA_001433935.1/GCA_001433935.1_IRGSP-1.0_genomic.db"
# gff_db = input.gff(db_name)
# gffutils is powerful and useful for annotating features
# But:
# - do not provide exon/CDS rank
# - slow for creating database
# - do not infer when attributes are incomplete (e.g. from exon/CDS position)

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