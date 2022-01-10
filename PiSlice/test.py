import numpy as np
import pandas as pd

# Test the GC function
from PiSlice.popstatistics import gc
test = "ATCGTAGCTGGGCTAGCTGATGCGCGCG" # 18 GC on 28 nucleotides
gc(test) == 18/28
test = "ATCGTNGCNGGGCTAGCTGATGCGCGCG" # 18 GC on 26 nucleotides + 2 Ns not to consider
gc(test) == 18/26
gc("ATC") is np.nan
gc("") is np.nan

# Test the CpG function
from PiSlice.popstatistics import cpg
test = "ATCGTAGCTGGGCTAGCTGATGCGCGCG" # 4 CpG sites on 14 sites
cpg(test) == 4/14
test = "ATCGTAGCNNGGCTAGCTGATGCGCGCG" # 4 CpG sites on 13 sites + 1 site NN
cpg(test) == 4/13
import PiSlice.input as input
# Test short sequence (< 6bp)
cpg("ATCG") is np.nan
# Test on real data
fasta_file = "data/Oryza_sativa_GCA_001433935.1/GCA_001433935.1_IRGSP-1.0_genomic.fna.gz"
genome = input.fasta(fasta_file)
nb_chromosome = 5
windows = pd.DataFrame({
    'Chromosome': genome.references[0:nb_chromosome],
    'Start': [1] * nb_chromosome,
    'End': list(map(lambda x: len(genome.fetch(x)), genome.references[0:nb_chromosome]))
})
from PiSlice.popstatistics import piSlice
piSlice(windows, statistics=["cpg"], fasta=genome)

windows = pd.DataFrame({
    'Chromosome': genome.references[0:2:1],
    'Start': [1000, 1000],
    'End': [2000, 1005]
})
piSlice(windows, statistics=["cpg"], fasta=genome)
