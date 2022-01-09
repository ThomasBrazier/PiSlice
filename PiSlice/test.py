

# Test the CpG function
from popstatistics import cpg
test = "ATCGTAGCTGGGCTAGCTGATGCGCGCG" # 4 CpG sites on 14 sites = 0.29 CpG density
cpg(test) == 4/14
import input
import numpy as np
# Test short sequence (< 6bp)
cpg("ATCG") is np.nan
# Test on real data
fasta_file = "data/Oryza_sativa_GCA_001433935.1/GCA_001433935.1_IRGSP-1.0_genomic.fna.gz"
genome = input.fasta(fasta_file)
nb_chromosome = 5
import pandas as pd
windows = pd.DataFrame({
    'Chromosome': genome.references[0:nb_chromosome],
    'Start': [1] * nb_chromosome,
    'End': list(map(lambda x: len(genome.fetch(x)), genome.references[0:nb_chromosome]))
})
from popstatistics import piSlice
piSlice(windows, statistics=["cpg"], fasta=genome)

windows = pd.DataFrame({
    'Chromosome': genome.references[0:2:1],
    'Start': [1000, 1000],
    'End': [2000, 1005]
})
piSlice(windows, statistics=["cpg"], fasta=genome)
