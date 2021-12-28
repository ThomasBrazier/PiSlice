from cyvcf2 import VCF
from pysam import FastaFile
import numpy as np
import pandas
import gzip
from pandas import DataFrame
import re

import pandas as pd
import input
fasta_file = "data/Oryza_sativa_GCA_001433935.1/GCA_001433935.1_IRGSP-1.0_genomic.fna.gz"
#fasta_file = "data/Arabidopsis_thaliana_GCA_000001735.2/GCA_000001735.2_TAIR10.1_genomic.fna.gz"
fasta = input.fasta(fasta_file)
import popstatistics as pop
gff_file = "data/Oryza_sativa_GCA_001433935.1/GCA_001433935.1_IRGSP-1.0_genomic.gff.gz"
gff_file = "data/Arabidopsis_thaliana_GCA_000001735.2/GCA_000001735.2_TAIR10.1_genomic.gff.gz"
gff = input.gff(gff_file)
#gff = gff.set_index(["feature"])
#gff = gff[:2000]
attr = gff['attribute']
