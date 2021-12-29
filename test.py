from cyvcf2 import VCF
from pysam import FastaFile
import numpy as np
import pandas
import gzip
from pandas import DataFrame
import re

import input
fasta_file = "data/Oryza_sativa_GCA_001433935.1/GCA_001433935.1_IRGSP-1.0_genomic.fna.gz"
#fasta_file = "data/Arabidopsis_thaliana_GCA_000001735.2/GCA_000001735.2_TAIR10.1_genomic.fna.gz"
fasta = input.fasta(fasta_file)
import popstatistics as pop
gff_file = "data/Oryza_sativa_GCA_001433935.1/GCA_001433935.1_IRGSP-1.0_genomic.gff.gz"
gff_file = "data/Arabidopsis_thaliana_GCA_000001735.2/GCA_000001735.2_TAIR10.1_genomic.gff.gz"
#gff = input.gff(gff_file)
gff = input.read_gff(gff_file)
#gff = gff[:2000]
#gff
gff = gff.gff.parse_attributes(infer_rank=True)
gff
gff.gff.feature("gene")
gff.gff.region(1,20000,"CP002684.1")
gff.gff.parent("exon-gnl|JCVI|mRNA.AT1G01010.1-2")
gff.gff.children("gene-AT1G01020")
gff.gff.rank(1)

save_filename = "data/Arabidopsis_thaliana_GCA_000001735.2/GCA_000001735.2_TAIR10.1_genomic.csv.gz"
input.write_gff2csv(gff, save_filename)
gff = input.read_gff(save_filename)