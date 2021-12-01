# -*- coding: utf-8 -*-
"""
Input functions.

Functions in this module import .fasta, .gff, .vcf and coordinate files in appropriate formats
"""

from cyvcf2 import VCF
from pysam import FastaFile

def readVcf(vcf):
    """
    Read a vcf file and return a cyvcf2.VCF object
    Filename is given in argument
    """
    df = VCF(vcf)
    return df

def readFasta(fasta):
    """
    Read a fasta file and return a sequence object
    :param fasta:
    :return:
    """
    sequence = FastaFile(fasta)
    return sequence