# -*- coding: utf-8 -*-
"""
Input

Define classes for genomic objects like fasta, gff and vcf

Functions in this module import .fasta, .gff, .vcf and coordinate files in appropriate formats
"""

# TODO Create a global unified and consistent 'genomic' object that contain all data in slots: genome sequence, features, snps and associated metadata

from cyvcf2 import VCF
from pysam import FastaFile
import numpy as np
import pandas
import gzip
from pandas import DataFrame
import re

class vcf(VCF):
    """
    Read a vcf file and return a 'vcf' object
    inheriting the 'VCF' class from 'cyvcf2'
    Filename is given in argument when creating the 'vcf' object
    """

    def sample_individual(self, samples):
        """
        Sample individuals and modify the 'vcf' object
        """
        # TODO copy the vcf instead on changing it directly
        # Make a list individuals to sample
        # TODO implement methods to sample individuals by indexes or names
        # Subset individuals
        self.set_samples(samples)

    def sample_variant(self, chromosome, start, end):
        """
        Sample variants on a chromosome within start-end coordinates (bp)
         and return a generator object to iterate over the 'vcf' object
        """
        # Make a list of loci to sample
        # TODO
        # Subset loci
        snps = self(":".join([str(chromosome), "-".join([str(start), str(end)])]))
        return snps



class genotype:
    """
    Create 'genotype' object from a 'vcf' object which contains:
    * gt: genotypes as a 2D numpy array
    * samples: sample names
    * loci: locus names
    """
    def __init__(self, variants, chromosome, start, end):
        # Subset a region
        snps = variants.sample_variant(chromosome, start, end)
        # Get names of loci
        loci = []
        for variant in snps:
            loci.append(variant.ID)
        # Set dimensions
        df = []
        # For loop, because cannot unpack non-iterable cyvcf2.cyvcf2.Variant object
        snps = variants.sample_variant(chromosome, start, end)
        for variant in snps:
            df.append(list(variant.gt_types))
        df = np.array(df)
        self.gt = df
        self.samples = variants.samples
        self.loci = loci



class fasta(FastaFile):
    """
    Read a fasta file and return a 'fasta' sequence object
    inheriting the 'FastaFile' class from 'pysam'
    Filename is given in argument when creating the 'fasta' object
    """
    def sample_chromosome(self, chromosome):
        """
        Sample a whole chromosome sequence
        :param chromosome: str, the name of the chromosome
        :return: A DNA sequence of the whole chromosome (str)
        """
        # Chromosome can be either an integer (index) or a string (name)
        return self.fetch(chromosome)


    def sample_sequence(self, chromosome, start, end):
        """
        Sample a DNA sequence on a chromosome within start-end coordinates (bp)
        on a 1-offset (sequence begins at 1)
        :param chromosome: str, the name of the chromosome
        :param start: int, start position in the sequence
        :param end: int, end position in the sequence
        :return: A DNA sequence within boundaries start-end on a given chromosome (str)
        """
        # Chromosome can be either an integer (index) or a string (name)
        # start - 1 because we work in 1-bp offset while python is 0-bp offset
        return self.fetch(chromosome, start - 1, end)


    def sample_feature(self, gff, feature, chromosome, start, end):
        """
        Sample a list of DNA sequences corresponding to features (e.g. CDS)
        within a chromosome or start-end boundaries
        :param gff: gff, a Data Frame with genome annotation
        :param feature: str, the type of feature
        :param chromosome: str, the name of the chromosome
        :param start: int, start position in the sequence
        :param end: int, end position in the sequence
        :return: str, list of DNA sequences
        """



class gff(DataFrame):
    """
    Read a gff file and return a pandas 'DataFrame' object
    Filename is given in argument
    """
    def __init__(self, gff_file):
        file = gzip.open(gff_file, 'r')
        super(gff, self).__init__(pandas.read_csv(file, sep="\t", comment="#", low_memory=False,
                           names=["seqname", "source", "feature", "start", "end",
                                  "score", "strand", "frame", "attribute"]))
        file.close()






def gff_parse_attributes(gff):
    """
    Parse the column attributes of a gff
    :param gff: gff, a gff file based on Pandas DataFrame
    """
    attr = gff['attribute']
    # Parse first the available information in attributes
    gene_id = []
    mRNA_id = []
    for a in attr:
        try:
            parent = re.findall(";Parent=[A-Za-z0-9\\.\\-\\|]*;", a)[0]
        except:
            parent = ""
        parent = parent.replace(";", "")
        parent = parent.replace("Parent=", "")
        try:
            gene_id.append(parent.split("|")[0])
        except:
            gene_id.append(None)
        try:
            mRNA_id.append([id for id in parent.split("|") if "mRNA" in id][0])
        except:
            mRNA_id.append(None)
    # Infer exon parent from position or mRNA if information

    # Infer CDS parents from position

    # Infer rank from position or parent

    # Complete the data frame
    gff["gene_id"] = gene_id
    gff["mRNA_id"] = mRNA_id

    return(gff)

