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
        """ Sample a whole chromosome sequence """
        # Chromosome can be either an integer (index) or a string (name)
        if isinstance(chromosome, str):
            return self.fetch(chromosome)
        elif isinstance(chromosome, int):
            return self.fetch(self.references[chromosome-1])
        else:
            print("Chromosome must be integer (index) or character (name).")


    def sample_sequence(self, chromosome, start, end):
        """ Sample a DNA sequence on a chromosome within start-end coordinates (bp) """
        # Chromosome can be either an integer (index) or a string (name)
        if isinstance(chromosome, str):
            return self.fetch(chromosome, start, end)
        elif isinstance(chromosome, int):
            return self.fetch(self.references[chromosome-1], start, end)
        else:
            print("Chromosome must be integer (index) or character (name).")



def gff(gff_file):
    """
    Read a gff file and return a pandas 'DataFrame' object
    Filename is given in argument
    """
    # def __init__(self):
    #     file = open(gff_file, 'r', encoding="utf-8")
    #     self.df = gff.read_csv(gff_file, sep="\t", comment="#", low_memory=False,
    #                       names=["seqname", "source", "feature", "start", "end",
    #                              "score", "strand", "frame", "attribute"])
    #     file.close()
    file = gzip.open(gff_file, 'r')
    gff = pandas.read_csv(file, sep="\t", comment="#", low_memory=False,
                          names=["seqname", "source", "feature", "start", "end",
                                 "score", "strand", "frame", "attribute"])
    file.close()

    return gff

