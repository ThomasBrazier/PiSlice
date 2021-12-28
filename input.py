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
    :param parse: False, if True, then attributes are parsed and the relationships parent-children are found between gene-mRNA-exon-CDS
    """
    def __init__(self, gff_file, parse=True):
        file = gzip.open(gff_file, 'r')
        super(gff, self).__init__(pandas.read_csv(file, sep="\t", comment="#", low_memory=False,
                           names=["seqname", "source", "feature", "start", "end",
                                  "score", "strand", "frame", "attribute"]))
        file.close()
        # TODO Parse attributes
        if (parse):
            self = self.parse_attributes()
        else:
            self["id"] = None
            self["parent"] = None
            self["name"] = None

    def gff2pandas(self):
        """
        Display the gff as an original Pandas DataFrame
        """
        frame = {"seqname": self["seqname"], "source": self["source"], "feature": self["feature"],
                "start": self["start"], "end": self["end"], "score": self["score"], "strand": self["strand"],
                "frame": self["frame"], "attribute": self["attribute"], "id": self["id"],
                 "parent": self["parent"], "name": self["name"]}
        df = pandas.DataFrame(frame)
        return(df)

    def parse_attributes(self):
        """
        Parse the column attributes of a gff
        :param gff: gff, a gff file based on Pandas DataFrame
        :return: gff, a gff augmented with three columns for attibutes
        """
        gff = self
        attr = self['attribute']
        # Parse first the available information in the attribute field
        id = []
        parent = []
        name = []
        for i, a in attr.iteritems():
            try:
                id_name = re.findall("ID=[A-Za-z0-9\\.\\-\\|\\_]*;", a)[0]
            except:
                id_name = ""
            id_name = id_name.replace(";", "")
            id_name = id_name.replace("ID=", "")
            try:
                id.append(id_name)
            except:
                id.append(None)

            try:
                parent_term = re.findall(";Parent=[A-Za-z0-9\\.\\-\\|\\_]*;", a)[0]
            except:
                parent_term = ""
            parent_term = parent_term.replace(";", "")
            parent_term = parent_term.replace("Parent=", "")
            try:
                parent.append(parent_term)
            except:
                parent.append(None)

            try:
                name_term = re.findall(";Name=[A-Za-z0-9\\.\\-\\|\\_]*[;]*", a)[0]
            except:
                name_term = ""
            name_term = name_term.replace(";", "")
            name_term = name_term.replace("Name=", "")
            try:
                name.append(name_term)
            except:
                name.append(None)
        gff["id"] = id
        gff["parent"] = parent
        gff["name"] = name

        # Infer exon/CDS rank from position or parent
        # CDS inherits rank of its parent exon

        return(gff)


    def region(self, start, end, chromosome):
        """
        Return a new gff with only features (rows) overlapping the queried genomic region
        """
        return(subset)


    def parent(self, id):
        """
        Return a new gff with the parent of given ids
        """
        return(subset)


    def children(self, id):
        """
        Return a new gff with children of given ids
        """
        return(subset)


    def feature(self, type):
        """
        Return a new gff with only features of a given type
        """
        return(subset)


    def summary(selfself):
        """
        Print summary statistics on the gff (e.g. count of features, number of chromosomes and chromosomes length)
        """









# DEPRECATED BELOW THIS LINE


# class gff(gffutils.FeatureDB):
#     """
#     Read a gff file and return a gffutils database
#     GFF filename or db is given in argument
#     If GFF, create a new database
#     otherwise, if db provided, import existing database
#     """
#     def __init__(self, file_name, *args, **kwargs):
#         if ".gff" in file_name:
#             super(gff, self).__init__(gffutils.create_db(file_name, re.sub(".gff.gz", ".db", file_name),
#                                                          force=True, keep_order=True, *args, **kwargs))
#         elif ".db" in file_name:
#             print(file_name)
#             super(gff, self).__init__(FeatureDB(file_name, keep_order=True, *args, **kwargs))

    # def parse_attributes(self):
    #     """
    #     Parse the column attributes of a gff
    #     :param gff: gff, a gff file based on Pandas DataFrame
    #     :return: dict, attributes parsed in a dictionary
    #     """
    #     gff = self
    #     attr = self['attribute']
    #     # Parse first the available information in the attribute field
    #     id = []
    #     gene_id = []
    #     for i, a in attr.iteritems():
    #         try:
    #             id_name = re.findall("ID=[A-Za-z0-9\\.\\-\\|\\_]*;", a)[0]
    #         except:
    #             id_name = ""
    #         id_name = id_name.replace(";", "")
    #         id_name = id_name.replace("ID=", "")
    #         try:
    #             id.append(id_name)
    #         except:
    #             id.append(None)
    #
    #         try:
    #             parent = re.findall(";Parent=[A-Za-z0-9\\.\\-\\|\\_]*;", a)[0]
    #         except:
    #             parent = ""
    #         parent = parent.replace(";", "")
    #         parent = parent.replace("Parent=", "")
    #         try:
    #             gene_id.append(parent.split("|")[0])
    #         except:
    #             gene_id.append(None)
    #     gff["id"] = id
    #     gff["gene_id"] = gene_id
    #
    #     mRNA_id = []
    #     for i, a in attr.iteritems():
    #         try:
    #             if (("|" not in parent) & (gff["feature"][i] in ("mRNA", "CDS", "exon"))):
    #                 if (gff["gene_id"][i] not in [None, '']):
    #                     mRNA_id.append(gff["gene_id"][i])
    #                 else:
    #                     gff_gene = gff[gff["feature"] == "gene"]
    #                     mRNA_parent = gff_gene[((gff_gene["start"] <= gff["start"][i]) & (gff_gene["end"] >= gff["end"][i])
    #                                & (gff_gene["seqname"].isin([gff["seqname"][i]])))]
    #                     # Conflict with overlapping genes
    #                     # Take the gene in the row just before (assume ordered annotation)
    #                     mRNA_id.append(list(mRNA_parent["id"])[0])
    #                     #mRNA_id.append([x for x, y in zip(id, mRNA_parent) if y][0])
    #             elif ("|" in parent):
    #                 mRNA_parent = [id for id in parent.split("|") if "mRNA" in id][0]
    #                 mRNA_id.append(mRNA_parent)
    #             else:
    #                 mRNA_id.append(None)
    #         except:
    #             mRNA_id.append(None)
    #         gff["mRNA_id"] = mRNA_id
    #     # Correct mRNA ids that were not correctly parsed
    #     # TODO Optim: computation can be long
    #     for i, a in attr.iteritems():
    #         try:
    #             if (mRNA_id[i] == None):
    #                 # If not possible with attributes, infer mRNA parent from mRNA-gene position intersection
    #                 # mRNA parent is ID of the gene containing mRNA
    #                 mRNA_parent = ((gff["start"] <= gff["start"][i]) & (gff["end"] >= gff["end"][i])
    #                                             & (gff["seqname"] == gff["seqname"][i])
    #                                             & (gff["feature"] == "gene"))
    #                 mRNA_id[i] = [x for x, y in zip(id, mRNA_parent) if y][0]
    #         except:
    #             pass
    #     # Infer exon parent from position or mRNA if information available
    #
    #     # Infer CDS parents from position
    #
    #     # Infer rank from position or parent
    #     # CDS inherits rank of its parent exon
    #
    #     # Complete the data frame
    #     attributes = {}
    #     attributes["id"] = id
    #     attributes["gene_id"] = gene_id
    #     attributes["mRNA_id"] = mRNA_id
    #     attributes["transcript_rank"] = transcript_rank
    #     attributes["exon_id"] = exon_id
    #     attributes["cds_id"] = cds_id
    #     attributes["rank"] = rank
    #     return (attributes)