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
import mapply
import multiprocessing


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


# Extending Pandas with new methods for gff
# Accessible by the gff namespace
@pandas.api.extensions.register_dataframe_accessor("gff")
class GffAccessor:
    """
    Read a gff file and return a pandas 'DataFrame' object
    Filename is given in argument
    :param parse: False, if True, then attributes are parsed and the relationships parent-children are found between gene-mRNA-exon-CDS
    """
    def __init__(self, pandas_obj):
        #self._validate(pandas_obj)
        self._obj = pandas_obj

    # @staticmethod
    # def _validate(obj):
    #     # verify there is a column latitude and a column longitude
    #     if "latitude" not in obj.columns or "longitude" not in obj.columns:
    #         raise AttributeError("Must have 'latitude' and 'longitude'.")

    def parse_attributes(self, infer_rank=False, n_cpus=0):
        """
        Parse the column attributes of a gff
        :param gff: gff, a gff file based on Pandas DataFrame
        :return: gff, a gff augmented with three columns for attibutes
        """
        gff_obj = self._obj.copy(deep=True)
        if (n_cpus == 0):
            n_cpus = multiprocessing.cpu_count()
        # Parse first the available information in the attribute field
        id = [""] * gff_obj.shape[0]
        parent = [""] * gff_obj.shape[0]
        name = [""] * gff_obj.shape[0]
        for i, r in gff_obj.iterrows():
            a = r["attribute"]
            try:
                id_name = re.findall("ID=[A-Za-z0-9\\.\\-\\|\\_]*;", a)[0]
            except:
                id_name = ""
            id_name = id_name.replace(";", "")
            id_name = id_name.replace("ID=", "")
            try:
                id[i] = id_name
            except:
                id[i] = None

            try:
                parent_term = re.findall(";Parent=[A-Za-z0-9\\.\\-\\|\\_]*;", a)[0]
            except:
                parent_term = ""
            parent_term = parent_term.replace(";", "")
            parent_term = parent_term.replace("Parent=", "")
            try:
                parent[i] = parent_term
            except:
                parent[i] = None

            try:
                name_term = re.findall(";Name=[A-Za-z0-9\\.\\-\\|\\_]*[;]*", a)[0]
            except:
                name_term = ""
            name_term = name_term.replace(";", "")
            name_term = name_term.replace("Name=", "")
            try:
                name[i] = name_term
            except:
                name[i] = None
        gff_obj["id"] = id
        gff_obj["parent"] = parent
        gff_obj["name"] = name

        # Infer exon/CDS rank from position or parent
        # CDS inherits rank of its parent exon
        # rank = [0] * gff_obj.shape[0]
        def rank_inference(gff_obj, x):
            # Optim: do not create new objects
            filter_ = ((gff_obj["feature"] == x["feature"]) & (gff_obj["parent"] == x["parent"]))
            #set = gff_obj.gff.children(x["parent"])
            #set = set.gff.feature(x["feature"])
            #gff_obj[(gff_obj["feature"] == x["feature"]) & (gff_obj["parent"] == x["parent"])]["start"]
            if (x["strand"] == "+"):
            # Rank: how many sequences of the same feature and parent are before that one?
                rk = sum(bool(z) for z in [s <= int(x["start"]) for s in list(gff_obj[filter_]["start"])])
            elif (x["strand"] == "-"):
            # Rank: how many sequences of the same feature and parent are after that one?
            # Read in the opposite direction
                rk = sum(bool(z) for z in [s >= int(x["start"]) for s in list(gff_obj[filter_]["start"])])
            else:
                rk = 0
            return(rk)

        if (infer_rank):
            # # TODO optim at this step: long time for iterations
            # for i, a in gff_obj.iterrows():
            #     if (a["feature"] in ["exon", "CDS"]):
            #         set = gff_obj.gff.children(a["parent"])
            #         set = set.gff.feature(a["feature"])
            #         if (a["strand"] == "+"):
            #             # Rank: how many sequences of the same feature and parent are before that one?
            #             rk = sum(bool(x) for x in [s <= int(a["start"]) for s in list(set["start"])])
            #         elif (a["strand"] == "-"):
            #             # Rank: how many sequences of the same feature and parent are after that one?
            #             # Read in the opposite direction
            #             rk = sum(bool(x) for x in [s >= int(a["start"]) for s in list(set["start"])])
            #     else:
            #         rk = 0
            #     rank[i] = int(rk)
            # TODO vectorization
            mapply.init(n_workers=n_cpus)
            rank = gff_obj.mapply(lambda x: rank_inference(gff_obj, x) if x["feature"] in ["exon", "CDS"] else 0,
                                 axis=1)
        else:
            rank = [0] * gff_obj.shape[0]

        gff_obj["rank"] = rank
        return(gff_obj)


    def region(self, start, end, seq):
        """
        Return a new gff object with only features (rows) completely within the queried genomic region
        """
        subset = self._obj[(self._obj["start"] >= start) & (self._obj["end"] >= end) & (self._obj["seqname"] == seq)]
        return(subset)


    def parent(self, id):
        """
        Return a new gff object with the parents of given children ids
        """
        if (type(id) == str):
            id = [id]
        parent = list(set(self._obj[self._obj["id"].isin(id)]["parent"]))
        subset = self._obj[self._obj["id"].isin(parent)]
        return(subset)


    def children(self, id):
        """
        Return a new gff object with children of given parent ids
        """
        if (type(id) == str):
            id = [id]
        subset = self._obj[self._obj["parent"].isin(id)]
        return(subset)


    def feature(self, feature):
        """
        Return a new gff object with only features of a given feature or list of features
        """
        if (type(feature) == str):
            feature = [feature]
        subset = self._obj[self._obj["feature"].isin(feature)]
        return(subset)

    def rank(self, rank):
        """
        Return a new gff object with only exon/CDS of a given rank
        """
        subset = self._obj.gff.feature(["exon", "CDS"])
        if (type(rank) == int):
            rank = [rank]
        subset = subset[subset["rank"].isin(rank)]
        return(subset)


    def summary(self):
        """
        Print summary statistics on the gff (e.g. count of features, number of chromosomes and chromosomes length)
        """


def read_gff(gff_file, parse=False):
    if (".gff" in gff_file):
        file = gzip.open(gff_file, 'r')
        gff = pandas.read_csv(file, sep="\t", comment="#", low_memory=False,
                              names=["seqname", "source", "feature", "start", "end",
                                     "score", "strand", "frame", "attribute"])
        file.close()
        # Parse attributes
        if (parse):
            gff = gff.gff.parse_attributes()
        else:
            gff["id"] = None
            gff["parent"] = None
            gff["name"] = None
            gff["rank"] = None
    elif (".csv" in gff_file):
        file = gzip.open(gff_file, 'r')
        gff = pandas.read_csv(file, sep="\t", keep_default_na=False,na_values=['NaN'])
        file.close()
    return(gff)

def write_gff2csv(gff, filename):
    gff.to_csv(filename, sep="\t", compression="gzip", index=False)

# class gff(DataFrame):
#     """
#     Read a gff file and return a pandas 'DataFrame' object
#     Filename is given in argument
#     :param parse: False, if True, then attributes are parsed and the relationships parent-children are found between gene-mRNA-exon-CDS
#     """
#     def __init__(self, gff_file, parse=True):
#         file = gzip.open(gff_file, 'r')
#         super(gff, self).__init__(pandas.read_csv(file, sep="\t", comment="#", low_memory=False,
#                            names=["seqname", "source", "feature", "start", "end",
#                                   "score", "strand", "frame", "attribute"]))
#         file.close()
#         # Parse attributes
#         if (parse):
#             self = self.parse_attributes()
#         else:
#             self["id"] = None
#             self["parent"] = None
#             self["name"] = None
#             self["rank"] = None
#
#
#     def gff2pandas(self):
#         """
#         Display the gff as an original Pandas DataFrame
#         """
#         frame = {"seqname": self["seqname"], "source": self["source"], "feature": self["feature"],
#                 "start": self["start"], "end": self["end"], "score": self["score"], "strand": self["strand"],
#                 "frame": self["frame"], "attribute": self["attribute"], "id": self["id"],
#                  "parent": self["parent"], "name": self["name"], "rank": self["rank"]}
#         df = pandas.DataFrame(frame)
#         return(df)
#
#     def parse_attributes(self, infer_rank=False):
#         """
#         Parse the column attributes of a gff
#         :param gff: gff, a gff file based on Pandas DataFrame
#         :return: gff, a gff augmented with three columns for attibutes
#         """
#         gff = self
#         attr = self['attribute']
#         # Parse first the available information in the attribute field
#         id = []
#         parent = []
#         name = []
#         for i, a in attr.iteritems():
#             try:
#                 id_name = re.findall("ID=[A-Za-z0-9\\.\\-\\|\\_]*;", a)[0]
#             except:
#                 id_name = ""
#             id_name = id_name.replace(";", "")
#             id_name = id_name.replace("ID=", "")
#             try:
#                 id.append(id_name)
#             except:
#                 id.append(None)
#
#             try:
#                 parent_term = re.findall(";Parent=[A-Za-z0-9\\.\\-\\|\\_]*;", a)[0]
#             except:
#                 parent_term = ""
#             parent_term = parent_term.replace(";", "")
#             parent_term = parent_term.replace("Parent=", "")
#             try:
#                 parent.append(parent_term)
#             except:
#                 parent.append(None)
#
#             try:
#                 name_term = re.findall(";Name=[A-Za-z0-9\\.\\-\\|\\_]*[;]*", a)[0]
#             except:
#                 name_term = ""
#             name_term = name_term.replace(";", "")
#             name_term = name_term.replace("Name=", "")
#             try:
#                 name.append(name_term)
#             except:
#                 name.append(None)
#         gff["id"] = id
#         gff["parent"] = parent
#         gff["name"] = name
#
#         # Infer exon/CDS rank from position or parent
#         # CDS inherits rank of its parent exon
#         rank = [""] * gff.shape[0]
#         if (infer_rank):
#             # TODO optim at this step: long time for iterations
#             for i, a in gff.iterrows():
#                 if (a["feature"] in ["exon", "CDS"]):
#                     set = gff.children(a["parent"])
#                     set = gff.feature(a["feature"])
#                     if (a["strand"] == "+"):
#                         # Rank: how many sequences of the same feature and parent are before that one?
#                         rk = sum(bool(x) for x in [s <= int(a["start"]) for s in list(set["start"])])
#                     elif (a["strand"] == "-"):
#                         # Rank: how many sequences of the same feature and parent are after that one?
#                         # Read in the opposite direction
#                         rk = sum(bool(x) for x in [s >= int(a["start"]) for s in list(set["start"])])
#                 else:
#                     rk = None
#                 rank[i] = rk
#
#         gff["rank"] = rank
#         return(gff)
#
#
#     def rows(self, idx):
#         """
#         Return a new gff object with only rows in the index range
#         """
#         if (type(idx) == str):
#             idx = [idx]
#         subset = self.iloc[idx]
#         return(subset)
#
#     def region(self, start, end, seq):
#         """
#         Return a new gff object with only features (rows) completely within the queried genomic region
#         """
#         subset = self[(self["start"] >= start) & (self["end"] >= end) & (self["seqname"] == seq)]
#         return(subset)
#
#
#     def parent(self, id):
#         """
#         Return a new gff object with the parents of given children ids
#         """
#         if (type(id) == str):
#             id = [id]
#         parent = list(set(self[self["id"].isin(id)]["parent"]))
#         subset = self[self["id"].isin(parent)]
#         return(subset)
#
#
#     def children(self, id):
#         """
#         Return a new gff object with children of given parent ids
#         """
#         if (type(id) == str):
#             id = [id]
#         subset = self[self["parent"].isin(id)]
#         return(subset)
#
#
#     def feature(self, feature):
#         """
#         Return a new gff object with only features of a given feature or list of features
#         """
#         if (type(feature) == str):
#             feature = [feature]
#         subset = self[self["feature"].isin(feature)]
#         return(subset)
#
#
#     def summary(self):
#         """
#         Print summary statistics on the gff (e.g. count of features, number of chromosomes and chromosomes length)
#         """









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