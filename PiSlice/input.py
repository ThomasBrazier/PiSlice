# -*- coding: utf-8 -*-
"""
Input

Define classes for genomic objects like fasta, gff and vcf

Functions in this module import .fasta, .gff, .vcf and coordinate files in appropriate formats
"""

# TODO Create a global unified and consistent 'genomic' object that contain all data in slots: genome sequence, features, snps and associated metadata
import pandas as pd
from cyvcf2 import VCF
import pandas
import numpy as np
import gzip
import re
import mapply
import multiprocessing
import intervaltree
from multiprocess import Pool
from Bio import SeqIO


# Use read_vcf of scikit-allel and do not develop a new class
# KISS
# TODO Implement vcf import by sgkit
# DEPRECATED
# class vcf(VCF):
#     """
#     Read a vcf file and return a 'vcf' object
#     inheriting the 'VCF' class from 'cyvcf2'
#     Filename is given in argument when creating the 'vcf' object
#     """
#
#     def sample_individual(self, samples):
#         """
#         Sample individuals and modify the 'vcf' object
#         """
#         # TODO copy the vcf instead on changing it directly
#         # Make a list individuals to sample
#         # TODO implement methods to sample individuals by indexes or names
#         # Subset individuals
#         self.set_samples(samples)
#
#     def sample_variant(self, chromosome, start, end):
#         """
#         Sample variants on a chromosome within start-end coordinates (bp)
#          and return a generator object to iterate over the 'vcf' object
#         """
#         # Make a list of loci to sample
#         # TODO
#         # Subset loci
#         snps = self(":".join([str(chromosome), "-".join([str(start), str(end)])]))
#         return snps



# class genotype:
#     """
#     Create 'genotype' object from a 'vcf' object which contains:
#     * gt: genotypes as a 2D numpy array
#     * samples: sample names
#     * loci: locus names
#     """
#     def __init__(self, variants, chromosome, start, end):
#         # Subset a region
#         snps = variants.sample_variant(chromosome, start, end)
#         # Get names of loci
#         loci = []
#         for variant in snps:
#             loci.append(variant.ID)
#         # Set dimensions
#         df = []
#         # For loop, because cannot unpack non-iterable cyvcf2.cyvcf2.Variant object
#         snps = variants.sample_variant(chromosome, start, end)
#         for variant in snps:
#             df.append(list(variant.gt_types))
#         df = np.array(df)
#         self.gt = df
#         self.samples = variants.samples
#         self.loci = loci

# TODO Testing the fasta.fetch() with testing data
# Test Driven Dev
# DONE Replace pysam FastaFile by a dictionnary with (chr name: sequence string)
# Assert if new class reproduce results of the pysam.fetch
# TODO parallelize pool.apply()
class fasta():
    """
    Read a fasta file and return a 'fasta' sequence object
    inheriting the 'FastaFile' class from 'pysam'
    Filename is given in argument when creating the 'fasta' object
    """
    # TODO Constructor for the fasta class
    def __init__(self, filename):
        f = gzip.open(filename, "rt")
        self.seq = {rec.id: str(rec.seq) for rec in SeqIO.parse(f, "fasta")}

    def summary(self):
        """
        Return a summary of the fasta file. Number of sequences, names and sequence length.
        """

    def length(self):
        """
        Return the length of each sequence
        """
        return [len(x) for x in self.sequence()]

    def seqname(self):
        """
        Return sequence names in a fasta object
        """
        return list(self.seq.keys())

    def sequence(self):
        """
        Return sequence names in a fasta object
        """
        return list(self.seq.values())

    def sample_chromosome(self, chromosome):
        """
        Sample a whole chromosome sequence
        :param chromosome: str, the name of the chromosome
        :return: A DNA sequence of the whole chromosome (str)
        """
        # Chromosome can be either an integer (index) or a string (name)
        return self.seq[str(chromosome)]


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
        # Verify that start-end positions are not inverted (start < end)
        start = min(start, end)
        end = max(start, end)
        return self.seq[str(chromosome)][(start - 1):end:]

    def sample_sequence_masked(self, chromosome, start, end, mask):
        """
        Sample a DNA sequence on a chromosome within start-end coordinates (bp)
        on a 1-offset (sequence begins at 1) and mask regions given in 'mask'
        :param chromosome: str, the name of the chromosome
        :param start: int, start position in the sequence
        :param end: int, end position in the sequence
        :param mask: int, a list of tuples giving start-end coordinates of the region to mask
        :return: A DNA sequence within boundaries start-end on a given chromosome (str)
        """
        # Expand masked intervals
        mask = [(x - 1, y + 1) for x,y in mask]
        # Take care of Null interval objects
        coord = intervaltree.IntervalTree.from_tuples([(start, end)])
        [coord.chop(x, y) for x, y in tuple(mask) if x != y]

        seq = [self.sample_sequence(chromosome, x[0], x[1]) for x in coord.items()]

        return seq

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

    def parse_attributes(self, infer_rank=False, parse_introns=False, parse_utr=False, n_cpus=8, verbose=True):
        """
        Parse the column attributes of a gff
        :param gff: gff, a gff file based on Pandas DataFrame
        :return: gff, a gff augmented with three columns for attibutes
        """
        gff_obj = self._obj
        gff_obj["id"] = None
        gff_obj["parent"] = None
        gff_obj["name"] = None
        gff_obj["gene_biotype"] = None
        #gff_obj = gff_obj.reset_index() # Reset index to get continuous row indexes to iterate
        if (n_cpus == 0):
            n_cpus = multiprocessing.cpu_count()
        sensible_cpus = mapply.parallel.sensible_cpu_count()
        # Parse first the available information in the attribute field
        if verbose:
            print("Parsing attributes")

        # TODO Vectorize
        attr = gff_obj["attribute"].astype(str)
        def find_attribute(s, category="Name"):
            """
            Find the attribute with the 'category' in a string
            """
            try:
                attribute = re.findall("[;]*" + category + "=[A-Za-z0-9\\.\\-\\|\\_\\:]*[;]*", s)[0]
            except:
                attribute = ""
            attribute = attribute.replace(";", "")
            attribute = attribute.replace(category + "=", "")
            return(attribute)

        id_attr = list(map(lambda x: find_attribute(x, category="ID"), attr))
        try:
            gff_obj["id"] = id_attr
        except:
            gff_obj["id"] = None
        name_attr = list(map(lambda x: find_attribute(x, category="Name"), attr))
        try:
            gff_obj["name"] = name_attr
        except:
            gff_obj["name"] = None
        parent_attr = list(map(lambda x: find_attribute(x, category="Parent"), attr))
        try:
            gff_obj["parent"] = parent_attr
        except:
            gff_obj["parent"] = None
        biotype_attr = list(map(lambda x: find_attribute(x, category="gene_biotype"), attr))
        try:
            gff_obj["gene_biotype"] = biotype_attr
        except:
            gff_obj["gene_biotype"] = None


        # for i, r in gff_obj.iterrows():
        #     a = r["attribute"]
        #     try:
        #         id_name = re.findall("ID=[A-Za-z0-9\\.\\-\\|\\_\\:]*;", a)[0]
        #     except:
        #         id_name = ""
        #     id_name = id_name.replace(";", "")
        #     id_name = id_name.replace("ID=", "")
        #     try:
        #         gff_obj.loc[i,"id"] = id_name
        #     except:
        #         gff_obj.loc[i,"id"] = None
        #
        #     try:
        #         parent_term = re.findall(";Parent=[A-Za-z0-9\\.\\-\\|\\_\\:]*[;]*", a)[0]
        #     except:
        #         parent_term = ""
        #     parent_term = parent_term.replace(";", "")
        #     parent_term = parent_term.replace("Parent=", "")
        #     try:
        #         gff_obj.loc[i,"parent"] = parent_term
        #     except:
        #         gff_obj.loc[i,"parent"] = None
        #
        #     try:
        #         name_term = re.findall(";Name=[A-Za-z0-9\\.\\-\\|\\_\\:]*[;]*", a)[0]
        #     except:
        #         name_term = ""
        #     name_term = name_term.replace(";", "")
        #     name_term = name_term.replace("Name=", "")
        #     try:
        #         gff_obj.loc[i,"name"] = name_term
        #     except:
        #         gff_obj.loc[i,"name"] = None
        #
        #     try:
        #         biotype_term = re.findall(";gene_biotype=[A-Za-z0-9\\.\\-\\|\\_]*[;]*", a)[0]
        #     except:
        #         biotype_term = ""
        #     biotype_term = biotype_term.replace(";", "")
        #     biotype_term = biotype_term.replace("gene_biotype=", "")
        #     try:
        #         gff_obj.loc[i,"gene_biotype"] = biotype_term
        #     except:
        #         gff_obj.loc[i,"gene_biotype"] = None

        # Change "transcript" feature to "mRNA" -> consistency
        gff_obj.loc[gff_obj["feature"] == "transcript", "feature"] = "mRNA"


        if (parse_introns):
            if verbose:
                print("Parsing introns")
            # Infer introns features
            # Introns are sequences between two consecutive exons within the same gene
            # They have a rank to infer subsequently
            def intron(gff, gene_id):
                """
                Return a gff data frame with intron features of a given gene
                :param gff: gff, the original gff dataframe
                :param gene_id: str, a gen id
                :return: pandas, a gff dataframe with introns features
                """
                subset = gff.gff.children(gene_id).copy()
                subset = subset.gff.feature("exon").copy()
                if (type(subset) == pd.core.frame.DataFrame):
                    gff_introns = subset[:-1].copy()
                    intron_start = list(subset["end"])
                    intron_end = list(subset["start"])
                    intron_start.sort()
                    intron_end.sort()
                    del intron_start[-1]
                    del intron_end[0]
                    intron_start = [x + 1 for x in intron_start]
                    intron_end = [x + 1 for x in intron_end]
                    # intron_start = subset["end"].sort_values()
                    # intron_start = intron_start.drop(index=intron_start.index[-1])
                    # intron_end = subset["start"].sort_values()
                    # intron_end = intron_end.drop(index=intron_end.index[0])
                    # gff_introns = subset.copy(deep=True)
                    # gff_introns = gff_introns.drop(index=gff_introns.index[-1])
                    gff_introns["start"] = intron_start
                    gff_introns["end"] = intron_end
                    gff_introns = gff_introns.replace(['exon'], 'intron')
                    gff_introns["id"] = gff_introns["id"].str.replace("exon", "intron", regex=True)
                    return (gff_introns)
                else:
                    pass

            p = Pool(min(sensible_cpus, n_cpus))
            # List of exon parents
            list_parent = list(gff_obj.gff.feature("exon")["parent"].unique())
            #res = list(map(lambda x: intron(gff_obj, x), list_parent))
            res = list(map(lambda x: intron(gff_obj, x), list_parent))
            # Append the gff of new introns to the dataset
            if verbose:
                print("Append new introns")
            gff_obj = gff_obj.append(res)
            #gff_obj = pd.concat([gff_obj, res])

        # Infer exon/CDS rank from position or parent
        # CDS inherits rank of its parent exon
        # rank = [0] * gff_obj.shape[0]
        def rank_inference(gff_obj, x):
            # Optim: do not create new objects
            #x = gff_obj.iloc[idx,:]
            gff_subset = gff_obj[(gff_obj["feature"] == x["feature"]) & (gff_obj["parent"] == x["parent"])].copy()
            #filter_ = ((gff_obj["feature"] == x["feature"]) & (gff_obj["parent"] == x["parent"]))
            #set = gff_obj.gff.children(x["parent"])
            #set = set.gff.feature(x["feature"])
            #gff_obj[(gff_obj["feature"] == x["feature"]) & (gff_obj["parent"] == x["parent"])]["start"]
            if (x["strand"] == "+"):
            # Rank: how many sequences of the same feature and parent are before that one?
                rk = sum(bool(z) for z in [s <= int(x["start"]) for s in list(gff_subset["start"])])
            elif (x["strand"] == "-"):
            # Rank: how many sequences of the same feature and parent are after that one?
            # Read in the opposite direction
                rk = sum(bool(z) for z in [s >= int(x["start"]) for s in list(gff_subset["start"])])
            else:
                rk = 0
            return(rk)

        if (infer_rank):
            if verbose:
                print("Inferring the rank of exons/CDS/introns")
            # TODO vectorization
            mapply.init(n_workers=min(sensible_cpus, n_cpus))
            if verbose:
                rank = gff_obj.mapply(lambda x: rank_inference(gff_obj, x) if x["feature"] in ["exon", "CDS", "intron"] else 0,
                                 axis=1)
            else:
                rank = gff_obj.mapply(
                    lambda x: rank_inference(gff_obj, x) if x["feature"] in ["exon", "CDS", "intron"] else 0,
                    axis=1, progressbar=False)
        else:
            rank = [0] * gff_obj.shape[0]
        gff_obj["rank"] = rank

        if (parse_utr):
            # Infer UTR features
            # UTR sequences are beginning of the first exon or the end of the last exon (exon - CDS)
            # According to gff specifications, "UTRs, splice sites and translational start and stop sites.
            # These are implied by the combination of exon and CDS and do not need to be explicitly annotated as part of the canonical gene"
            if verbose:
                print("Parsing UTR regions")

            def utr_transcript(gff_obj, mrna):
                """
                Return UTR of a single transcript
                :param gff_obj: a GFF where to store data
                :param mrna: id of a gene or mRNA to get a single transcript
                """
                transcript = gff_obj.gff.children(mrna, all=True)
                if ("CDS" not in transcript["feature"].unique()):
                    return (None)
                strand = gff_obj.gff.id(mrna)["strand"]
                if (strand.unique().item() == "-"):
                    # TSS is the end of first CDS
                    tss = max(transcript.loc[transcript["feature"] == "CDS", "end"])
                    tts = min(transcript.loc[transcript["feature"] == "CDS", "start"])
                    # UTR5: all exons ending after TSS, minus first CDS
                    # min start become TSS + 1
                    # UTR3: all exons starting before TTS, minus last CDS
                    #
                    utr5 = transcript.loc[((transcript["end"] > tss) & (transcript["feature"] == "exon")), :].copy()
                    try:
                        utr5.loc[utr5["start"] == min(utr5["start"]), "start"] = tss + 1
                        utr5["feature"] = "utr5"
                    except ValueError:
                        pass
                    utr3 = transcript.loc[((transcript["start"] < tts) & (transcript["feature"] == "exon")),:].copy()
                    try:
                        utr3.loc[utr3["end"] == max(utr3["end"]), "end"] = tts - 1
                        utr3["feature"] = "utr3"
                    except ValueError:
                        pass
                elif (strand.unique().item() == "+"):
                    # TSS is the start of first CDS; TTS is the end of last CDS
                    tss = min(transcript.loc[transcript["feature"] == "CDS", "start"])
                    tts = max(transcript.loc[transcript["feature"] == "CDS", "end"])
                    # UTR5: all exons starting before TSS, max end minus start of first CDS
                    # UTR3: all exons ending after TTS, min start minus end of last CDS
                    utr5 = transcript.loc[((transcript["start"] < tss) & (transcript["feature"] == "exon")),:].copy()
                    try:
                        utr5.loc[utr5["end"] == max(utr5["end"]), "end"] = tss - 1
                        utr5["feature"] = "utr5"
                    except ValueError:
                        pass
                    utr3 = transcript.loc[((transcript["end"] > tts) & (transcript["feature"] == "exon")), :].copy()
                    try:
                        utr3.loc[utr3["start"] == max(utr3["start"]), "start"] = tts + 1
                        utr3["feature"] = "utr3"
                    except ValueError:
                        pass
                #utr = utr5
                #utr = utr.append(utr3)
                utr = pd.concat([utr5, utr3])
                return (utr)

            def utr_parse(gff_obj, gene):
                """
                Parse UTRs of each transcript in a gene and return a list of lists
                """
                children = gff_obj.gff.children(gene, all=True)
                # Take care of non coding genes
                if ("CDS" not in children["feature"].unique()):
                    return(None)
                if ("mRNA" in children["feature"].unique()):
                    ids = children.loc[(children["feature"] == "mRNA"), "id"]
                else:
                    ids = gene
                res = list(map(lambda x: utr_transcript(gff_obj, x), ids))
                clean_res = list(filter(lambda x: x is not None, res))
                try:
                    res_utr = pd.concat(clean_res)
                    res_utr["rank"] = 0
                    return (res_utr)
                except ValueError:
                    pass


            #utrs = list(map(lambda x: utr_parse(gff_obj, x), list_genes))
            #gff_obj = gff_obj.append(utrs)

            # 5'UTR in first exon, 3'UTR in last exon, UTRs inherit attributes from their parent exon
            # Subset first and last exons for each gene/mRNA
            list_genes = list(gff_obj.loc[(gff_obj["feature"] == "gene"), "id"])
            list_genes = list(filter(lambda x: x is not None, list_genes))
            #list_genes = gff_obj.loc[(gff_obj["feature"] == "gene"), "id"]
            #utrs = [utr_parse(gff_obj, x) for x in list_genes]
            p = Pool(min(sensible_cpus, n_cpus))
            utrs = list(p.map(lambda x: utr_parse(gff_obj, x), list_genes))
            #utrs = list(p.map(lambda x: utr_parse(gff_obj, x), list_genes))
            #mapply.init(n_workers=n_cpus)
            #utrs = list_genes.mapply(lambda x: utr_parse(gff_obj, x))
            #utrs = pd.concat(utrs)
            # Clean the list for NoneType
            # TODO better handling of this error: utr_parse return NoneType
            clean_utrs = list(filter(lambda x: x is not None, utrs))
            if verbose:
                print("Append new UTRs")
            #clean_utrs = [x for x in utrs if x is not None]
            gff_obj = gff_obj.append(clean_utrs)
            #clean_utrs = pd.Series(clean_utrs)
            #clean_utrs = pd.DataFrame(clean_utrs)
            #gff_obj = pd.concat([gff_obj, clean_utrs])

        gff_obj.start = gff_obj.start.astype(int, errors='ignore') # Leave NA values as they are
        gff_obj.end = gff_obj.end.astype(int, errors='ignore')
        return(gff_obj)


    def region(self, start, end, seq):
        """
        Return a new gff object with only features (rows) completely within the queried genomic region
        """
        subset = self._obj[(self._obj["start"] >= start) & (self._obj["end"] <= end) & (self._obj["seqname"] == seq)]
        return(subset)

    def id(self, id):
        """
        Return a new gff object with the lines of id
        """
        if (type(id) == str):
            id = [id]
        subset = self._obj[self._obj["id"].isin(id)]
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


    def children(self, id, all=True):
        """
        Return a new gff object with children of given parent ids
        """
        if (type(id) == str):
            id = [id]
        subset = self._obj[self._obj["parent"].isin(id)]
        # Second order children, e.g. exons children of mRNA
        if all:
            recursive_children = self._obj.gff.children(list(subset["id"]), all=False)
            # subset = subset.append(recursive_children) ## Deprecated in pandas since 1.4.0; use concat() instead
            subset = pd.concat([subset, recursive_children])
            # Add one more level if children are found
            # if len(recursive_children) > 0:
            #     recursive_children = self._obj.gff.children(list(subset["id"]), all=False)
            #     subset = subset.append(recursive_children)
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


def read_gff(gff_file, parse=False, parse_introns=False, parse_utr=False, infer_rank=False, n_cpus=8):
    if (".gff" in gff_file):
        file = gzip.open(gff_file, 'r')
        gff = pandas.read_csv(file, sep="\t", comment="#", low_memory=False,
                              names=["seqname", "source", "feature", "start", "end",
                                     "score", "strand", "frame", "attribute"])
        file.close()
        # Parse attributes
        if (parse):
            gff = gff.gff.parse_attributes(parse_introns=parse_introns, parse_utr=parse_utr, infer_rank=infer_rank, n_cpus=n_cpus)
        else:
            gff["id"] = None
            gff["parent"] = None
            gff["name"] = None
            gff["gene_biotype"] = None
            gff["rank"] = None
    elif (".csv" in gff_file):
        file = gzip.open(gff_file, 'r')
        gff = pandas.read_csv(file, sep="\t", keep_default_na=False, na_values=['NaN'])
        gff = gff.astype({"seqname": str})
        file.close()
    return(gff)

def write_gff2csv(gff, filename):
    # Preserve column order
    col = gff.columns
    gff[col].to_csv(filename, sep="\t", compression="gzip", index=False, na_rep="NA", header=True, columns=col)
