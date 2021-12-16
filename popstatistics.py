# -*- coding: utf-8 -*-
"""
Population genomics statistics.

Functions in this module are used to estimate population genomics statistics along a sequence.
"""

import pandas as pd
from Bio.Seq import Seq
import input
from itertools import compress
import numpy as np

def piSlice(windows, statistics=[""], *args, **kwargs):
    """
    The main function to return a data frame of population genomics statistics for a list of genomic windows.
    :param windows: DataFrame, a 3 columns pandas data frame with chromosome, start, end
    :param statistics: str, a list of statistics to compute
    :param **fasta: fasta, a fasta object with multiple fasta sequences
    :param **gff: DataFrame, a gff object
    :param **vcf: vcf, a vcf object
    :return: DataFrame, a data frame with population statistics for each windows
    """
    fasta = kwargs.get("fasta", "")
    gff = kwargs.get("gff", "")
    vcf = kwargs.get("vcf", "")
    # TODO A progress bar
    # Header
    print("Number of windows:", len(windows.index))
    print("Chromosomes are", " ".join(windows.Chromosome.unique()))

    if "gene_count" in statistics:
        print("Process number of genes")
        estimates = windows.apply(lambda x: gene_count(gff,
                                             x["Chromosome"],
                                             x["Start"],
                                             x["End"]),
                             axis=1)
        windows["gene_count"] = estimates

    if "snp_count" in statistics:
        print("Process number of SNPs")
        estimates = windows.apply(lambda x: snp_count(vcf,
                                             x["Chromosome"],
                                             x["Start"],
                                             x["End"]),
                             axis=1)
        windows["snp_count"] = estimates

    if "gc" in statistics:
        print("Process GC content")
        # Sample sequences
        # Sample all sequences from chromosomes and start-end positions
        list_seq = list(windows.apply(lambda x: fasta.sample_sequence(x["Chromosome"], x["Start"], x["End"]), axis=1))
        # Take care of short sequences (< 6bp) that introduce errors below
        # TODO must return a windows data frame of same size as input, with NA values
        length_seq = list(map(lambda x: len(x), list_seq))
        list_seq = list(compress(list_seq, list(map(lambda x: int(x) > 6, length_seq))))
        windows = windows[list(map(lambda x: int(x) > 6, length_seq))]

        # Compute GC content
        estimates = list(map(lambda x: gc(x), list_seq))
        # Add column for statistics
        windows["gc"] = estimates


    if "gc_cds" in statistics:
        print("Process GC content with codon positions")
        # Compute GC content
        estimates = windows.apply(lambda x: gc_cds(fasta,
                                             gff,
                                             x["Chromosome"],
                                             x["Start"],
                                             x["End"]),
                             axis=1)
        list_gc = [item[0] for item in estimates]
        list_gc1 = [item[1] for item in estimates]
        list_gc2 = [item[2] for item in estimates]
        list_gc3 = [item[3] for item in estimates]
        # Add column for statistics
        windows["gc_cds"] = list_gc
        windows["gc1"] = list_gc1
        windows["gc2"] = list_gc2
        windows["gc3"] = list_gc3

    return windows

def gene_count(gff, chromosome, start, end):
    """
    Count the number of genes beginning in the window (number of start positions).
    :param gff: DataFrame, a gff file with gene annotations
    :param chromosome: str, Chromosome name
    :param start: int, Start position of the sequence
    :param end: int, End position of the sequence
    :return: int, number of genes in the gff window
    """
    gene_count = gff[(gff['seqname'] == str(chromosome)) &
               (gff['start'] >= int(start)) &
               (gff['start'] < int(end)) &
               (gff['feature'] == "gene")]
    gene_count = len(gene_count)

    return gene_count



def snp_count(vcf, chromosome, start, end):
    """
    Count the number of snps in the window.
    :param vcf: vcf, a vcf file with SNPs and their genomic positions
    :param chromosome: str, Chromosome name
    :param start: int, Start position of the sequence
    :param end: int, End position of the sequence
    :return: int, number of snps in the vcf window
    """
    snp_count = vcf.sample_variant(str(chromosome), int(start), int(end))
    snp_count = sum(1 for item in snp_count)

    return snp_count



def gc(sequence):
    """
    Estimate the fraction of G+C bases in a DNA sequence.
    It reads a DNA sequence and count the number of G+C bases divided by the total number of bases.
    GC = (G+C)/(G+C+A+T)
    :param sequence: str, A string containing a DNA sequence
    :return: float, Numeric value of the GC proportion in the sequence
    """
    # Make sequence uppercase for simple computation
    sequence = sequence.upper()
    base_a = sequence.count("A")
    base_c = sequence.count("C")
    base_g = sequence.count("G")
    base_t = sequence.count("T")
    try:
        gc_content = (base_g + base_c)/(base_a + base_c + base_g + base_t)
    except ZeroDivisionError:
        gc_content = np.NaN
    # Do not use the GC calculation from Biopython
    # Because it does not deal with 'N' nucleotides
    # gc_content = GC(sequence)/100
    return gc_content


# TODO GC exact computation to account for ambiguous nucleotides S(G or C)


# TODO Test gc_cds for GC1, GC2, GC3 contents
def gc_cds(fasta, gff, chromosome, start, end):
    """
    Estimate the fraction of G+C bases within CDS at codon positions 1, 2 and 3.
    Use a list of CDS features (start, end, frame, phase) to subset a list of DNA sequences
    and estimate GC content at each position.
    :param fasta: str, A fasta object with the same coordinates as the gff
    :param gff: DataFrame, A gff data frame
    :param chromosome: str, Chromosome name
    :param start: int, Start position of the sequence
    :param end: int, End position of the sequence
    :return: Numeric values of the global GC proportion in the sequence and
    GC proportion at each codon position in the sequence
    """
    # TODO - gc_cds return GC for one sequence only
    # TODO - and use function piSlice to treat multiple sequences
    # GC of the full CDS is the GC content of all CDS regions, without subsetting by rank
    gc_content = gc_cds_rank(fasta, gff, chromosome, start, end, rank="all")
    return gc_content


def gc_cds_rank(fasta, gff, chromosome, start, end, rank=1):
    """
    Estimate the fraction of G+C bases within CDS of a given rank (i.e. position) in the gene
    at codon positions 1, 2 and 3.
    Use a list of CDS features (start, end, frame, phase) and a rank number to subset a list of DNA sequences
    and estimate GC content at each position.

    This function computes directly GC content of a given rank.
    An alternative strategy would be to subset the gff object by rank before computing with gc_cds() directly.

    :param fasta: str, A fasta object with the same coordinates as the gff
    :param gff: DataFrame, A gff data frame
    :param chromosome: str, Chromosome name
    :param start: int, Start position of the sequence
    :param end: int, End position of the sequence
    :param rank: int, rank (position) of the CDS in the gene
    :return: Numeric values of the global GC proportion in the sequence and
    GC proportion at each codon position in the sequence
    """
    # Subset features
    # exons contain UTR that can alter the frame shift
    # It is preferable to estimate GC content on CDS
    if rank == "all":
        feat = gff[(gff['seqname'] == str(chromosome)) &
                   (gff['start'] >= int(start)) &
                   (gff['end'] <= int(end)) &
                   (gff['feature'] == "CDS")]
    else:
        feat = gff[(gff['seqname'] == str(chromosome)) &
                   (gff['start'] >= int(start)) &
                   (gff['end'] <= int(end)) &
                   (gff['feature'] == "CDS") &
                   (gff['rank'] == int(rank))]
    # Sample sequences
    # Sample all sequences from chromosomes and start-end positions
    # Subset a list of DNA sequences according to features positions
    # TODO verify the indexing 0-1 offset - A priori we are not good for maize, GC3 lower than GC1
    list_seq = list(feat.apply(lambda x: fasta.sample_sequence(x["seqname"], x["start"], x["end"]), axis=1))

    # Take care of short sequences (< 6bp) that introduce errors below
    length_seq = list(feat.apply(lambda x: int(x['end']) - int(x['start']), axis=1))
    feat = feat[list(map(lambda x: int(x) > 6, length_seq))]
    list_seq = list(compress(list_seq, list(map(lambda x: int(x) > 6, length_seq))))

    # debug purpose
    # Verify that 1rst exon begins by start codon after reverse complement and frame shift
    # list(map(lambda x: x[0:3], list_seq))
    # Ok, it does verify

    # Strand of the feature
    # Reverse the DNA sequence if strand == "-"
    strand = list(feat.apply(lambda x: x['strand'], axis=1))
    for i, seq in enumerate(list_seq):
        if strand[i] == "-":
            list_seq[i] = seq[::-1]
            # list_seq[i] = str(Seq(seq).reverse_complement())

    # Phase of CDS features
    # Remove 0, 1 or 2 bp at the beginning
    frame = list(feat.apply(lambda x: x['frame'], axis=1))
    for i, seq in enumerate(list_seq):
        list_seq[i] = seq[int(frame[i])::]

    # Split in three vectors of codon position
    codons = "".join(map(lambda x: x[::], list_seq))
    codon1 = "".join(map(lambda x: x[0::3], list_seq))
    codon2 = "".join(map(lambda x: x[1::3], list_seq))
    codon3 = "".join(map(lambda x: x[2::3], list_seq))
    # Estimate GC content at each codon position
    gc123 = gc(codons)
    gc1 = gc(codon1)
    gc2 = gc(codon2)
    gc3 = gc(codon3)
    gc_content = (gc123, gc1, gc2, gc3)
    return gc_content



def gc1(fasta, gff, chromosome, start, end):
    gc1 = gc_cds(fasta, gff, chromosome, start, end)[1]
    return gc1


def gc2(fasta, gff, chromosome, start, end):
    gc2 =  gc_cds(fasta, gff, chromosome, start, end)[2]
    return gc2


def gc3(fasta, gff, chromosome, start, end):
    gc = gc_cds(fasta, gff, chromosome, start, end)[3]
    return gc3




def pi(polymorphism):
    """
    Compute the Nucleotide diversity Pi at a given site (window) in a population, as described by Nei and Li in 1979.
    :param polymorphim:
    :return:

    Reference:
    Nei, M.; Masatoshi Nei; Wen-Hsiung Li (October 1, 1979).
    "Mathematical Model for Studying Genetic Variation in Terms of Restriction Endonucleases". PNAS. 76 (10): 5269–73.
    """
    polymorphism

    return pi