# -*- coding: utf-8 -*-
"""
Population genomics statistics.

Functions in this module are used to estimate population genomics statistics along a sequence.
"""

import pandas as pd


def piSlice(fasta, gff, windows, statistics):
    """
    The main function to return a data frame of population genomics statistics for a list of genomic windows.
    :param fasta: fasta, a fasta object with multiple fasta sequences
    :param gff: gff, a gff object
    :param windows: DataFrame, a 3 columns pandas data frame with chromosome, start, end
    :param statistics: str, a list of statistics to compute
    :return: DataFrame, a data frame with population statistics for each windows
    """
    # Create pandas data frame from the windows input
    popData = windows
    if "gcpos" in statistics:
        estimates = list(popData.apply(lambda x: gcpos(fasta.sample_chromosome(x["Chromosome"]),
                                           gff,
                                           x["Chromosome"],
                                           int(x["Start"]),
                                           int(x["End"])), axis=1))
        list_gc = [item["gc"] for item in estimates]
        list_gc1 = [item["gc1"] for item in estimates]
        list_gc2 = [item["gc2"] for item in estimates]
        list_gc3 = [item["gc3"] for item in estimates]
        # Add column for statistics
        popData["gc"] = list_gc
        popData["gc1"] = list_gc1
        popData["gc2"] = list_gc2
        popData["gc3"] = list_gc3

    return popData


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
    gc_content = (base_g + base_c)/(base_a + base_c + base_g + base_t)

    return gc_content


# TODO GC exact computation


# TODO Test GCpos for GC1, GC2, GC3 contents
def gcpos(sequence, features, chromosome, start, end):
    """
    Estimate the fraction of G+C bases within CDS at codon positions 1, 2 and 3.
    Use a list of CDS features (start, end, frame, phase) to subset a list of DNA sequences
    and estimate GC content at each position.
    :param sequence: str, A string containing the complete DNA sequence with the same coordinates as the gff
    :param features: DataFrame, A gff data frame
    :param chromosome: str, Chromosome name
    :param start: int, Start position of the sequence
    :param end: int, End position of the sequence
    :return: Numeric values of the global GC proportion in the sequence and
    GC proportion at each codon position in the sequence
    """
    # Subset features
    # exons contain UTR that can alter the frame shift
    # It is preferable to estimate GC content on CDS
    feat = features[(features['seqname'] == str(chromosome)) &
                    (features['start'] >= start) &
                    (features['end'] <= end) &
                    (features['feature'] == "CDS")]

    # Subset a list of DNA sequences according to features positions
    # Beware of bp offset, GFF begins numbering at 1 (1-base offset) while Python begins at 0
    # TODO verify the indexing 0-1 offset - A priori we are good, GC3 higher in rice, GC2 lower, as expected
    list_seq = list(feat.apply(lambda x: sequence[x['start'] - 1:x['end']], axis=1))

    # Strand of the feature
    # Reverse str if strand == "-"
    strand = list(feat.apply(lambda x: x['strand'], axis=1))
    for i, seq in enumerate(list_seq):
        if strand[i] == "-":
            list_seq[i] = seq[::-1]

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
    gc_content = {'gc':gc123, 'gc1':gc1, 'gc2':gc2, 'gc3':gc3}

    return gc_content



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