#0 -*- coding: utf-8 -*-
"""
Population genomics statistics.

Functions in this module are used to estimate population genomics statistics along a sequence.
"""

import pandas as pd
import numpy as np
import re
import intervaltree
import warnings

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

# DONE Factorize gene_length, exon_length, intron_length to a generic feature_length function
def feature_length(gff, chromosome, start, end, feature="gene"):
    """
    Estimate the mean length (bp) of a feature in the window
    :param gff: DataFrame, a gff file with gene annotations
    :param chromosome: str, Chromosome name
    :param start: int, Start position of the sequence
    :param end: int, End position of the sequence
    :param feature: str, the type of feature to sample
    :return: int, mean feature length
    """
    feat_count = gff[(gff['seqname'] == str(chromosome)) &
                     (gff['start'] >= int(start)) &
                     (gff['start'] < int(end)) &
                     (gff['feature'] == feature)]
    feat_len = feat_count['end'] - feat_count['start'] + 1
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        feat_len = np.mean(feat_len)
    return feat_len


def max_rank(gff, gene_id):
    """
    Return the max rank for a given gene id. Gff must be parsed for ranks
    """
    # Get second order children (mRNA and exons)
    children = gff.gff.children(gene_id, all=True)
    #children2 = gff.gff.children(children1["id"])
    #frames = [children1, children2]
    #result = pd.concat(frames)
    max_rank = np.max(children["rank"])
    return(max_rank)


def gene_nbexons(gff, chromosome, start, end):
    """
    Estimate the mean number of exons in genes
    :param gff: DataFrame, a gff file with gene annotations, must be parsed before
    :param chromosome: str, Chromosome name
    :param start: int, Start position of the sequence
    :param end: int, End position of the sequence
    :return: int, mean number of exons per gene
    """
    genes = gff[(gff['seqname'] == str(chromosome)) &
               (gff['start'] >= int(start)) &
               (gff['start'] < int(end))].copy()
    # TODO Parse only if ranks have not been inferred
    # genes = genes.gff.parse_attributes(infer_rank=True, verbose=False)
    # Max rank for each gene
    list_genes = genes['id'][genes['feature'] == "gene"]
    gene_nbexons = list(list_genes.apply(lambda x: max_rank(genes, x)))
    # mean = 0 if ranks have not been inferred before
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        gene_nbexons = np.mean(gene_nbexons)
    return(gene_nbexons)

def gene_density(gff, chromosome, start, end):
    """
    Estimate gene density in the window (between 0 and 1)
    :param gff: DataFrame, a gff file with gene annotations
    :param chromosome: str, Chromosome name
    :param start: int, Start position of the sequence
    :param end: int, End position of the sequence
    :return: int, gene density
    """
    gene_count = gff[(gff['seqname'] == str(chromosome)) &
               (gff['start'] >= int(start)) &
               (gff['start'] < int(end)) &
               (gff['feature'] == "gene")]
    gene_len = gene_count['end'] - gene_count['start'] + 1
    gene_len = np.sum(gene_len)
    gene_density = gene_len/(end - start + 1)
    return gene_density



def gc(sequence, min_bp=6):
    """
    Estimate the fraction of G+C bases in a DNA sequence.
    It reads a DNA sequence and count the number of G+C bases divided by the total number of bases.
    GC = (G+C)/(G+C+A+T)
    :param sequence: str, A string containing a DNA sequence
    :return: float, Numeric value of the GC proportion in the sequence
    """
    if isinstance(sequence, str):
        # sequence = str(sequence)
        if len(sequence) > min_bp:
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
        else:
            gc_content = np.NaN
        return gc_content
    else:
        return np.NaN


def gc_count(sequence):
    """
    Count the number of G+C nucleotides in a DNA sequence.
    :param sequence: str, A string containing a DNA sequence
    """
    if isinstance(sequence, str):
        sequence = sequence.upper()
        base_c = sequence.count("C")
        base_g = sequence.count("G")
        gc_count = base_g + base_c
        return gc_count
    else:
        return np.NaN

def at_count(sequence):
    """
    Count the number of A+T nucleotides in a DNA sequence.
    :param sequence: str, A string containing a DNA sequence
    """
    if isinstance(sequence, str):
        sequence = sequence.upper()
        base_a = sequence.count("A")
        base_t = sequence.count("T")
        at_count = base_a + base_t
        return at_count
    else:
        return np.NaN

def missing_nucleotide(sequence):
    """
    Count the number of missing nucleotides in a DNA sequence.
    :param sequence: str, A string containing a DNA sequence
    """
    if isinstance(sequence, str):
        sequence = sequence.upper()
        missing_nuc = sequence.count("N")
        return missing_nuc
    else:
        return np.NaN


# TODO GC exact computation to account for ambiguous nucleotides S(G or C)


# TODO Test gc_cds for GC1, GC2, GC3 contents
def gc_codon(fasta, gff, chromosome, start, end, min_bp=6):
    """
    Estimate the fraction of G+C bases within CDS at codon positions 1, 2 and 3.
    Use a list of CDS features (start, end, frame, phase) to subset a list of DNA sequences
    and estimate GC content at each position.
    :param fasta: str, A fasta object with the same coordinates as the gff
    :param gff: DataFrame, A gff data frame
    :param chromosome: str, Chromosome name
    :param start: int, Start position of the sequence
    :param end: int, End position of the sequence
    :param min_bp: int, the minimal number of nucleotides to consider a sequence
    :return: Numeric values of the global GC proportion in the sequence and
    GC proportion at each codon position in the sequence
    """
    # Subset features
    # exons contain UTR that can alter the frame shift
    # It is preferable to estimate GC content on CDS
    feat = gff[(gff['seqname'] == str(chromosome)) &
               (gff['start'] >= int(start)) &
               (gff['end'] <= int(end)) &
               (gff['feature'] == "CDS")]

    if (feat.shape[0] > 0):
        # Sample all sequences from chromosomes and start-end positions
        # Subset a list of DNA sequences according to features positions
        list_seq = list(feat.apply(lambda x: fasta.sample_sequence(x["seqname"], x["start"], x["end"]), axis=1))
        # Take care of short sequences (typically < 6bp) that introduce errors below
        # Remove sequences shorter than the required number of nucleotides
        # list_seq = list(map(lambda x: x.upper(), list_seq))
        length_seq = list(map(lambda x: len(re.findall("[ATCGatcg]", str(x))), list_seq))
        # Reduce the dataset
        feat = feat.loc[list(map(lambda x: int(x) > min_bp, length_seq))]
        list_seq = list(feat.apply(lambda x: fasta.sample_sequence(x["seqname"], x["start"], x["end"]), axis=1))

        if (feat.shape[0] > 0):
            # Merge overlapping coordinates to estimate CDS proportion properly
            # in case of splicing variants and overlapping sequences
            # Refactoring: use sample_sequence_masked to get CDS regions masked then p = 1 - q/l
            # where p = proportion of CDS and q = length of sequence after masking CDS and l = total length of sequence
            # Masking regions
            mask = [(x, y) for x, y in zip(list(feat.start), list(feat.end))]
            # Sample sequences
            seq = fasta.sample_sequence_masked(chromosome, start, end, mask)
            try:
                cds_proportion = 1 - abs(len(seq) / (end - start + 1))
            except ZeroDivisionError:
                cds_proportion = np.NaN

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
            gc123 = gc(codons, min_bp=min_bp)
            gc1 = gc(codon1, min_bp=min_bp)
            gc2 = gc(codon2, min_bp=min_bp)
            gc3 = gc(codon3, min_bp=min_bp)
        else:
            gc123 = np.NaN
            gc1 = np.NaN
            gc2 = np.NaN
            gc3 = np.NaN
            cds_proportion = np.NaN
    else:
        gc123 = np.NaN
        gc1 = np.NaN
        gc2 = np.NaN
        gc3 = np.NaN
        cds_proportion = np.NaN

    gc_content = (gc123, gc1, gc2, gc3, cds_proportion)
    return gc_content


def gc_noncoding(fasta, gff, chromosome, start, end, min_bp=6):
    """
    Estimate the fraction of G+C bases within non-coding sequences.
    Use a list of CDS features (start, end, frame, phase) to subset a list of non-coding DNA sequences
    :param fasta: str, A fasta object with the same coordinates as the gff
    :param gff: DataFrame, A gff data frame
    :param chromosome: str, Chromosome name
    :param start: int, Start position of the sequence
    :param end: int, End position of the sequence
    :param min_bp: int, the minimal number of nucleotides to consider a sequence
    :return: int, a tuple with the GC content in non-coding sequences and the proportion of non-coding sequence in the window
    """
    feat = gff[(gff['seqname'] == str(chromosome)) &
               (gff['start'] >= int(start)) &
               (gff['end'] <= int(end)) &
               (gff['feature'] == "CDS")]
    if (feat.shape[0] == 0):
        noncoding_seq = fasta.sample_sequence(chromosome, start, end)
        noncoding_prop = 1
        gc_noncoding = gc(noncoding_seq, min_bp=min_bp)
    elif (feat.shape[0] > 0):
        # Masking regions
        mask = [(x,y) for x,y in zip(list(feat.start), list(feat.end))]
        # Sample sequences
        seq = fasta.sample_sequence_masked(chromosome, start, end, mask)
        gc_noncoding = gc(seq)
        try:
            noncoding_prop = len(seq)/(end-start)
        except ZeroDivisionError:
            noncoding_prop = np.NaN
    else:
        gc_noncoding = np.NaN
        noncoding_prop = np.NaN
    return (gc_noncoding, noncoding_prop)


def gc_intergenic(fasta, gff, chromosome, start, end, min_bp=6):
    """
    Estimate the fraction of G+C bases within intergenic sequences.
    Use a list of gene features (start, end) to subset a list of intergenic DNA sequences
    :param fasta: str, A fasta object with the same coordinates as the gff
    :param gff: DataFrame, A gff data frame
    :param chromosome: str, Chromosome name
    :param start: int, Start position of the sequence
    :param end: int, End position of the sequence
    :param min_bp: int, the minimal number of nucleotides to consider a sequence
    :return: int, a tuple with the GC content in intergenic sequences and the proportion of intergenic sequence in the window
    """
    feat = gff[(gff['seqname'] == str(chromosome)) &
               (gff['start'] >= int(start)) &
               (gff['end'] <= int(end)) &
               (gff['feature'] == "gene")]
    if (feat.shape[0] == 0):
        noncoding_seq = fasta.sample_sequence(chromosome, start, end)
        intergenic_prop = 1
        gc_intergenic = gc(noncoding_seq, min_bp=min_bp)
    elif (feat.shape[0] > 0):
        # Masking regions
        mask = [(x,y) for x,y in zip(list(feat.start), list(feat.end))]
        # Sample sequences
        seq = fasta.sample_sequence_masked(chromosome, start, end, mask)
        gc_intergenic = gc(seq)
        try:
            intergenic_prop = len(seq)/(end-start)
        except ZeroDivisionError:
            intergenic_prop = np.NaN
    else:
        gc_intergenic = np.NaN
        intergenic_prop = np.NaN
    return (gc_intergenic, intergenic_prop)



def gc_intron(fasta, gff, chromosome, start, end, min_bp=6, splicing_strategy="merge"):
    """
    Estimate the fraction of G+C bases within intron sequences.
    Use a list of intron features (start, end) to subset a list of intron DNA sequences
    :param fasta: str, A fasta object with the same coordinates as the gff
    :param gff: DataFrame, A gff data frame
    :param chromosome: str, Chromosome name
    :param start: int, Start position of the sequence
    :param end: int, End position of the sequence
    :param splicing_strategy: int, the minimal number of nucleotides to consider a sequence
    :return: int, a tuple with the GC content in intron sequences and the proportion of intron sequence in the window
    """
    feat = gff[(gff['seqname'] == str(chromosome)) &
               (gff['start'] >= int(start)) &
               (gff['end'] <= int(end)) &
               (gff['feature'] == "intron")]
    if (feat.shape[0] == 0):
        gc_intron = np.NaN
        intron_prop = np.NaN
    # If only one feature, Dataframe is transformed in Series
    elif (isinstance(feat, pd.Series)):
        list_start = [min(feat["start"], feat["end"])]
        list_end = [max(feat["start"], feat["end"])]
        list_seq = [fasta.sample_sequence(chromosome, x, y) for x,y in zip(list_start, list_end)]
        # Sample sequences
        seq = "".join(list_seq)
        gc_intron = gc(seq, min_bp)
        try:
            intron_prop = len(seq)/(end-start)
        except ZeroDivisionError:
            intron_prop = np.NaN
    elif (feat.shape[0] > 0):
        if (splicing_strategy == "merge"):
            list_start = [x[1] for x in feat["start"].items()]
            list_end = [x[1] for x in feat["end"].items()]
            intervals = [(min(x, y), max(x, y)) for x, y in zip(list_start, list_end) if x != y]
            merge_splicing = intervaltree.IntervalTree.from_tuples(intervals)
            list_start = [x.begin for x in merge_splicing]
            list_end = [x.end for x in merge_splicing]
            # if ((len(list_start) > 1) & (len(list_end) > 1)):
            #     intervals = [(x,y) for x,y in zip([list_start], [list_end])]
            #     merge_splicing = intervaltree.IntervalTree.from_tuples(intervals)
            #     list_start = [x.begin for x in merge_splicing]
            #     list_end = [x.end for x in merge_splicing]
            # else:
            #     # Inverse coordinates if sequence is on "-" strand (i.e. start > end)
            #     list_start = min(list_start, list_end)
            #     list_end = min(list_start, list_end)
        list_seq = [fasta.sample_sequence(chromosome, x, y) for x,y in zip(list_start, list_end)]
        # Sample sequences
        seq = "".join(list_seq)
        gc_intron = gc(seq, min_bp)
        try:
            intron_prop = len(seq)/(end-start)
        except ZeroDivisionError:
            intron_prop = np.NaN
    else:
        gc_intron = np.NaN
        intron_prop = np.NaN
    return (gc_intron, intron_prop)



def gc1(fasta, gff, chromosome, start, end):
    gc1 = gc_codon(fasta, gff, chromosome, start, end)[1]
    return gc1


def gc2(fasta, gff, chromosome, start, end):
    gc2 = gc_codon(fasta, gff, chromosome, start, end)[2]
    return gc2


def gc3(fasta, gff, chromosome, start, end):
    gc3 = gc_codon(fasta, gff, chromosome, start, end)[3]
    return gc3

def gc3exon1(fasta, gff, chromosome, start, end, min_bp=6):
    gffexon1 = gff.loc[((gff["rank"] == 0) | (gff["rank"] == 1))]
    gc3exon1 = gc_codon(fasta, gffexon1, chromosome, start, end, min_bp=min_bp)[3]
    return gc3exon1


def cpg(sequence):
    """"
    Estimate the CpG density as the number of CG sites divided by the total number of sites
    (i.e. total number of nucleotides divided by two)
    :param sequence: str, a fasta sequence
    :return: float, a CpG density
    """
    if (isinstance(sequence, str)):
        if len(sequence) > 6:
            sequence = sequence.upper()
            if "CG" in sequence:
                seq_len = len(re.findall("[ATCGatcg]", sequence))
                cpg_density = sequence.count('CG')/(seq_len/2)
            else:
                cpg_density = 0
        else:
            cpg_density = np.NaN
        return(cpg_density)
    else:
        return(np.NaN)
