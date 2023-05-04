#0 -*- coding: utf-8 -*-
"""
Population genomics statistics estimated from a vcf file.
"""

import allel
import numpy as np


def snp_count(vcf, chrom, start, stop):
    """
    Count the number of SNPs from a vcf file within a specified interval start:stop
    Use the scikit-allel package and vcf format
    :vcf: a sckit-allel vcf format
    :chrom: string, the chromosome name
    :start: int, start position, +1 index
    :stop: int, stop position, +1 index
    :return: snp density (snp/bp)
    """
    pos = vcf["variants/POS"]
    pos = pos[vcf["variants/CHROM"] == chrom]
    pos = pos[(pos >= start) & (pos <= stop)]
    snpcount = len(pos)
    return(snpcount)


def snp_count_combination(vcf, chrom, start, stop, base1, base2):
    """
    Count the number of SNPs corresponding to a given combination of two nucleotides
    from a vcf file within a specified interval start:stop
    Use the scikit-allel package and vcf format
    :vcf: a sckit-allel vcf format
    :chrom: string, the chromosome name
    :start: int, start position, +1 index
    :stop: int, stop position, +1 index
    :base1: a nucleotide, either A, C, T, G
    :base2: the other nucleotide
    :return: Number of AG SNPs
    """
    pos = vcf["variants/POS"]
    pos = pos[vcf["variants/CHROM"] == chrom]

    snp_ref = vcf["variants/REF"]
    snp_ref = snp_ref[vcf["variants/CHROM"] == chrom]
    snp_ref = snp_ref[(pos >= start) & (pos <= stop)]

    snp_alt = vcf["variants/ALT"]
    snp_alt = snp_alt[vcf["variants/CHROM"] == chrom]
    snp_alt = snp_alt[(pos >= start) & (pos <= stop)]

    ref1 = (snp_ref == str(base1))
    ref2 = (snp_ref == str(base2))
    alt1 = [str(base1) in i for i in snp_alt]
    alt2 = [str(base2) in i for i in snp_alt]
    try:
        snpcount = sum(ref1 & alt1) + sum(ref2 & alt2)
    except TypeError:
        snpcount = np.NaN
    return(snpcount)


def snp_count_at(vcf, chrom, start, stop):
    """
    Count the number of AT SNPs from a vcf file within a specified interval start:stop
    Use the scikit-allel package and vcf format
    :vcf: a sckit-allel vcf format
    :chrom: string, the chromosome name
    :start: int, start position, +1 index
    :stop: int, stop position, +1 index
    :return: Number of AT SNPs
    """
    snpcount = snp_count_combination(vcf, chrom, start, stop, "A", "T")
    return(snpcount)


def snp_count_gc(vcf, chrom, start, stop):
    """
    Count the number of GC SNPs from a vcf file within a specified interval start:stop
    Use the scikit-allel package and vcf format
    :vcf: a sckit-allel vcf format
    :chrom: string, the chromosome name
    :start: int, start position, +1 index
    :stop: int, stop position, +1 index
    :return: Number of GC SNPs
    """
    snpcount = snp_count_combination(vcf, chrom, start, stop, "G", "C")
    return(snpcount)

def snp_count_ag(vcf, chrom, start, stop):
    """
    Count the number of AG SNPs from a vcf file within a specified interval start:stop
    Use the scikit-allel package and vcf format
    :vcf: a sckit-allel vcf format
    :chrom: string, the chromosome name
    :start: int, start position, +1 index
    :stop: int, stop position, +1 index
    :return: Number of AG SNPs
    """
    snpcount = snp_count_combination(vcf, chrom, start, stop, "A", "G")
    return(snpcount)

def snp_count_ac(vcf, chrom, start, stop):
    """
    Count the number of AC SNPs from a vcf file within a specified interval start:stop
    Use the scikit-allel package and vcf format
    :vcf: a sckit-allel vcf format
    :chrom: string, the chromosome name
    :start: int, start position, +1 index
    :stop: int, stop position, +1 index
    :return: Number of AC SNPs
    """
    snpcount = snp_count_combination(vcf, chrom, start, stop, "A", "C")
    return(snpcount)

def snp_count_gt(vcf, chrom, start, stop):
    """
    Count the number of GT SNPs from a vcf file within a specified interval start:stop
    Use the scikit-allel package and vcf format
    :vcf: a sckit-allel vcf format
    :chrom: string, the chromosome name
    :start: int, start position, +1 index
    :stop: int, stop position, +1 index
    :return: Number of GT SNPs
    """
    snpcount = snp_count_combination(vcf, chrom, start, stop, "G", "T")
    return(snpcount)


def snp_density(vcf, chrom, start, stop):
    """
    Estimate the SNP density from a vcf file within a specified interval start:stop
    Use the scikit-allel package and vcf format
    :vcf: a sckit-allel vcf format
    :chrom: string, the chromosome name
    :start: int, start position, +1 index
    :stop: int, stop position, +1 index
    :return: snp density (snp/bp)
    """
    snps = snp_count(vcf, chrom, start, stop)
    seqlen = stop - start
    snpdens = snps/seqlen
    return(snpdens)


def pi(vcf, chrom, start, stop):
    """
    Estimate the nucleotide diversity Pi from a vcf file within a specified interval start:stop
    Use the scikit-allel package and vcf format
    Return the raw value of Pi over the genomic interval
    Divide by the length of interval to get Pi per site
    :vcf: a sckit-allel vcf format
    :chrom: string, the chromosome name
    :start: int, start position, +1 index
    :stop: int, stop position, +1 index
    :return: float, pi, nucleotide diversity value
    """
    chrom = str(chrom)
    g = allel.GenotypeArray(vcf["calldata/GT"])
    g = g[vcf["variants/CHROM"] == chrom]
    try:
        ac = g.count_alleles()
        pos = vcf["variants/POS"]
        pos = pos[vcf["variants/CHROM"] == chrom]
        tmp = pos[pos >= start]
        tmp = tmp[tmp <= stop]
        if len(tmp) > 0:
            if start > stop:
                pi = allel.sequence_diversity(pos, ac, start=stop, stop=start)
            else:
                pi = allel.sequence_diversity(pos, ac, start=start, stop=stop)
        else:
            pi = np.NaN
    except ValueError:
        pi = np.NaN

    return(pi)


def theta_watterson(vcf, chrom, start, stop):
    """
    Estimate the Watterson's Theta from a vcf file within a specified interval start:stop
    Use the scikit-allel package and vcf format
    :vcf: a sckit-allel vcf format
    :chrom: string, the chromosome name
    :start: int, start position, +1 index
    :stop: int, stop position, +1 index
    :return: float, theta_watterson, nucleotide diversity value
    """
    chrom = str(chrom)
    g = allel.GenotypeArray(vcf["calldata/GT"])
    g = g[vcf["variants/CHROM"] == chrom]
    ac = g.count_alleles()
    pos = vcf["variants/POS"]
    pos = pos[vcf["variants/CHROM"] == chrom]
    tmp = pos[pos >= start]
    tmp = tmp[tmp <= stop]
    if len(tmp) > 0:
        if start > stop:
            theta_watterson = allel.watterson_theta(pos, ac, start=stop, stop=start)
        else:
            theta_watterson = allel.watterson_theta(pos, ac, start=start, stop=stop)
    else:
        theta_watterson = np.NaN
    return(theta_watterson)

def tajima_d(vcf, chrom, start, stop):
    """
    Estimate the Watterson's Theta from a vcf file within a specified interval start:stop
    Use the scikit-allel package and vcf format
    :vcf: a sckit-allel vcf format
    :chrom: string, the chromosome name
    :start: int, start position, +1 index
    :stop: int, stop position, +1 index
    :return: float, theta_watterson, nucleotide diversity value
    """
    chrom = str(chrom)
    g = allel.GenotypeArray(vcf["calldata/GT"])
    g = g[vcf["variants/CHROM"] == chrom]
    ac = g.count_alleles()
    pos = vcf["variants/POS"]
    pos = pos[vcf["variants/CHROM"] == chrom]
    tmp = pos[pos >= start]
    tmp = tmp[tmp <= stop]
    if len(tmp) > 0:
        if start > stop:
            tajima_d = allel.tajima_d(ac, pos, start=stop, stop=start)
        else:
            tajima_d = allel.tajima_d(ac, pos, start=start, stop=stop)
    else:
        tajima_d = np.NaN
    return(tajima_d)

