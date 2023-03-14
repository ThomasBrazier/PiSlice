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
    pos = vcf["variants/POS"]

    A_ref = pos[vcf["variants/REF"] == "A"]
    T_ref = pos[vcf["variants/REF"] == "T"]
    A_alt = pos[["A" in i for i in vcf["variants/ALT"]]]
    T_alt = pos[["T" in i for i in vcf["variants/ALT"]]]
    AT_pos = list(set(A_ref) & set(T_alt)) + list(set(T_ref) & set(A_alt))
    AT_pos.sort()

    pos = pos[vcf["variants/CHROM"] == chrom]

    tmp = list(set(pos) & set(AT_pos))
    tmp.sort()

    tmp = tmp[(tmp >= start) & (tmp <= stop)]
    snpcount = len(tmp)
    return(snpcount)


def snp_count_gc(vcf, chrom, start, stop):
    """
    Count the number of GC SNPs from a vcf file within a specified interval start:stop
    Use the scikit-allel package and vcf format
    :vcf: a sckit-allel vcf format
    :chrom: string, the chromosome name
    :start: int, start position, +1 index
    :stop: int, stop position, +1 index
    :return: Number of AT SNPs
    """
    pos = vcf["variants/POS"]

    G_ref = pos[vcf["variants/REF"] == "G"]
    C_ref = pos[vcf["variants/REF"] == "C"]
    G_alt = pos[["G" in i for i in vcf["variants/ALT"]]]
    C_alt = pos[["C" in i for i in vcf["variants/ALT"]]]
    GC_pos = list(set(G_ref) & set(C_alt)) + list(set(C_ref) & set(G_alt))
    GC_pos.sort()

    pos = pos[vcf["variants/CHROM"] == chrom]

    tmp = list(set(pos) & set(GC_pos))
    tmp.sort()

    tmp = tmp[(tmp >= start) & (tmp <= stop)]
    snpcount = len(tmp)
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

