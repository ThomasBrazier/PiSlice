#0 -*- coding: utf-8 -*-
"""
Population genomics statistics estimated from a vcf file.
"""

import allel


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
    ac = g.count_alleles()
    pos = vcf["variants/POS"]
    pos = pos[vcf["variants/CHROM"] == chrom]
    pi = allel.sequence_diversity(pos, ac, start=start, stop=stop)
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
    theta_watterson = allel.watterson_theta(pos, ac, start=start, stop=stop)
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
    tajima_d = allel.tajima_d(ac, pos, start=start, stop=stop)
    return(tajima_d)

