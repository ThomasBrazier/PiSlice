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
    :start: int, start position, 0 index
    :stop: int, stop position, 0 index
    :return: snp density (snp/bp)
    """
    pos = vcf["variants/POS"]
    pos = pos[vcf["variants/CHROM"] == chrom]
    pos = pos[pos >= start]
    pos = pos[pos <= stop]
    snpcount = len(pos)
    return(snpcount)


def snp_density(vcf, chrom, start, stop):
    """
    Estimate the SNP density from a vcf file within a specified interval start:stop
    Use the scikit-allel package and vcf format
    :vcf: a sckit-allel vcf format
    :chrom: string, the chromosome name
    :start: int, start position, 0 index
    :stop: int, stop position, 0 index
    :return: snp density (snp/bp)
    """
    snps = snp_count(vcf, chrom, start, stop)
    seqlen = stop - start
    snpdens = snps/seqlen
    return(snpdens)


def pi(vcf, chrom, start, stop):
    """
    Estimate the nucleotide diversity from a vcf file within a specified interval start:stop
    Use the scikit-allel package and vcf format
    :vcf: a sckit-allel vcf format
    :chrom: string, the chromosome name
    :start: int, start position, 0 index
    :stop: int, stop position, 0 index
    :return: pi, nucleotide diversity value
    """
    chrom = str(chrom)
    g = allel.GenotypeArray(vcf["calldata/GT"])
    g = g[vcf["variants/CHROM"] == chrom]
    ac = g.count_alleles()
    pos = vcf["variants/POS"]
    pos = pos[vcf["variants/CHROM"] == chrom]
    pi = allel.sequence_diversity(pos, ac, start=start, stop=stop)
    return(pi)




