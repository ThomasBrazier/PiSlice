#0 -*- coding: utf-8 -*-
"""
Estimate statistics from sequence alignments
"""

import egglib

def create_align(fasta, vcf, chromosome, start, end, ploidy=2):
    """
    Create multi-samples alignment from a reference fasta and a vcf
    :fasta: str, a fasta reference file
    :vcf: a sckit-allel vcf format
    :chromosome: string, the chromosome name
    :start: int, start position, +1 index
    :stop: int, stop position, +1 index
    :return: list, a list of tuples of type [("sample_name","sequence")]
    """
    fas = fasta.sample_sequence(chromosome, start, end)
    samples = vcf["samples"]
    pos = vcf["variants/POS"]
    gt = vcf["calldata/GT"]
    gt = gt[(pos >= start) & (pos <= end)]
    ref = vcf["variants/REF"]
    ref = ref[(pos >= start) & (pos <= end)]
    alt = vcf["variants/ALT"]
    alt = alt[(pos >= start) & (pos <= end)]
    # Index of nucleotides to change, +1 base
    idx = pos[(pos >= start) & (pos <= end)]
    # shifted to match sequence
    idx = idx - start + 1
    align = list()
    for j, sam in enumerate(samples):
        # Take care of diploids - ploidy level
        for p in range(1, (ploidy + 1)):
            # reference sequence
            seq = fas
            # change nucleotides
            for i, x in enumerate(idx):
                allel = gt[:,samples == samples[j]]
                allel = allel[i][0][(p - 1)]
                if allel == '0':
                    base = ref[i]
                elif (allel == '1') or (allel == '2') or (allel == '3'):
                    base = alt[i][int(allel) - 1]
                else:
                    base = "N"
                if base == "":
                    base = "N"
                temp = list(seq)
                temp[int(x) - 1] = str(base)
                seq = "".join(temp)
            sample = sam + "-" + str(p)
            align.append((sample, seq))
    return(align)


def codon_align(fasta, vcf, gff, chromosome, start, end, ploidy=2):
    """
    Create multi-samples alignment of coding seauences
    :fasta: str, a fasta reference file
    :vcf: a sckit-allel vcf format
    :gff: a gff object with CDS information
    :chromosome: string, the chromosome name
    :start: int, start position, +1 index
    :stop: int, stop position, +1 index
    :return: list, a list of tuples of type [("sample_name","full coding sequence")]
    """
    # get CDS coordinates in the region

    # get sequences for each interval (CDS part)

    # clean up start and stop codons

    # check size (multiple of 3 because codons take three nucleotides)

    # concatenate CDS sequences (exons) to get the full sequence
