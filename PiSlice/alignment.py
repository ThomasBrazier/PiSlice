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
    Create multi-samples alignment of coding sequences
    Concatenate the coding sequences over the given region
    Use it on gene/mRNA GFF features to get a transcript
    Use it on exons GFF features to get individual spliced CDS parts
    :fasta: str, a fasta reference file
    :vcf: a sckit-allel vcf format
    :gff: a gff object with CDS information
    :chromosome: string, the chromosome name
    :start: int, start position, +1 index
    :stop: int, stop position, +1 index
    :return: list, a list of tuples of type [("sample_name","full coding sequence")]
    """
    # get CDS coordinates in the region
    pos = gff.gff.feature("CDS")
    pos = pos.gff.region(start, end, chromosome)

    # get sequences for each interval (CDS part)
    # a list of align objects
    # Iterating over multiple columns - differing data type
    cdsparts = [align.create_align(fasta, vcf, row[0], row[1], row[2], ploidy=2) for row in zip(pos["seqname"], pos["start"], pos["end"]) if len(row) > 0]

    # Take care of strand
    # Reverse '-' strand
    strand = pos["strand"]
    for idx, aln in enumerate(cdsparts):
        if strand.iloc[idx] == "-":
            aln = [(sample, seq[::-1]) for sample, seq in aln]
            cdsparts[idx] = aln

    # Remove None values
    # cdsparts = [c for c in cdsparts if c is not None]

    # Frame shift (phase feature)
    frame = pos["frame"]
    for idx, aln in enumerate(cdsparts):
        aln = [(sample, seq[int(frame.iloc[idx]):]) for sample, seq in aln]
        cdsparts[idx] = aln

    # clean up start and stop codons
    # for each CDS part get 5' and 3' codons
    # remove them if 5' == 'ATG'
    # or 3' == ['TAA', 'TAG', 'TGA']
    for idx, aln in enumerate(cdsparts):
        aln = [(sample, seq[3:]) if seq[:3] == "ATG" else (sample, seq) for sample, seq in aln]
        aln = [(sample, seq[:-3]) if (seq[-3:] == "TAA" or seq[-3:] == "TAG" or seq[-3:] == "TGA") else (sample, seq) for sample, seq in aln]
        cdsparts[idx] = aln

    # check size of the CDS part (multiple of 3 because codons take three nucleotides)
    # Remove CDS parts with wrong size
    # for idx, aln in enumerate(cdsparts):
    #     aln = [(sample, seq) for sample, seq in aln if len(seq) % 3 == 0]
    #     cdsparts[idx] = aln

    # concatenate CDS sequences (exons) to get the full sequence over the region
    tmp = [[('a', 'aaaa'), ('b', 'ggg')], [('a', 'ccc'), ('b', 'ttt')]]
    nsamples = len(tmp)
    concat = list()
    for i in range(0, nsamples):
        concat.append(([n[0] for n in [x[i] for x in tmp]][0], ''.join([x[1] for x in [aln[i] for aln in tmp]])))

    # Check the full CDS sequence is a multiple of 3

    return(cds)
