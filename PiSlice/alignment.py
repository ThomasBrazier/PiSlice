#0 -*- coding: utf-8 -*-
"""
Estimate statistics from sequence alignments
"""

import egglib
import warnings
from itertools import chain
import pandas as pd

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
                allel = str(allel)
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
            if len(seq) > 0 and isinstance(seq, str):
                align.append((sample, seq))
    return(align)


def transcript_align(fasta, vcf, gff, geneid, ploidy=2):
    """
    Create multi-samples alignment of a single transcript
    Requires a gene/mRNA GFF features to get a transcript
    :fasta: str, a fasta reference file
    :vcf: a sckit-allel vcf format
    :gff: a gff object with CDS information
    :geneid: string, a gene ID according to the GFF
    :return: list, a list of tuples of type [("sample_name","full coding sequence")]
    """
    # keep a single transcript - the first one
    mrna = gff.gff.children(geneid)
    if len(mrna.gff.feature("mRNA")) > 0:
        m = mrna.gff.feature("mRNA")
        m = m.iloc[0]
        m = m["id"]
        id = m
    else:
        id = geneid

    # get mRNA/CDS coordinates in the region
    pos = gff.gff.children(id)
    pos = pos.gff.feature("CDS")

    # get sequences for each interval (CDS part)
    # a list of align objects
    # Iterating over multiple columns - differing data type
    cdsparts = [create_align(fasta, vcf, row[0], row[1], row[2], ploidy=ploidy) for row in
                zip(pos["seqname"], pos["start"], pos["end"])]

    if len(cdsparts) > 0:
        # concatenate CDS sequences (exons) to get the full transcript
        # tmp = [[(1, 'aaaa'), ('b', 'ggg'), ('c', 'agtc')], [(1, 'ccc'), ('b', 'ttt'), ('c', 'agct')]]
        tmp = cdsparts
        nsamples = len(tmp[0])
        concat = list()
        for i in range(0, nsamples):
            concat.append(([n[0] for n in [x[i] for x in tmp]][0], ''.join([x[1] if isinstance(x[1], str) else "" for x in [aln[i] for aln in tmp]])))

        # Take care of strand
        # Reverse '-' strand
        strand = pos["strand"].iloc[0]
        for idx, aln in enumerate(cdsparts):
            if strand == "-":
                aln = [(sample, seq[::-1]) for sample, seq in aln]
                cdsparts[idx] = aln

        # Frame shift (phase feature)
        # While taking all CDS parts and concatenating them, no need to do a phase shift
        # See https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md for details on phase
        # frame = pos["frame"]
        # for idx, aln in enumerate(cdsparts):
        #     aln = [(sample, seq[int(frame.iloc[idx]):]) for sample, seq in aln]
        #     cdsparts[idx] = aln

        # clean up start and stop codons
        # for each CDS part get 5' and 3' codons
        # remove them if 5' == 'ATG'
        # or 3' == ['TAA', 'TAG', 'TGA']
        for idx, aln in enumerate(cdsparts):
            aln = [(sample, seq[3:]) if seq[:3] == "ATG" else (sample, seq) for sample, seq in aln]
            aln = [(sample, seq[:-3]) if (seq[-3:] == "TAA" or seq[-3:] == "TAG" or seq[-3:] == "TGA") else (sample, seq)
                   for sample, seq in aln]
            cdsparts[idx] = aln

        # Check the full CDS sequence is a multiple of 3
        checksize = list()
        goodsize = [True if len(seq) % 3 == 0 else False for sample, seq in concat]
        checksize.append(goodsize)
        if any(False in item for item in checksize):
            warnings.warn("The CDS is not a multiple of 3 in gene:" + geneid + "!")
            cdsparts = []
    else:
        cdsparts = []

    return(cdsparts)


def codon_align(fasta, vcf, gff,  chromosome, start, end, ploidy=2):
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
    # Process alignment for a single exon or any window shorter than gene.
    #  Detect with the gff if start and end are those of a single exon/CDS.
    #  Make frame shift and reduce to multiple of three
    # Detect from start:end coordinates if the sequence is a single exon/CDS
    window = gff.gff.region(start, end, chromosome)
    window = window.gff.feature("CDS")
    if window.shape[0] == 1:
        mrna = create_align(fasta, vcf, chromosome, start, end, ploidy=ploidy)
        phase = window["frame"]
        mrna = [(i, x[int(phase):]) for i, x in mrna]
        # keep a length multiple of 3
        mrna = [(i, x[0:(len(x) - (len(x) % 3))]) for i, x in mrna]
        concat = mrna
    else:
        # keep only complete CDS/transcripts
        gene = gff.gff.feature("gene")
        gene = gene.gff.region(start, end, chromosome)
        if isinstance(gene["id"], str):
            geneid = gene["id"]
        elif isinstance(gene["id"], pd.Series):
            geneid = list(gene["id"])
        else:
            geneid = list()

        mrna = list()
        [mrna.append(transcript_align(fasta, vcf, gff, str(i), ploidy=ploidy)) for i in geneid]
        mrna = list(chain.from_iterable(mrna))

        # concatenate transcripts to get the full sequence over the region
        # tmp = [[('a', 'aaaa'), ('b', 'ggg'), ('c', 'agtc')], [('a', 'ccc'), ('b', 'ttt'), ('c', 'agct')]]
        tmp = mrna
        nsamples = len(tmp[0])
        concat = list()
        for i in range(0, nsamples):
            concat.append(([n[0] for n in [x[i] for x in tmp]][0], ''.join([x[1] for x in [aln[i] for aln in tmp]])))

    # Check the full CDS sequence is a multiple of 3
    checksize = list()
    goodsize = [True if len(seq) % 3 == 0 else False for sample, seq in concat]
    checksize.append(goodsize)
    if any(False in item for item in checksize):
        warnings.warn("The CDS is not a multiple of 3 in gene:", geneid, "!")

    return(concat)


def pi_alignment(fasta, vcf, chromosome, start, end, ploidy=2, max_missing=0.05):
    """
    Estimate Pi over a given genomic region.
    Beware that missing data is inferred ONLY for polymorphic sites from the vcf,
    thus statistics on a genomic region can be underestimated due to the low number of missing data
    For a discussion on this issue see https://onlinelibrary.wiley.com/doi/10.1111/1755-0998.13326
    S is the number of polymorphic sites
    Pi is the nucleotide diversity
    lseff is the number of sites used for analysis (excluding those with either too many missing data or too many alleles)
    nseff is the average number of used samples among included sites
    Divide by lseff to get the value of Pi per site
    :fasta: str, a fasta reference file
    :vcf: a sckit-allel vcf format
    :chromosome: string, the chromosome name
    :start: int, start position, +1 index
    :stop: int, stop position, +1 index
    :ploidy: int, ploidy level
    :max_missing: float, max proportion of missing data per site
    :return: dict, S, Pi, lseff, nseff
    """
    # get the sequence of the region
    sequence = create_align(fasta, vcf, chromosome, start, end, ploidy=ploidy)

    # convert to Egglib object
    aln = egglib.Align.create(sequence, alphabet=egglib.alphabets.DNA)

    # compute stats
    cs = egglib.stats.ComputeStats()
    cs.add_stats('Pi', 'S', 'lseff', 'nseff')
    stats = cs.process_align(aln, max_missing=max_missing)
    stats

    return(stats)



def pi_coding(fasta, vcf, gff, chromosome, start, end, ploidy=2, max_missing=0.05):
    """
    Estimate PiN and PiS (nucleotide diversity in non-synonymous and synonymous positions)
    over a given coding region.
    Beware that missing data is inferred ONLY for polymorphic sites from the vcf,
    thus statistics on a genomic region can be underestimated due to the low number of missing data
    For a discussion on this issue see https://onlinelibrary.wiley.com/doi/10.1111/1755-0998.13326
    S is the number of polymorphic sites
    Pi is the nucleotide diversity in the coding region (total, divide by lseff to get the per site value)
    Pi0 is Pi for non-synonymous polymorphic positions only, 0-fold degenerated (per site, divided by the number of NS sites)
    Pi4 is Pi for synonymous polymorphic positions only, 4-fold degenerated (per site, divided by the number of S sites)
    lseff is the number of sites used for analysis (excluding those with either too many missing data or too many alleles)
    nseff is the average number of used samples among included sites
    nNS is the count of non-synonymous positions
    nS is the count of synonymous positions
    npolNS is the count of polymorphic non-synonymous positions
    npolS is the count of polymorphic synonymous positions
    :fasta: str, a fasta reference file
    :vcf: a sckit-allel vcf format
    :gff: a gff object with CDS information
    :chromosome: string, the chromosome name
    :start: int, start position, +1 index
    :stop: int, stop position, +1 index
    :ploidy: int, ploidy level
    :max_missing: float, max proportion of missing data per site
    :return: dict, S, Pi, Pi4, Pi0, lseff, nseff, nS, nNS, npolS, npolNS
    """
    # get the CDS of the region
    cds = codon_align(fasta, vcf, gff, chromosome, start, end, ploidy=ploidy)

    # convert to Egglib object
    aln = egglib.Align.create(cds, alphabet=egglib.alphabets.DNA)

    # compute stats
    cs = egglib.stats.ComputeStats()
    cs.add_stats('Pi', 'S', 'lseff', 'nseff')
    stats_region = cs.process_align(aln, max_missing=max_missing)
    stats_region

    # get synonymous/non-synonymous positions
    codons = aln
    codons.to_codons()
    synnonsyn = egglib.stats.CodingDiversity()
    synnonsyn.process(codons, code=1, max_missing=max_missing)

    # count number of sites of each class
    npolS = synnonsyn.num_pol_S
    npolNS = synnonsyn.num_pol_NS
    nS = synnonsyn.num_sites_S
    nNS = synnonsyn.num_sites_NS

    # estimate PiN/PiS per polymorphic site
    pos_S = [int(i) for i in synnonsyn.positions_S]
    pos_NS = [int(i) for i in synnonsyn.positions_NS]

    sub_S = aln.extract(pos_S)
    cs = egglib.stats.ComputeStats()
    cs.add_stats('Pi')
    stats_S = cs.process_align(sub_S, max_missing=max_missing)
    sub_NS = aln.extract(pos_NS)
    cs = egglib.stats.ComputeStats()
    cs.add_stats('Pi')
    stats_NS = cs.process_align(sub_NS, max_missing=max_missing)

    # create output
    stats_synnonsyn = {
        "Pi4": stats_S["Pi"],
        "Pi0": stats_NS["Pi"],
        "nS": nS,
        "nNS": nNS,
        "npolS": npolS,
        "npolNS": npolNS
    }
    stats = dict(stats_region, **stats_synnonsyn)

    return(stats)