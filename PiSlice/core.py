#0 -*- coding: utf-8 -*-
"""
The PiSlice core function
"""

import mapply
import multiprocessing
import PiSlice.diversity as div
import PiSlice.nucleotide as nuc
import PiSlice.alignment as align


def piSlice(windows, statistics=[""], min_bp=6, splicing_strategy="merge", n_cpus=6, *args, **kwargs):
    """
    The main function to return a data frame of population genomics statistics for a list of genomic windows.
    :param windows: DataFrame, a pandas data frame (can be gff) with at least three columns: seqname, start, end
    :param statistics: str, a list of statistics to compute
    :param **fasta: fasta, a fasta object with multiple fasta sequences
    :param **gff: DataFrame, a gff object
    :param **vcf: vcf, a vcf object
    :return: DataFrame, a data frame with population statistics for each windows
    """
    fasta = kwargs.get("fasta", "")
    gff = kwargs.get("gff", "")
    vcf = kwargs.get("vcf", "")
    #pandarallel.initialize(nb_workers=n_cpus, progress_bar=True)

    # Function to subset sequences in the fasta file
    def make_dataset(windows, fasta):
        # Sample sequences
        # Sample all sequences from chromosomes and start-end positions
        list_seq = list(windows.apply(lambda x: fasta.sample_sequence(x["seqname"], x["start"], x["end"]), axis=1))
        return(list_seq)

    # TODO A progress bar
    # Header
    print("Number of windows:", len(windows.index))
    print("Chromosomes are", " ".join(windows.seqname.unique()))
    if (n_cpus == 0):
        n_cpus = multiprocessing.cpu_count()
    sensible_cpus = mapply.parallel.sensible_cpu_count()
    mapply.init(n_workers=min(sensible_cpus, n_cpus))

    if "gene_count" in statistics:
        print("Process number of genes")
        estimates = windows.mapply(lambda x: nuc.gene_count(gff,
                                             x["seqname"],
                                             x["start"],
                                             x["end"]),
                             axis=1)
        windows["gene_count"] = estimates

    if "gene_length" in statistics:
        print("Process mean gene length (bp)")
        estimates = windows.mapply(lambda x: nuc.feature_length(gff,
                                             x["seqname"],
                                             x["start"],
                                             x["end"],
                                             feature="gene"),
                             axis=1)
        windows["gene_length"] = estimates

    if "exon_length" in statistics:
        print("Process mean exon length (bp)")
        estimates = windows.mapply(lambda x: nuc.feature_length(gff,
                                             x["seqname"],
                                             x["start"],
                                             x["end"],
                                             feature="exon"),
                             axis=1)
        windows["exon_length"] = estimates

    if "intron_length" in statistics:
        print("Process mean intron length (bp)")
        estimates = windows.mapply(lambda x: nuc.feature_length(gff,
                                             x["seqname"],
                                             x["start"],
                                             x["end"],
                                             feature="intron"),
                             axis=1)
        windows["intron_length"] = estimates

    if "gene_nbexons" in statistics:
        print("Process the mean number of exons")
        estimates = windows.mapply(lambda x: nuc.gene_nbexons(gff,
                                             x["seqname"],
                                             x["start"],
                                             x["end"]),
                             axis=1)
        windows["gene_nbexons"] = estimates

    if "gene_density" in statistics:
        print("Process gene density")
        estimates = windows.mapply(lambda x: nuc.gene_density(gff,
                                             x["seqname"],
                                             x["start"],
                                             x["end"]),
                             axis=1)
        windows["gene_density"] = estimates

    if "snp_count" in statistics:
        print("Process number of SNPs")
        estimates = windows.mapply(lambda x: div.snp_count(vcf,
                                             x["seqname"],
                                             x["start"],
                                             x["end"]),
                             axis=1)
        windows["snp_count"] = estimates

    if "snp_density" in statistics:
        print("Process SNP density")
        estimates = windows.mapply(lambda x: div.snp_density(vcf,
                                             x["seqname"],
                                             x["start"],
                                             x["end"]),
                             axis=1)
        windows["snp_density"] = estimates

    if "pi" in statistics:
        print("Process nucleotide diveristy (Pi)")
        estimates = windows.mapply(lambda x: div.pi(vcf,
                                             x["seqname"],
                                             x["start"],
                                             x["end"]),
                             axis=1)
        windows["pi"] = estimates

        if "theta_watterson" in statistics:
            print("Process Theta Watterson")
            estimates = windows.mapply(lambda x: div.theta_watterson(vcf,
                                                        x["seqname"],
                                                        x["start"],
                                                        x["end"]),
                                       axis=1)
            windows["theta_watterson"] = estimates

    if "tajima_d" in statistics:
        print("Process Tajima's D")
        estimates = windows.mapply(lambda x: div.tajima_d(vcf,
                                                                 x["seqname"],
                                                                 x["start"],
                                                                 x["end"]),
                                   axis=1)
        windows["tajima_d"] = estimates

    if "gc" in statistics:
        print("Process GC content")
        list_seq = make_dataset(windows, fasta)
        # Compute GC content
        estimates = list(map(lambda x: nuc.gc(x), list_seq))
        # Add column for statistics
        windows["gc"] = estimates

    if "gc_noncoding" in statistics:
        print("Process non-coding GC content")
        estimates = windows.mapply(lambda x: nuc.gc_noncoding(fasta,
                                             gff,
                                             x["seqname"],
                                             x["start"],
                                             x["end"],
                                             min_bp=min_bp),
                             axis=1)
        list_gc = [item[0] for item in estimates]
        list_density = [item[1] for item in estimates]
        # Add column for statistics
        windows["gc_noncoding"] = list_gc
        windows["noncoding_proportion"] = list_density

    if "gc_intergenic" in statistics:
        print("Process intergenic GC content")
        estimates = windows.mapply(lambda x: nuc.gc_intergenic(fasta,
                                             gff,
                                             x["seqname"],
                                             x["start"],
                                             x["end"],
                                             min_bp=min_bp),
                             axis=1)
        list_gc = [item[0] for item in estimates]
        list_density = [item[1] for item in estimates]
        # Add column for statistics
        windows["gc_intergenic"] = list_gc
        windows["intergenic_proportion"] = list_density

    if "gc_intron" in statistics:
        print("Process intron GC content")
        estimates = windows.mapply(lambda x: nuc.gc_intron(fasta,
                                             gff,
                                             x["seqname"],
                                             x["start"],
                                             x["end"],
                                             min_bp=min_bp,
                                             splicing_strategy=splicing_strategy),
                             axis=1)
        list_gc = [item[0] for item in estimates]
        list_density = [item[1] for item in estimates]
        # Add column for statistics
        windows["gc_intron"] = list_gc
        windows["intron_proportion"] = list_density


    if "gc_codon" in statistics:
        print("Process GC content with codon positions")
        # Compute GC content
        estimates = windows.mapply(lambda x: nuc.gc_codon(fasta,
                                             gff,
                                             x["seqname"],
                                             x["start"],
                                             x["end"],
                                             min_bp=min_bp),
                             axis=1)
        list_gc = [item[0] for item in estimates]
        list_gc1 = [item[1] for item in estimates]
        list_gc2 = [item[2] for item in estimates]
        list_gc3 = [item[3] for item in estimates]
        list_cds_proportion = [item[4] for item in estimates]
        # Add column for statistics
        windows["gc_codon"] = list_gc
        windows["gc1"] = list_gc1
        windows["gc2"] = list_gc2
        windows["gc3"] = list_gc3
        windows["cds_proportion"] = list_cds_proportion

    if "gc3exon1" in statistics:
        print("Process GC3 first exon")
        estimates = windows.mapply(lambda x: nuc.gc3exon1(fasta,
                                                     gff,
                                                     x["seqname"],
                                                     x["start"],
                                                     x["end"],
                                                     min_bp=min_bp),
                                  axis=1)
        windows["gc3_exon1"] = estimates

    if "cpg" in statistics:
        print("Process CpG densities")
        list_seq = make_dataset(windows, fasta)
        # Compute CpG density
        estimates = list(map(lambda x: nuc.cpg(x), list_seq))
        # Add column for statistics
        windows["cpg"] = estimates

    if "seq" in statistics:
        print("Estimating nucleotide diversity (Pi)")
        sequences = list(map(lambda x: fasta.sample_sequence(windows.loc[x, "seqname"],
                                              windows.loc[x, "start"],
                                              windows.loc[x, "end"]),
                             windows.index))
        windows["seq"] = sequences

    if "pi" in statistics:
        print("Estimating nucleotide diversity (Pi) from a VCF")
        stat = list(map(lambda x: div.pi(vcf, windows.loc[x, "seqname"],
                                            windows.loc[x, "start"],
                                            windows.loc[x, "end"]),
                             windows.index))
        windows["Pi"] = stat

    if "theta_watterson" in statistics:
        print("Estimating Theta Watterson from a VCF")
        stat = list(map(lambda x: div.theta_watterson(vcf, windows.loc[x, "seqname"],
                                            windows.loc[x, "start"],
                                            windows.loc[x, "end"]),
                             windows.index))
        windows["Theta_Watterson"] = stat

    if "tajima_d" in statistics:
        print("Estimating Tajima's D from a VCF")
        stat = list(map(lambda x: div.tajima_d(vcf, windows.loc[x, "seqname"],
                                            windows.loc[x, "start"],
                                            windows.loc[x, "end"]),
                             windows.index))
        windows["Tajima_D"] = stat

    if "pi_alignment" in statistics:
        print("Estimating nucleotide diversity (Pi) from sequence alignment")
        stats = list(map(lambda x: align.pi_alignment(fasta, vcf,
                                                          windows.loc[x, "seqname"],
                                                          windows.loc[x, "start"],
                                                          windows.loc[x, "end"],
                                                          ploidy=2,
                                                          max_mising=0.05),
                             windows.index))
        windows["Pi.alignment"] = stats["Pi"]
        windows["S.alignment"] = stats["S"]
        windows["lseff.alignment"] = stats["lseff"]
        windows["nseff.alignment"] = stats["nseff"]

    if "pi_coding" in statistics:
        print("Estimating nucleotide diversity (Pi/PiN/PiS) in coding sequences")
        stats = list(map(lambda x: align.pi_coding(fasta, vcf,
                                                          windows.loc[x, "seqname"],
                                                          windows.loc[x, "start"],
                                                          windows.loc[x, "end"],
                                                          ploidy=2,
                                                          max_mising=0.05),
                             windows.index))
        windows["Pi.coding"] = stats["Pi"]
        windows["S.coding"] = stats["S"]
        windows["lseff.coding"] = stats["lseff"]
        windows["nseff.coding"] = stats["nseff"]
        windows["Pi4"] = stats["Pi4"]
        windows["Pi0"] = stats["Pi0"]
        windows["nS"] = stats["nS"]
        windows["nNS"] = stats["nNS"]
        windows["npolS"] = stats["npolS"]
        windows["npolNS"] = stats["npolNS"]


    return windows
