import pandas as pd

from PiSlice.input import fasta

def test_open_fasta():
    fasta_file = "fasta.fa.gz"
    assert (fasta(fasta_file))

def test_fasta_length():
    fasta_file = "fasta.fa.gz"
    genome = fasta(fasta_file)
    assert (genome.length() == [980, 2170])

def test_seqname():
    fasta_file = "fasta.fa.gz"
    genome = fasta(fasta_file)
    assert (genome.seqname() == ['CP002684.1', 'CP002685.1'])

def test_sequence():
    fasta_file = "fasta.fa.gz"
    genome = fasta(fasta_file)
    assert (genome.sequence() == ['CCCTAAACCCTAAACCCTAAACCCTAAACCTCTGAATCCTTAATCCCTAAATCCCTAAATCTTTAAATCCTACATCCATGAATCCCTAAATACCTAATTCCCTAAACCCGAAACCGGTTTCTCTGGTTGAAAATCATTGTGTATATAATGATAATTTTATCGTTTTTATGTAATTGCTTATTGTTGTGTGTAGATTTTTTAAAAATATCATTTGAGGTCAATACAAATCCTATTTCTTGTGGTTTTCTTTCCTTCACTTAGCTATGGATGGTTTATCTTCATTTGTTATATTGGATACAAGCTTTGCTACGATCTACATTTGGGAATGTGAGTCTCTTATTGTAACCTTAGGGTTGGTTTATCTCAAGAATCTTATTAATTGTTTGGACTGTTTATGTTTGGACATTTATTGTCATTCTTACTCCTTTGTGGAAATGTTTGTTCTATCAATTTATCTTTTGTGGGAAAATTATTTAGTTGTAGGGATGAAGTCTTTCTTCGTTGTTGTTACGCTTGTCATCTCATCTCTCAATGATATGGGATGGTCCTTTAGCATTTATTCTGAAGTTCTTCTGCTTGATGATTTTATCCTTAGCCAAAAGGATTGGTGGTTTGAAGACACATCATATCAAAAAAGCTATCGCCTCGACGATGCTCTATTTCTATCCTTGTAGCACACATTTTGGCACTCAAAAAAGTATTTTTAGATGTTTGTTTTGCTTCTTTGAAGTAGTTTCTCTTTGCAAAATTCCTCTTTTTTTAGAGTGATTTGGATGATTCAAGACTTCTCGGTACTGCAAAGTTCTTCCGCCTGATTAATTATCCATTTTACCTTTGTCGTAGATATTAGGTAATCTGTAAGTCAACTCATATACAACTCATAATTTAAAATAAAATTATGATCGACACACGTTTACACATAAAATCTGTAAATCAACTCATATACCCGTTATTCCCACAATCATATGCTTTCTAAAAGC', 'NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNGAATTCGTCGACCAGGACGGCGGAATGCTCGACCAGGACGATGAATGGGCGATGAAAATCTATCGGGTTAGAGGAATGGTCGACCGGGTCCGAGAATTCGTCGACCAGGACGAGGAGTGGTCGAGGATTTGTCGACCAGGAGTTGAAATCGTCGACCGGGTCCGAGAATTCGTCGACCAGGACGGCGGAACCCTCGACCAGGACGATGAATGGGCGATGAAAATCTATCGGGTTCGAGGAATGGTCGACCAGAAGTTGGAATCGTCGACCGGGTCCGAGAATTCGTCGACCAGGCCGAGGAGTGGTCGAGGATTTGTCGACCAGGAGTTGAAATCGTCGACCGGGACCGAGAATTCGTCGACCAGGACGGCGGAACCCTCGACCAGGACGATGAATGGGCGATGAAAATCTATCGGGTTCGAGGAATGGTCGACCAGGGGTTGAAATCGTCGACCAGGTCCGAGACTTCATCGACCGGGTCCGAGGATTCGTCGACCAGGACGGCCGGATGTCCGAGAAAAAAAAATGTTGCCGAATAACTTTCGAAAATCATTGGATATGATGCAATGTTTTGTGATCGAATCTCTTAAAATACATCAATAAAGAGTTTAGGATGTCAAGTTTGCATCAAATATGCCCACGGAGCCCCAACTAGACCATGAAAATCCGTTGTTGTATCAGGTCAAATGACCTAGCTAGAGGTGTCAGAAAATTATGTAAATTTACCAGAAAATAGGATTTAGTATCCTTATGATGCATGCTAAAAAGAATTTTCAAATTCCAAGTATTTCTTTTTTTTTGGCACCGGTGTCTCCTCAGACATTTCAATGTCTGTTGGTGCCAAGAGGGAAAAGGGCTATTAAGCTATATAGGGGGGTGGGTGTTGAGGGAGTCTGGGCAGTCCGTGGGAACCCCCTTTTTCGGTTCGGACTTGGGTAGCGATCGAGGGATGGTATCGGATATCGACACGAGGAATGACCGACCGTCCGGCCGCCGGGATTTTCGCCGGAAAACTTTTTCGGCGACTTTTCCGGCGATCGGTTTTGTTGCCTTTTTCCGAGTTTTCTCAGCAGTTCTCGGACAAAAACTGATGAATCGTCGAGGAGAATGAGCTTGCCTTGCGTGGGCTGCCATTAGTTCTTCGAGGCGTTAGGGTGGCGGCGGTATAAA'])

def test_sample_chromosome():
    fasta_file = "fasta.fa.gz"
    genome = fasta(fasta_file)
    assert (genome.sample_chromosome('CP002684.1') == 'CCCTAAACCCTAAACCCTAAACCCTAAACCTCTGAATCCTTAATCCCTAAATCCCTAAATCTTTAAATCCTACATCCATGAATCCCTAAATACCTAATTCCCTAAACCCGAAACCGGTTTCTCTGGTTGAAAATCATTGTGTATATAATGATAATTTTATCGTTTTTATGTAATTGCTTATTGTTGTGTGTAGATTTTTTAAAAATATCATTTGAGGTCAATACAAATCCTATTTCTTGTGGTTTTCTTTCCTTCACTTAGCTATGGATGGTTTATCTTCATTTGTTATATTGGATACAAGCTTTGCTACGATCTACATTTGGGAATGTGAGTCTCTTATTGTAACCTTAGGGTTGGTTTATCTCAAGAATCTTATTAATTGTTTGGACTGTTTATGTTTGGACATTTATTGTCATTCTTACTCCTTTGTGGAAATGTTTGTTCTATCAATTTATCTTTTGTGGGAAAATTATTTAGTTGTAGGGATGAAGTCTTTCTTCGTTGTTGTTACGCTTGTCATCTCATCTCTCAATGATATGGGATGGTCCTTTAGCATTTATTCTGAAGTTCTTCTGCTTGATGATTTTATCCTTAGCCAAAAGGATTGGTGGTTTGAAGACACATCATATCAAAAAAGCTATCGCCTCGACGATGCTCTATTTCTATCCTTGTAGCACACATTTTGGCACTCAAAAAAGTATTTTTAGATGTTTGTTTTGCTTCTTTGAAGTAGTTTCTCTTTGCAAAATTCCTCTTTTTTTAGAGTGATTTGGATGATTCAAGACTTCTCGGTACTGCAAAGTTCTTCCGCCTGATTAATTATCCATTTTACCTTTGTCGTAGATATTAGGTAATCTGTAAGTCAACTCATATACAACTCATAATTTAAAATAAAATTATGATCGACACACGTTTACACATAAAATCTGTAAATCAACTCATATACCCGTTATTCCCACAATCATATGCTTTCTAAAAGC')


def test_sample_sequence():
    fasta_file = "fasta.fa.gz"
    genome = fasta(fasta_file)
    assert (genome.sample_sequence('CP002684.1',5,14) == 'AAACCCTAAA')

def test_sample_sequence_masked():
    fasta_file = "fasta.fa.gz"
    genome = fasta(fasta_file)
    mask = [(1, 4), (10, 12)]
    assert (genome.sample_sequence_masked('CP002684.1', 1, 14, mask) == ['AAACC', 'AA'])

def test_read_gff():
    from PiSlice.input import read_gff
    import pandas as pd
    gff_file = "GFF_test.gff.gz"
    gff = read_gff(gff_file, parse=False, parse_introns=False, parse_utr=False, infer_rank=False, n_cpus=8)
    assert (gff.shape == (85, 14)) & (isinstance(gff, pd.DataFrame))

def test_region():
    from PiSlice.input import read_gff
    gff_file = "GFF_test.gff.gz"
    gff = read_gff(gff_file, parse=False, parse_introns=False, parse_utr=False, infer_rank=False, n_cpus=8)
    gff = gff.gff.region(1, 5000, "NC_003070.9")
    assert (gff["start"].min() >= 1) & (gff["end"].max() <= 5000) & (gff["seqname"].unique() == "NC_003070.9")


def test_id():
    from PiSlice.input import read_gff
    gff_file = "GFF_test.gff.gz"
    gff = read_gff(gff_file, parse=True, parse_introns=False, parse_utr=False, infer_rank=False, n_cpus=8)
    gff = gff.gff.id("gene-AT1G01010")
    assert (gff["id"].unique() == "gene-AT1G01010")


def test_parent():
    from PiSlice.input import read_gff
    gff_file = "GFF_test.gff.gz"
    gff = read_gff(gff_file, parse=True, parse_introns=False, parse_utr=False, infer_rank=False, n_cpus=8)
    single_parent = gff.gff.parent("exon-NM_099983.2-1")
    assert (isinstance(single_parent, pd.DataFrame)) & (single_parent["id"].unique() == "rna-NM_099983.2")

def test_children():
    from PiSlice.input import read_gff
    gff_file = "GFF_test.gff.gz"
    gff = read_gff(gff_file, parse=True, parse_introns=False, parse_utr=False, infer_rank=False, n_cpus=8)
    single_child = gff.gff.children("gene-AT1G01010", all=False)
    recursive_child = gff.gff.children("gene-AT1G01010", all=True)
    assert (isinstance(single_child, pd.DataFrame)) &\
           (single_child["id"].unique() == "rna-NM_099983.2") &\
           (isinstance(recursive_child, pd.DataFrame)) &\
           (recursive_child.shape[0] == 13)


def test_feature():
    from PiSlice.input import read_gff
    gff_file = "GFF_test.gff.gz"
    gff = read_gff(gff_file, parse=True, parse_introns=False, parse_utr=False, infer_rank=False, n_cpus=8)
    mrna = gff.gff.feature("mRNA")
    exon = gff.gff.feature("exon")
    assert (isinstance(mrna, pd.DataFrame)) & (mrna["feature"].unique() == "mRNA") & \
           (isinstance(exon, pd.DataFrame)) & (exon["feature"].unique() == "exon")


def test_rank():
    from PiSlice.input import read_gff
    gff_file = "GFF_test.gff.gz"
    gff = read_gff(gff_file, parse=True, parse_introns=False, parse_utr=False, infer_rank=True, n_cpus=8)
    r1 = gff.gff.rank(1)
    assert (isinstance(r1, pd.DataFrame)) & (r1["rank"].unique() == 1)