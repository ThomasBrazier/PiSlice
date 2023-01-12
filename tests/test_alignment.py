def test_create_align():
    import PiSlice.alignment as align
    import PiSlice.input as input
    import numpy as np
    # Simulate a vcf with three diploid individuals and four SNPs
    vcf = dict([
        ("samples", np.array(['A','B','C'])),
        ("calldata/GT", np.array([[[0,0], [0,1], ['.',1]],
                                   [[0,0], [0,0], [0,0]],
                                   [[0,0], [0,1], [1,1]],
                                   [[0,0], [0,1], [2,1]]])),
        ("variants/ALT", np.array([['A', '', ''], ['G', '', ''], ['T', '', ''], ['C', '', '']])),
        ("variants/CHROM", np.array(['1', '1', '1', '1'])),
        ("variants/FILTER_PASS", np.array([True, True, True, True])),
        ("variants/ID", np.array(['.', '.', '.', '.'])),
        ("variants/POS", np.array([1, 4, 10, 15])),
        ("variants/QUAL", np.array([40, 40, 40, 40])),
        ("variants/REF", np.array(['T', 'A', 'A', 'G']))
    ])
    fasta_file = "test.fna.gz"
    fasta = input.fasta(fasta_file)
    expected = [('A-1', 'TCCAAAACCATAAAGCCTAA'),
                ('A-2', 'TCCAAAACCATAAAGCCTAA'),
                ('B-1', 'TCCAAAACCATAAAGCCTAA'),
                ('B-2', 'ACCAAAACCTTAAACCCTAA'),
                ('C-1', 'NCCAAAACCTTAAANCCTAA'),
                ('C-2', 'ACCAAAACCTTAAACCCTAA')]

    chromosome = "1"
    start = 1
    end = 20
    aln = align.create_align(fasta, vcf, chromosome, start, end, ploidy=2)
    aln

    assert(aln == expected)
