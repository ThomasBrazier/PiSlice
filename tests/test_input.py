from PiSlice.input import fasta

def test_fasta():
    assert (cpg("ATTANNCGCGGCCCGATA") == 0.375)
