from PiSlice.popstatistics import cpg
from PiSlice.popstatistics import gc
import numpy as np

def test_cpg():
    assert (cpg("ATTANNCGCGGCCCGATA") == 0.375)



def test_gc():
    test1 = "ATCGTAGCTGGGCTAGCTGATGCGCGCG"  # 18 GC on 28 nucleotides
    test2 = "ATCGTNGCNGGGCTAGCTGATGCGCGCG"  # 18 GC on 26 nucleotides + 2 Ns not to consider

    assert ((gc(test1) == 18 / 28) & (gc(test2) == 18 / 26) & (gc("ATC") is np.nan) & (gc("") is np.nan))

