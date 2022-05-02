from pysam import FastaFile
import gzip
from Bio import SeqIO
import time
import statistics
import numpy as np

# Set of parameters
nboot = 1000
chromosome = "CP002684.1"
start = 100
end = 5000
# Use Pysam or dict structure for fasta sequences?
filename = "PiSlice/data/fasta_test.fna.bgz"
genome = FastaFile(filename)

boot = np.zeros(nboot)
for i in range(nboot):
    #print(i)
    t1 = time.time()
    a = genome.fetch(chromosome, start - 1, end)
    t2 = time.time()
    t = t2 - t1
    boot[i] = t
statistics.mean(boot)


# Dict
filename = "PiSlice/data/fasta_test.fna.gz"
f = gzip.open(filename, "rt")
seq = {rec.id: str(rec.seq) for rec in SeqIO.parse(f, "fasta")}

boot = np.zeros(nboot)
for i in range(nboot):
    #print(i)
    t1 = time.time()
    a = seq[str(chromosome)][start - 1:end:]
    t2 = time.time()
    t = t2 - t1
    boot[i] = t
statistics.mean(boot)

# Dict is 100 times faster than Pysam
# See also https://peerj.com/preprints/970v1/