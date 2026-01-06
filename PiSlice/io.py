"""
Input/output functions

Uses the polars-bio package and the polars framework as backend

Read FASTA, GFF, VCF and import them in genomic intervals 'polars-bio' objects
"""

import polars_bio as pb


class fasta():
    """
    Read a fasta file and return a 'fasta' sequence intervals object
    using the 'polars-bio' backend
    
    :param filename: string, the path to the 'fasta' file
    """
    def __init__(self, filename):
        f = pb.read_fasta(filename)
        self = f

    def summary(self):
        """
        Returns a summary of the fasta file. Number of sequences, names and sequence length.
        """

    def length(self):
        """
        Returns the length of each sequence (list)
        """
        out = [print(len(x)) for x in self["sequence"]]
        return out

    def seqname(self):
        """
        Returns sequence names of a fasta object (list)
        """
        return list(self.seq.keys())

    def sequence(self):
        """
        Returns sequences of a fasta object (list)
        """
        return list(self.seq.values())

    def sample_chromosome(self, chromosome):
        """
        Sample a whole chromosome sequence
        :param chromosome: str, the name of the chromosome
        :return: A DNA sequence of the whole chromosome (str)
        """
        # Chromosome can be either an integer (index) or a string (name)
        return self.seq[str(chromosome)]

    def sample_sequence(self, chromosome, start, end):
        """
        Sample a DNA sequence on a chromosome within start-end coordinates (bp)
        on a 1-offset (sequence begins at 1)
        :param chromosome: str, the name of the chromosome
        :param start: int, start position in the sequence
        :param end: int, end position in the sequence
        :return: A DNA sequence within boundaries start-end on a given chromosome (str)
        """
        # Chromosome can be either an integer (index) or a string (name)
        # start - 1 because we work in 1-bp offset while python is 0-bp offset
        # Verify that start-end positions are not inverted (start < end)
        start = min(start, end)
        end = max(start, end)
        if end > start:
            res = self.seq[str(chromosome)][(start - 1):end:]
        else:
            res = np.NaN
        return res

    def sample_sequence_masked(self, chromosome, start, end, mask):
        """
        Sample a DNA sequence on a chromosome within start-end coordinates (bp)
        on a 1-offset (sequence begins at 1) and mask regions given in 'mask'
        :param chromosome: str, the name of the chromosome
        :param start: int, start position in the sequence
        :param end: int, end position in the sequence
        :param mask: int, a list of tuples giving start-end coordinates of the region to mask
        :return: A DNA sequence within boundaries start-end on a given chromosome (str)
        """
        # Expand masked intervals
        mask = [(x - 1, y + 1) for x, y in mask]
        # Take care of Null interval objects
        if end > start:
            coord = intervaltree.IntervalTree.from_tuples([(start, end)])
            [coord.chop(x, y) for x, y in tuple(mask) if x != y]
            seq = [self.sample_sequence(chromosome, x[0], x[1]) for x in coord.items()]
        else:
            seq = ""

        return seq

    def sample_feature(self, gff, feature, chromosome, start, end):
        """
        Sample a list of DNA sequences corresponding to features (e.g. CDS)
        within a chromosome or start-end boundaries
        :param gff: gff, a Data Frame with genome annotation
        :param feature: str, the type of feature
        :param chromosome: str, the name of the chromosome
        :param start: int, start position in the sequence
        :param end: int, end position in the sequence
        :return: str, list of DNA sequences
        """