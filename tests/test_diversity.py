def test_snp_count():
    import allel
    import PiSlice.diversity as div
    vcf_path = "test.vcf.gz"
    vcf = allel.read_vcf(vcf_path)
    chrom = "1"
    start = 1000
    stop = 20000
    snp_number = div.snp_count(vcf, chrom, start, stop)
    assert(snp_number == 65)

def test_snp_count():
    import allel
    import PiSlice.diversity as div
    vcf_path = "test.vcf.gz"
    vcf = allel.read_vcf(vcf_path)
    chrom = "1"
    start = 1000
    stop = 20000
    snp_density = div.snp_density(vcf, chrom, start, stop)
    assert(snp_density == 0.0034210526315789475)
