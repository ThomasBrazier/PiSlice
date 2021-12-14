#============================================================================#
# Loading environment ----
#============================================================================#
# clear global environment: remove all variables
rm(list=ls(all=TRUE))
# Loading packages
library(rstudioapi)
library(ggplot2)
# Get the directory of the file & set working directory
wd=dirname(rstudioapi::getSourceEditorContext()$path)
setwd(wd)

#============================================================================#
# Look assembly reports for chromosome names and length ----
#============================================================================#
# *_assembly_report.txt in ncbi directory
# Sequence-Name	Sequence-Role	Assigned-Molecule	Assigned-Molecule-Location/Type	GenBank-Accn	Relationship	RefSeq-Accn	Assembly-Unit	Sequence-Length	UCSC-style-name

#============================================================================#
# Import gff ----
#============================================================================#
gff_arabidopsis = read.table(gzfile("data/Arabidopsis_thaliana_GCA_000001735.2/GCA_000001735.2_TAIR10.1_genomic.gff.gz"),
           header = FALSE, sep ="\t", comment.char = "#",
           col.names = c("seqid", "source", "feature", "start", "end", "score", "strand", "phase", "attributes"))
unique(gff_arabidopsis$seqid)
# Trim chromosomes
gff_arabidopsis = subset(gff_arabidopsis, gff_arabidopsis$seqid %in% c("CP002684.1", "CP002685.1", "CP002686.1", "CP002687.1", "CP002688.1"))
unique(gff_arabidopsis$seqid)


gff_oryza= read.table(gzfile("data/Oryza_sativa_GCA_001433935.1/GCA_001433935.1_IRGSP-1.0_genomic.gff.gz"),
                             header = FALSE, sep ="\t", comment.char = "#",
                      col.names = c("seqid", "source", "feature", "start", "end", "score", "strand", "phase", "attributes"))
unique(gff_oryza$seqid)
# Trim chromosomes
gff_oryza = subset(gff_oryza, gff_oryza$seqid %in% c("AP014957.1", "AP014958.1", "AP014959.1", "AP014960.1", "AP014961.1",
                                                     "AP014962.1", "AP014963.1", "AP014964.1", "AP014965.1", "AP014966.1",
                                                     "AP014967.1", "AP014968.1"))
unique(gff_oryza$seqid)


#============================================================================#
# Extract 100kb windows ----
#============================================================================#
gff = gff_arabidopsis
gff = gff_oryza
# Get list of chromosomes
chr = unique(gff$seqid)
df = data.frame(chromosome = character(), start = integer(), end = integer())
for (c in chr) {
  start = seq(1, max(gff$end[gff$seqid == c], na.rm = TRUE), 100000)
  end = start[-1] - 1
  start = start[-length(start)]
  df = rbind(df, data.frame(c, start, end))
}

write.table(df, file = "data/Arabidopsis_thaliana_GCA_000001735.2/Athaliana_100kb.csv",
            col.names = F, row.names = F, quote = F, sep = "\t")
write.table(df, file = "data/Oryza_sativa_GCA_001433935.1/Osativa_100kb.csv",
            col.names = F, row.names = F, quote = F, sep = "\t")

#============================================================================#
# Extract CDS first exon ----
#============================================================================#
gff = gff_arabidopsis
gff = gff_oryza
# Get list of chromosomes
chr = unique(gff$seqid)
df = data.frame(chromosome = character(), start = integer(), end = integer())
for (c in chr) {
  # Get only CDS first exon
  
  start = seq(1, max(gff$end[gff$seqid == c], na.rm = TRUE), 100000)
  end = start[-1] - 1
  start = start[-length(start)]
  df = rbind(df, data.frame(c, start, end))
}

write.table(df, file = "data/Arabidopsis_thaliana_GCA_000001735.2/Athaliana_CDSexon1.csv",
            col.names = F, row.names = F, quote = F, sep = "\t")
write.table(df, file = "data/Oryza_sativa_GCA_001433935.1/Osativa_CDSexon1.csv",
            col.names = F, row.names = F, quote = F, sep = "\t")



