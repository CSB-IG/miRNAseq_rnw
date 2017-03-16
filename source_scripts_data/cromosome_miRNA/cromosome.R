

library(Gviz)


start <- scan(file="4.txt", what="")
start <- as.numeric(start)
end <- scan(file="5.txt", what="")
end <- as.numeric(end)
mir <- scan(file="7.txt", what="")
strand <- scan(file="6.txt", what="")

#######

# from = 100109655, to = 101560422
pdf("1.pdf", width = 10, height = 1)
ideoTrack <- IdeogramTrack(genome = "hg19", chromosome = "chr14")
plotTracks(ideoTrack, from = 100000000, to = 101571000, showBandId = TRUE, cex.bands = 0.5)
dev.off()

pdf("2.pdf", width = 10, height = 1)
axisTrack <- GenomeAxisTrack()
plotTracks(axisTrack, from = 100000000, to = 101571000, add53 = TRUE, add35 = TRUE, littleTicks = TRUE)
dev.off()

pdf("7.pdf", width = 10, height = 5)
E <- AnnotationTrack(range = NULL, start= start, end = end, group = mir, strand = strand, chromosome = "chr14", genome = "hg19", name="AlignedReadTrack")
plotTracks(E, groupAnnotation = "group", from = 100000000, to = 101571000)
dev.off()

a <-c(ideoTrack, axisTrack, E)

pdf("8.pdf", width = 10, height = 5)
plotTracks(a, groupAnnotation = "group", from = 100000000, to = 101571000)
dev.off()

pdf("9.pdf", width = 10, height = 5)
grtrack <- GeneRegionTrack(geneModels, genome = "hg19", chromosome = "chr14", name = "Gene Model", transcriptAnnotation = "symbol", background.title = "brown")
plotTracks(grtrack, from = 100000000, to = 102571000)
dev.off()


chr14:100210000–101150000


pdf("10.pdf", width = 10, height = 5)
biomTrack <- BiomartGeneRegionTrack(genome = "hg19", chromosome = "chr14", start = 100000000, end = 102571000, name = "ENSEMBL", symbol = c("DLK1", "MEG3", "RTL1", "DIO3"), )
plotTracks(biomTrack)
dev.off()





