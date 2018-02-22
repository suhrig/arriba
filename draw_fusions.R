#!/usr/bin/env Rscript

# parse command-line parameters
args <- commandArgs(trailingOnly=T)
parseBooleanParameter <- function(parameter, args, default) {
	arg <- sub(paste0("^--", parameter, "="), "", args[which(grepl(paste0("^--", parameter, "="), args))], perl=T)
	if (length(arg) == 0) {
		return(default)
	} else if (arg == "TRUE" || arg == "T") {
		return(T)
	} else if (arg == "FALSE" || arg == "F") {
		return(F)
	} else {
		stop(paste0("Invalid argument to --", parameter))
	}
}
parseStringParameter <- function(parameter, args, default="") {
	result <- sub(paste0("^--", parameter, "="), "", args[tail(which(grepl(paste0("^--", parameter, "="), args)), 1)], perl=T)
	return(ifelse(length(result) == 0, default, result))
}
if (any(grepl("^--help", args)))
	stop("usage: draw_fusions.R --annotation=annotation.gtf --fusions=fusions.tsv --output=output.pdf [--alignments=Aligned.out.bam] [--assemblyVersion=hg19] [--squishIntrons=TRUE] [--printExonLabels=TRUE] [--printStats=TRUE] [--printFusionTranscript=TRUE] [--pdfPaper=a4r] [--pdfWidth=11] [--pdfHeight=7]")
exonsFile <- parseStringParameter("annotation", args)
if (file.access(exonsFile) == -1)
	stop(sprintf("Exon annotation file (%s) does not exist", exonsFile))
fusionsFile <- parseStringParameter("fusions", args)
if (file.access(fusionsFile) == -1)
	stop(sprintf("Fusions file (%s) does not exist", fusionsFile))
outputFile <- parseStringParameter("output", args)
if (outputFile == "")
	stop("Output file not specified")
alignmentsFile <- parseStringParameter("alignments", args)
if (alignmentsFile != "") {
	if (file.access(alignmentsFile) == -1)
		stop(sprintf("Alignments file (%s) does not exist", alignmentsFile))
	if (!suppressPackageStartupMessages(require(GenomicAlignments)))
		stop("Package 'GenomicAlignments' must be installed when '--alignments' is used")
}
assemblyVersion <- parseStringParameter("assemblyVersion", args)
ideograms <- NULL
if (assemblyVersion != "") {
	if (!suppressPackageStartupMessages(require(rtracklayer)) || !suppressPackageStartupMessages(require(biovizBase)))
		stop("Packages 'rtracklayer' and 'biovizBase' must be installed when '--assemblyVersion' is used")
	ideograms <- as.data.frame(getIdeogram(assemblyVersion, cytoband = TRUE))
}
squishIntrons <- parseBooleanParameter("squishIntrons", args, T)
printExonLabels <- parseBooleanParameter("printExonLabels", args, T)
printStats <- parseBooleanParameter("printStats", args, T)
printFusionTranscript <- parseBooleanParameter("printFusionTranscript", args, T)
pdfPaper <- parseStringParameter("pdfPaper", args, "a4r")
pdfWidth <- as.integer(parseStringParameter("pdfWidth", args, "11"))
pdfHeight <- as.integer(parseStringParameter("pdfHeight", args, "7"))

# read fusions
fusions <- read.table(fusionsFile, stringsAsFactors=F, sep="\t", header=T, comment.char="", quote="")
colnames(fusions)[colnames(fusions) == "X.gene1"] <- "gene1"
fusions$contig1 <- sub(":.*", "", fusions$breakpoint1)
fusions$breakpoint1 <- as.numeric(sub(".*:", "", fusions$breakpoint1, perl=T))
fusions$contig2 <- sub(":.*", "", fusions$breakpoint2)
fusions$breakpoint2 <- as.numeric(sub(".*:", "", fusions$breakpoint2, perl=T))

pdf(outputFile, onefile=T, paper=pdfPaper, width=pdfWidth, height=pdfHeight, title=fusionsFile)
if (nrow(fusions) == 0) {
	plot(0, 0, type="l", xaxt="n", yaxt="n", xlab="", ylab="")
	text(0, 0, "Error: empty input file\n")
	dev.off()
	quit("no")
}

# convenience functions to add/remove "chr" prefix
addChr <- function(contig) {
	ifelse(contig == "MT", "chrM", paste0("chr", contig))
}
removeChr <- function(contig) {
	sub("chr", "", sub("chrM", "MT", contig))
}

# read exon annotation
message("Loading annotation")
exons <- read.table(exonsFile, header=F, sep="\t", comment.char="#", quote="", stringsAsFactors=F)[,c(1, 3, 4, 5, 7, 9)]
colnames(exons) <- c("contig", "type", "start", "end", "strand", "attributes")
exons <- exons[exons$type %in% c("exon", "CDS"),]
exons$geneName <- gsub(".*gene_name \"?([^;\"]+)\"?;.*", "\\1", exons$attributes)
exons$transcript <- gsub(".*transcript_id \"?([^;\"]+)\"?;.*", "\\1", exons$attributes)
exons$exonNumber <- ifelse(printExonLabels & grepl("exon_number ", exons$attributes), gsub(".*exon_number \"?([^;\"]+)\"?;.*", "\\1", exons$attributes), "")
exons$contig <- removeChr(exons$contig)

# prepare ideogram data
if (!is.null(ideograms)) {
	ideograms$seqnames <- removeChr(as.character(ideograms$seqnames))
	ideograms$gieStain <- as.character(ideograms$gieStain)
	ideograms <- ideograms[order(ideograms$seqnames, ideograms$start, ideograms$end),]
}

# insert dummy annotations for dummy genes
if (any(grepl(",", fusions$gene1) | grepl(",", fusions$gene2))) {
	intergenicBreakpoints <- rbind(
		setNames(fusions[grepl(",", fusions$gene1),c("gene1", "contig1", "breakpoint1")], c("gene", "contig", "breakpoint")),
		setNames(fusions[grepl(",", fusions$gene2),c("gene2", "contig2", "breakpoint2")], c("gene", "contig", "breakpoint"))
	)
	exons <- rbind(exons, data.frame(
		contig=intergenicBreakpoints$contig,
		type="exon",
		start=intergenicBreakpoints$breakpoint-5000,
		end=intergenicBreakpoints$breakpoint+5000,
		strand="",
		attributes="",
		geneName=intergenicBreakpoints$gene,
		transcript=intergenicBreakpoints$gene,
		exonNumber="intergenic"
	))
}
exons <- unique(exons[order(exons$start),])

drawCurlyBrace <- function(left, right, top, bottom, tip) {
	arcSize <- (top-bottom)/2
	lines(left+(1+cos(seq(pi/2,pi,len=10)))*arcSize, bottom+sin(seq(pi/2,pi,len=10))*arcSize)
	lines(c(left+arcSize, tip-arcSize), c(bottom+arcSize, bottom+arcSize))
	lines(tip+(cos(seq(0,pi/2,len=10))-1)*arcSize, top-sin(seq(0,pi/2,len=10))*arcSize)
	lines(tip+(1+cos(seq(pi/2,pi,len=10)))*arcSize, top-sin(seq(pi/2,pi,len=10))*arcSize)
	lines(c(tip+arcSize, right-arcSize), c(bottom+arcSize, bottom+arcSize))
	lines(right+(cos(seq(0,pi/2,len=10))-1)*arcSize, bottom+sin(seq(0,pi/2,len=10))*arcSize)
}

drawIdeogram <- function(left, right, top, bottom, ideograms, contig, breakpoint, scale) {
	# define design of ideogram
	bandColors <- c(gneg="#ffffff", gpos25="#bbbbbb", gpos50="#888888", gpos75="#444444", gpos100="#000000", acen="#ff0000", stalk="#0000ff")
	ideograms$color <- bandColors[ideograms$gieStain]
	curlyBraceHeight <- 0.03
	# extract bands of given contig
	bands <- ideograms[ideograms$seqnames==contig,]
	if (nrow(bands) == 0)
		stop(paste("Giemsa bands of contig", contig, "not found"))
	# scale width of ideogram to fit inside given region
	bands$left <- (bands$start - min(bands$start)) / max(bands$end) * (right-left-0.01-curlyBraceHeight*2) * scale
	bands$right <- (bands$end - min(bands$start)) / max(bands$end) * (right-left-0.01-curlyBraceHeight*2) * scale
	# center ideograms
	bands$left <- (right+left-max(bands$right))/2 + bands$left
	bands$right <- (right+left-max(bands$right))/2 + bands$right
	# draw curly braces
	tip <- min(bands$left) + (max(bands$right)-min(bands$left)) / (max(bands$end)-min(bands$end)) * breakpoint
	drawCurlyBrace(left, right, bottom+curlyBraceHeight, bottom, tip)
	# draw title of chromosome
	text((max(bands$right)+min(bands$left))/2, top+0.04, paste("chromosome", contig))
	# draw name of band
	bandName <- bands[which(bands$start <= breakpoint & bands$end >= breakpoint), "name"]
	text(tip, top+0.01, bandName, cex=0.75)
	# draw start of chromosome
	polygon(
		bands[1,"left"] + (1+cos(seq(pi/2,1.5*pi,len=10))) * (bands[1,"right"]-bands[1,"left"]),
		(top+bottom+curlyBraceHeight+0.01)/2 + sin(seq(pi/2,1.5*pi,len=10)) * (top-bottom-curlyBraceHeight-0.01)/2,
		col=bands[1,"color"]
	)
	# draw bands
	centromereStart <- T
	for (band in 2:(nrow(bands)-1)) {
		if (bands[band,"gieStain"] != "acen") {
			rect(bands[band,"left"], bottom+curlyBraceHeight+0.01, bands[band,"right"], top, col=bands[band,"color"])
		} else { # draw centromere
			if (centromereStart) {
				polygon(c(bands[band,"left"], bands[band,"right"], bands[band,"left"]), c(bottom+curlyBraceHeight+0.01, (top+bottom+curlyBraceHeight+0.01)/2, top), col=bands[band,"color"])
				centromereStart <- F
			} else {
				polygon(c(bands[band,"right"], bands[band,"left"], bands[band,"right"]), c(bottom+curlyBraceHeight+0.01, (top+bottom+curlyBraceHeight+0.01)/2, top), col=bands[band,"color"])
			}
		}
	}
	# draw end of chromosome
	band <- nrow(bands)
	polygon(
		bands[band,"right"] - (1+cos(seq(1.5*pi,pi/2,len=10))) * (bands[band,"right"]-bands[band,"left"]),
		(top+bottom+curlyBraceHeight+0.01)/2 + sin(seq(1.5*pi,pi/2,len=10)) * (top-bottom-curlyBraceHeight-0.01)/2,
		col=bands[band,"color"]
	)
}

drawCoverage <- function(left, right, y, coverage, start, end) {
	# draw coverage as bars
	if (!is.null(coverage)) {
		coverageData <- as.numeric(coverage[IRanges(start, end)])
		for (position in 1:length(coverageData))
			rect(left+(position-1)/(end-start)*(right-left), y, left+position/(end-start)*(right-left), y+coverageData[position]*0.1, col="grey", border=NA)
	}
}

drawStrands <- function(left, right, y, col, strand) {
	# draw strand
	lines(c(left, right), c(y, y), col=col, lwd=2)
	# indicate orientation
	if (strand %in% c("+", "-")) {
		if (right - left > 0.01)
			for (i in seq(left+0.01, right-0.01, by=sign(right-left-2*0.01)*0.01))
				arrows(i, y, i+0.001*ifelse(strand=="+", 1, -1), y, col=col, length=0.08, lwd=2, angle=60)
	}
}

drawExon <- function(left, right, y, color, title, type) {
	if (title != "dummy") {
		if (type == "CDS") {
			# draw coding regions as thicker bars
			rect(left, y+0.05, right, y+0.025, col=color, border=NA)
			rect(left, y-0.05, right, y-0.025, col=color, border=NA)
		} else {
			rect(left, y+0.025, right, y-0.025, col=color, border=NA)
			text((left+right)/2, y, title, cex=0.75)
		}
	}
}

findExons <- function(exons, gene, direction, contig, breakpoint) {
	# look for exon with breakpoint as splice site
	transcripts <- exons[exons$geneName == gene & exons$contig == contig & exons$type == "exon" & (direction == "downstream" & abs(exons$end - breakpoint) <= 2 | direction == "upstream" & abs(exons$start - breakpoint) <= 2),"transcript"]
	candidateExons <- exons[exons$transcript %in% transcripts,]
	# if none was found search for transcripts which encompass the breakpoint
	if (nrow(candidateExons) == 0) {
		candidateExons <- exons[exons$geneName == gene & exons$contig == contig,]
		if (nrow(candidateExons) > 0) {
			transcriptStart <- aggregate(candidateExons$start, by=list(candidateExons$transcript), min)
			rownames(transcriptStart) <- transcriptStart[,1]
			transcriptEnd <- aggregate(candidateExons$end, by=list(candidateExons$transcript), max)
			rownames(transcriptEnd) <- transcriptEnd[,1]
			candidateExons <- candidateExons[transcriptStart[candidateExons$transcript,2] <= breakpoint & transcriptEnd[candidateExons$transcript,2] >= breakpoint,]
		}
	}
	# find the consensus transcript, if there are multiple hits
	if (length(unique(candidateExons$transcript)) > 1) {
		consensusTranscript <-
			ifelse(grepl("appris_principal_1", candidateExons$attributes), 10,
			ifelse(grepl("appris_principal_2", candidateExons$attributes), 9,
			ifelse(grepl("appris_principal_3", candidateExons$attributes), 8,
			ifelse(grepl("appris_principal_4", candidateExons$attributes), 7,
			ifelse(grepl("appris_principal_5", candidateExons$attributes), 6,
			ifelse(grepl("appris_principal", candidateExons$attributes), 5,
			ifelse(grepl("appris_alternative_1", candidateExons$attributes), 4,
			ifelse(grepl("appris_alternative_2", candidateExons$attributes), 3,
			ifelse(grepl("appris_alternative", candidateExons$attributes), 2,
			ifelse(grepl("CCDS", candidateExons$attributes), 1,
			0
		))))))))))
		candidateExons <- candidateExons[consensusTranscript == max(consensusTranscript),]
	}
	# use the transcript with the longest coding sequence, if there are still multiple hits
	if (length(unique(candidateExons$transcript)) > 1) {
		codingSequenceLength <- ifelse(candidateExons$type == "CDS", candidateExons$end - candidateExons$start, 0)
		totalCodingSequenceLength <- aggregate(codingSequenceLength, by=list(candidateExons$transcript), sum)
		rownames(totalCodingSequenceLength) <- totalCodingSequenceLength[,1]
		candidateExons <- candidateExons[totalCodingSequenceLength[candidateExons$transcript,2] == max(totalCodingSequenceLength[,2]),]
	}
	# use the transcript with the longest overall sequence, if there are still multiple hits
	if (length(unique(candidateExons$transcript)) > 1) {
		exonLength <- candidateExons$end - candidateExons$start
		totalExonLength <- aggregate(exonLength, by=list(candidateExons$transcript), sum)
		rownames(totalExonLength) <- totalExonLength[,1]
		candidateExons <- candidateExons[totalExonLength[candidateExons$transcript,2] == max(totalExonLength[,2]),]
	}
	# if there are still multiple hits, select the first one
	candidateExons <- candidateExons[candidateExons$transcript == head(unique(candidateExons$transcript), 1),]
	return(candidateExons)
}

for (fusion in 1:nrow(fusions)) {

	message(paste0("Drawing fusion #", fusion, ": ", fusions[fusion,"gene1"], ":", fusions[fusion,"gene2"]))

	exons1 <- findExons(exons, fusions[fusion,"gene1"], fusions[fusion,"direction1"], fusions[fusion,"contig1"], fusions[fusion,"breakpoint1"])
	if (nrow(exons1) == 0) {
		plot(0, 0, type="l", xaxt="n", yaxt="n", xlab="", ylab="")
		text(0, 0, paste0("Error: exon coordinates of ", fusions[fusion,"gene1"], " not found in\n", exonsFile))
		next
	}
	exons2 <- findExons(exons, fusions[fusion,"gene2"], fusions[fusion,"direction2"], fusions[fusion,"contig2"], fusions[fusion,"breakpoint2"])
	if (nrow(exons2) == 0) {
		plot(0, 0, type="l", xaxt="n", yaxt="n", xlab="", ylab="")
		text(0, 0, paste0("Error: exon coordinates of ", fusions[fusion,"gene2"], " not found in\n", exonsFile))
		next
	}

	# compute coverage from alignments file
	coverage1 <- NULL
	coverage2 <- NULL
	if (alignmentsFile != "") {
		# function which reads alignments from BAM file with & without "chr" prefix
		readCoverage <- function(alignmentsFile, contig, start, end) {
			coverageData <- tryCatch(
				{
					alignments <- readGAlignments(alignmentsFile, param=ScanBamParam(which=GRanges(contig, IRanges(start, end))))
					coverage(alignments)[[contig]]
				},
				error=function(e) {
					alignments <- readGAlignments(alignmentsFile, param=ScanBamParam(which=GRanges(addChr(contig), IRanges(start, end))))
					coverage(alignments)[[addChr(contig)]]
				}
			)
			if (exists("alignments")) rm(alignments)
			return(coverageData)
		}
		# get coverage track
		coverage1 <- readCoverage(alignmentsFile, fusions[fusion,"contig1"], min(exons1$start), max(exons1$end))
		coverage2 <- readCoverage(alignmentsFile, fusions[fusion,"contig2"], min(exons2$start), max(exons2$end))
		# normalize coverage
		coverageNormalization <- max(max(coverage1), max(coverage2))
		coverage1 <- coverage1/coverageNormalization
		coverage2 <- coverage2/coverageNormalization
	}

	# insert dummy exons, if breakpoints are outside the gene (e.g., in UTRs)
	# this avoids plotting artifacts
	breakpoint1 <- fusions[fusion,"breakpoint1"]
	breakpoint2 <- fusions[fusion,"breakpoint2"]
	if (breakpoint1 < min(exons1$start)) {
		exons1 <- rbind(c(exons1[1,"contig"], "dummy", breakpoint1-1000, breakpoint1-1000, exons1[1,"strand"], "", "dummy", exons1[1,"transcript"], ""), exons1)
	} else if (breakpoint1 > max(exons1$end)) {
		exons1 <- rbind(exons1, c(exons1[1,"contig"], "dummy", breakpoint1+1000, breakpoint1+1000, exons1[1,"strand"], "", "dummy", exons1[1,"transcript"], ""))
	}
	if (breakpoint2 < min(exons2$start)) {
		exons2 <- rbind(c(exons2[1,"contig"], "dummy", breakpoint2-1000, breakpoint2-1000, exons2[1,"strand"], "", "dummy", exons2[1,"transcript"], ""), exons2)
	} else if (breakpoint2 > max(exons2$end)) {
		exons2 <- rbind(exons2, c(exons2[1,"contig"], "dummy", breakpoint2+1000, breakpoint2+1000, exons2[1,"strand"], "", "dummy", exons2[1,"transcript"], ""))
	}
	exons1$start <- as.integer(exons1$start)
	exons1$end <- as.integer(exons1$end)
	exons2$start <- as.integer(exons2$start)
	exons2$end <- as.integer(exons2$end)

	exons1$left <- exons1$start
	exons1$right <- exons1$end
	exons2$left <- exons2$start
	exons2$right <- exons2$end

	if (squishIntrons) {
		# hide introns in gene1
		cumulativeIntronLength <- 0
		previousExonEnd <- -200
		for (exon in 1:nrow(exons1)) {
			if (breakpoint1 > previousExonEnd+1 && breakpoint1 < exons1[exon,"left"])
				breakpoint1 <- (breakpoint1-previousExonEnd) / (exons1[exon,"left"]-previousExonEnd) * 200 + previousExonEnd - cumulativeIntronLength
			if (exons1[exon,"left"] > previousExonEnd)
				cumulativeIntronLength <- cumulativeIntronLength + exons1[exon,"left"] - previousExonEnd - 200
			if (breakpoint1 >= exons1[exon,"left"] && breakpoint1 <= exons1[exon,"right"]+1)
				breakpoint1 <- breakpoint1 - cumulativeIntronLength
			previousExonEnd <- exons1[exon,"right"]
			exons1[exon,"left"] <- exons1[exon,"left"] - cumulativeIntronLength
			exons1[exon,"right"] <- exons1[exon,"right"] - cumulativeIntronLength
		}

		# hide introns in gene2
		cumulativeIntronLength <- 0
		previousExonEnd <- -200
		for (exon in 1:nrow(exons2)) {
			if (breakpoint2 > previousExonEnd+1 && breakpoint2 < exons2[exon,"left"])
				breakpoint2 <- (breakpoint2-previousExonEnd) / (exons2[exon,"left"]-previousExonEnd) * 200 + previousExonEnd - cumulativeIntronLength
			if (exons2[exon,"left"] > previousExonEnd)
				cumulativeIntronLength <- cumulativeIntronLength + exons2[exon,"left"] - previousExonEnd - 200
			if (breakpoint2 >= exons2[exon,"left"] && breakpoint2 <= exons2[exon,"right"]+1)
				breakpoint2 <- breakpoint2 - cumulativeIntronLength
			previousExonEnd <- exons2[exon,"right"]
			exons2[exon,"left"] <- exons2[exon,"left"] - cumulativeIntronLength
			exons2[exon,"right"] <- exons2[exon,"right"] - cumulativeIntronLength
		}
	} else { # don't squish introns
		# shift exon coordinates to align the gene to the left border of the plot
		exons1$right <- exons1$right - min(exons1$left)
		breakpoint1 <- breakpoint1 - min(exons1$left)
		exons1$left <- exons1$left - min(exons1$left)
		exons2$right <- exons2$right - min(exons2$left)
		breakpoint2 <- breakpoint2 - min(exons2$left)
		exons2$left <- exons2$left - min(exons2$left)
	}

	# scale exon sizes to 1
	scalingFactor <- max(exons1$right) + max(exons2$right)
	exons1$left <- exons1$left / scalingFactor
	exons1$right <- exons1$right / scalingFactor
	exons2$left <- exons2$left / scalingFactor
	exons2$right <- exons2$right / scalingFactor
	breakpoint1 <- breakpoint1 / scalingFactor
	breakpoint2 <- breakpoint2 / scalingFactor

	# shift gene2 to the right of gene1 with a little bit of padding
	gene2Offset <- max(exons1$right)+0.05

	# center fusion horizontally
	fusionOffset1 <- (max(exons1$right)+gene2Offset)/2 - ifelse(fusions[fusion,"direction1"] == "downstream", breakpoint1, max(exons1$right)-breakpoint1)
	fusionOffset2 <- fusionOffset1 + ifelse(fusions[fusion,"direction1"] == "downstream", breakpoint1, max(exons1$right)-breakpoint1)
	
	par(mar=c(0, 0, 0, 0))
	plot(0, 0, type="l", xlim=c(-0.12, 1.12), ylim=c(-0.1, 1), bty="n", xaxt="n", yaxt="n")

	# vertical coordinates of layers
	yIdeograms <- 0.95
	yGeneNames <- 0.5
	yBreakpointLabels <- 0.8
	yCoverage <- 0.675
	yExons <- 0.6
	yFusion <- 0.3
	yStats <- 0

	# draw ideograms
	if (!is.null(ideograms)) {
		# scale ideograms such that the chromosomes fit in the region of the transcript
		contig1Size <- max(ideograms[ideograms$seqnames==fusions[fusion,"contig1"], "end"])
		contig2Size <- max(ideograms[ideograms$seqnames==fusions[fusion,"contig2"], "end"])
		transcript1Size <- max(exons1$right)
		transcript2Size <- max(exons2$right)
		if (transcript1Size/contig1Size < transcript2Size/contig2Size) {
			scaleIdeogram1 <- 1
			scaleIdeogram2 <- (transcript1Size/contig1Size) / (transcript2Size/contig2Size)
		} else {
			scaleIdeogram1 <- (transcript2Size/contig2Size) / (transcript1Size/contig1Size)
			scaleIdeogram2 <- 1
		}
		drawIdeogram(0, max(exons1$right), yIdeograms+0.04, yIdeograms-0.04, ideograms, fusions[fusion,"contig1"], fusions[fusion,"breakpoint1"], scaleIdeogram1)
		drawIdeogram(gene2Offset, gene2Offset+max(exons2$right), yIdeograms+0.04, yIdeograms-0.04, ideograms, fusions[fusion,"contig2"], fusions[fusion,"breakpoint2"], scaleIdeogram2)
	}

	# draw gene names
	text(max(exons1$right)/2, yGeneNames, fusions[fusion,"gene1"])
	text(gene2Offset+max(exons2$right)/2, yGeneNames, fusions[fusion,"gene2"])

	# label breakpoints
	text(breakpoint1, yBreakpointLabels, paste0("breakpoint\n", fusions[fusion,"contig1"], ":", fusions[fusion,"breakpoint1"]), cex=0.75)
	text(gene2Offset+breakpoint2, yBreakpointLabels, paste0("breakpoint\n", fusions[fusion,"contig2"], ":", fusions[fusion,"breakpoint2"]), cex=0.75)

	# draw coverage axis
	lines(c(-0.02, -0.01, -0.01, -0.02), c(yCoverage, yCoverage, yCoverage+0.1, yCoverage+0.1))
	text(-0.025, yCoverage, "0", adj=c(1,0.5), cex=0.75)
	text(-0.025, yCoverage+0.1, coverageNormalization, adj=c(1,0.5), cex=0.75)
	text(-0.06, yCoverage+0.05, "Coverage", srt=90, cex=0.75)

	# plot coverage 1
	if (squishIntrons) {
		for (exon in 1:nrow(exons1))
			drawCoverage(exons1[exon,"left"], exons1[exon,"right"], yCoverage, coverage1, exons1[exon,"start"], exons1[exon,"end"])
	} else {
		drawCoverage(min(exons1$left), max(exons1$right), yCoverage, coverage1, min(exons1$start), max(exons1$end))
	}

	# plot coverage 2
	if (squishIntrons) {
		for (exon in 1:nrow(exons2))
			drawCoverage(gene2Offset+exons2[exon,"left"], gene2Offset+exons2[exon,"right"], yCoverage, coverage2, exons2[exon,"start"], exons2[exon,"end"])
	} else {
		drawCoverage(gene2Offset+min(exons2$left), gene2Offset+max(exons2$right), yCoverage, coverage2, min(exons2$start), max(exons2$end))
	}

	# plot gene 1
	drawStrands(0, max(exons1$right), yExons, "darkolivegreen4", exons1[1,"strand"])
	for (exon in 1:nrow(exons1))
		drawExon(exons1[exon,"left"], exons1[exon,"right"], yExons, "darkolivegreen2", exons1[exon,"exonNumber"], exons1[exon,"type"])

	# plot gene 2
	drawStrands(gene2Offset, gene2Offset+max(exons2$right), yExons, col="deepskyblue4", exons2[1,"strand"])
	for (exon in 1:nrow(exons2))
		drawExon(gene2Offset+exons2[exon,"left"], gene2Offset+exons2[exon,"right"], yExons, "deepskyblue2", exons2[exon,"exonNumber"], exons2[exon,"type"])

	# plot gene1 of fusion
	if (fusions[fusion,"direction1"] == "downstream") {
		# plot strands
		drawStrands(fusionOffset1, fusionOffset1+breakpoint1, yFusion, col="darkolivegreen4", exons1[1,"strand"])
		# plot exons
		for (exon in 1:nrow(exons1))
			if (exons1[exon,"start"] <= fusions[fusion,"breakpoint1"])
				drawExon(fusionOffset1+exons1[exon,"left"], fusionOffset1+min(breakpoint1, exons1[exon,"right"]), yFusion, "darkolivegreen2", exons1[exon,"exonNumber"], exons1[exon,"type"])
		# plot trajectories
		lines(c(0, 0, fusionOffset1), c(yExons+0.05, yExons-0.05, yFusion+0.05), col="red", lty=2)
		lines(c(breakpoint1, breakpoint1, fusionOffset1+breakpoint1), c(yBreakpointLabels-0.02, yExons-0.05, yFusion+0.05), col="red", lty=2)
	} else if (fusions[fusion,"direction1"] == "upstream") {
		# plot strands
		drawStrands(fusionOffset1, fusionOffset2, yFusion, col="darkolivegreen4", chartr("+-", "-+", exons1[1,"strand"]))
		# plot exons
		for (exon in 1:nrow(exons1))
			if (exons1[exon,"end"]+1 >= fusions[fusion,"breakpoint1"])
				drawExon(fusionOffset1+max(exons1$right)-exons1[exon,"right"], min(fusionOffset2, fusionOffset1+max(exons1$right)-exons1[exon,"left"]), yFusion, "darkolivegreen2", exons1[exon,"exonNumber"], exons1[exon,"type"])
		# plot trajectories
		lines(c(max(exons1$right), max(exons1$right), fusionOffset1), c(yExons+0.05, yExons-0.05, yFusion+0.05), col="red", lty=2)
		lines(c(breakpoint1, breakpoint1, fusionOffset1+max(exons1$right)-breakpoint1), c(yBreakpointLabels-0.02, yExons-0.05, yFusion+0.05), col="red", lty=2)
	}
	
	# plot gene2 of fusion
	if (fusions[fusion,"direction2"] == "downstream") {
		# plot strands
		drawStrands(fusionOffset2, fusionOffset2+breakpoint2, yFusion, col="deepskyblue4", chartr("+-", "-+", exons2[1,"strand"]))
		# plot exons
		for (exon in 1:nrow(exons2))
			if (exons2[exon,"start"] <= fusions[fusion,"breakpoint2"])
				drawExon(max(fusionOffset2, fusionOffset2+breakpoint2-exons2[exon,"right"]), fusionOffset2+breakpoint2-exons2[exon,"left"], yFusion, "deepskyblue2", exons2[exon,"exonNumber"], exons2[exon,"type"])
		# plot trajectories
		lines(c(gene2Offset, gene2Offset, fusionOffset2+breakpoint2), c(yExons+0.05, yExons-0.05, yFusion+0.05), col="red", lty=2)
		lines(c(gene2Offset+breakpoint2, gene2Offset+breakpoint2, fusionOffset2), c(yBreakpointLabels-0.02, yExons-0.05, yFusion+0.05), col="red", lty=2)
	} else if (fusions[fusion,"direction2"] == "upstream") {
		# plot strands
		drawStrands(fusionOffset2, fusionOffset2+max(exons2$right)-breakpoint2, yFusion, col="deepskyblue4", exons2[1,"strand"])
		# plot exons
		for (exon in 1:nrow(exons2))
			if (exons2[exon,"end"]+1 >= fusions[fusion,"breakpoint2"])
				drawExon(max(fusionOffset2, fusionOffset2+exons2[exon,"left"]-breakpoint2), fusionOffset2+exons2[exon,"right"]-breakpoint2, yFusion, "deepskyblue2", exons2[exon,"exonNumber"], exons2[exon,"type"])
		# plot trajectories
		lines(c(gene2Offset+max(exons2$right), gene2Offset+max(exons2$right), fusionOffset2+max(exons2$right)-breakpoint2), c(yExons+0.05, yExons-0.05, yFusion+0.05), col="red", lty=2)
		lines(c(gene2Offset+breakpoint2, gene2Offset+breakpoint2, fusionOffset2), c(yBreakpointLabels-0.02, yExons-0.05, yFusion+0.05), col="red", lty=2)
	}
	
	if (printStats) {
		# print statistics about supporting alignments
		text(0.5, yStats, paste("Split reads in", fusions[fusion,"gene1"], "=", fusions[fusion,"split_reads1"]), cex=0.75)
		text(0.5, yStats-0.03, paste("Split reads in", fusions[fusion,"gene2"], "=", fusions[fusion,"split_reads2"]), cex=0.75)
		text(0.5, yStats-0.06, paste("Discordant mates =", fusions[fusion,"discordant_mates"]), cex=0.75)
	}

	if (printFusionTranscript && fusions[fusion,"fusion_transcript"] != ".") {
		# print fusion transcript colored by gene of origin
		fusion_transcript1 <- gsub("\\|.*", "", fusions[fusion,"fusion_transcript"], perl=T)
		fusion_transcript1 <- substr(fusion_transcript1, max(1, nchar(fusion_transcript1)-30), nchar(fusion_transcript1))
		fusion_transcript2 <- gsub(".*\\|", "", fusions[fusion,"fusion_transcript"], perl=T)
		fusion_transcript2 <- substr(fusion_transcript2, 1, min(nchar(fusion_transcript2), 30))
		text(0.5, yStats-0.09, bquote("Fusion transcript = " * phantom(.(fusion_transcript1)) * phantom(.(fusion_transcript2))), cex=0.75)
		text(0.5, yStats-0.09, bquote(phantom("Fusion transcript = ") * .(fusion_transcript1) * phantom(.(fusion_transcript2))), col="darkolivegreen4", cex=0.75)
		text(0.5, yStats-0.09, bquote(phantom("Fusion transcript = ") * phantom(.(fusion_transcript1)) * .(fusion_transcript2)), col="deepskyblue4", cex=0.75)
	}
}

devNull <- dev.off()
