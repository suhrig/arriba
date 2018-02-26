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
	stop("usage: draw_fusions.R --annotation=annotation.gtf --fusions=fusions.tsv --output=output.pdf [--alignments=Aligned.out.bam] [--species=hg19] [--printIdeograms=FALSE] [--printCircos=FALSE] [--minConfidenceForCircosPlot=medium] [--squishIntrons=TRUE] [--printExonLabels=TRUE] [--printStats=TRUE] [--printFusionTranscript=TRUE] [--pdfPaper=a4r] [--pdfWidth=11] [--pdfHeight=7] [--color1=#e5a5a5] [--color2=#a7c4e5]")
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
species <- parseStringParameter("species", args)
printIdeograms <- parseBooleanParameter("printIdeograms", args, F)
printCircos <- parseBooleanParameter("printCircos", args, F)
if (printCircos)
	if (!suppressPackageStartupMessages(require(circlize)))
		stop("Package 'circlize' must be installed when '--printCircos' is used")
minConfidenceForCircosPlot <- parseStringParameter("minConfidenceForCircosPlot", args, "medium")
if (!(minConfidenceForCircosPlot %in% c("low", "medium", "high")))
	stop("Invalid argument to --minConfidenceForCircosPlot")
squishIntrons <- parseBooleanParameter("squishIntrons", args, T)
printExonLabels <- parseBooleanParameter("printExonLabels", args, T)
printStats <- parseBooleanParameter("printStats", args, T)
printFusionTranscript <- parseBooleanParameter("printFusionTranscript", args, T)
pdfPaper <- parseStringParameter("pdfPaper", args, "a4r")
pdfWidth <- as.integer(parseStringParameter("pdfWidth", args, "11"))
pdfHeight <- as.integer(parseStringParameter("pdfHeight", args, "7"))
color1 <- parseStringParameter("color1", args, "#e5a5a5")
color2 <- parseStringParameter("color2", args, "#a7c4e5")


# get darker variants of colors
getDarkColor <- function(color) {
	rgb(
		max(0,col2rgb(color)["red",]-100),
		max(0,col2rgb(color)["green",]-100),
		max(0,col2rgb(color)["blue",]-100),
		maxColorValue=255
	)
}
darkColor1 <- getDarkColor(color1)
darkColor2 <- getDarkColor(color2)

# read fusions
fusions <- read.table(fusionsFile, stringsAsFactors=F, sep="\t", header=T, comment.char="", quote="")
colnames(fusions)[colnames(fusions) %in% c("X.gene1", "strand1.gene.fusion.", "strand2.gene.fusion.")] <- c("gene1", "strand1", "strand2")
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

# prepare ideogram data
ideograms <- NULL
if (printIdeograms || printCircos) {
	# try to load from file or download from UCSC using circlize
	if (file.access(species) != -1) {
		ideograms <- read.table(species, header=F)
	} else {
		if (!suppressPackageStartupMessages(require(circlize)))
			stop("Package 'circlize' must be installed when '--species' is used without a local file")
		ideograms <- read.cytoband(species=species)$df
	}
	colnames(ideograms) <- c("contig", "start", "end", "name", "giemsa")
	ideograms$contig <- removeChr(as.character(ideograms$contig))
	ideograms$giemsa <- as.character(ideograms$giemsa)
	ideograms <- ideograms[order(ideograms$contig, ideograms$start, ideograms$end),]
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

# insert dummy annotations for dummy genes
if (any(grepl(",", fusions$gene1) | grepl(",", fusions$gene2))) {
	intergenicBreakpoints <- rbind(
		setNames(fusions[grepl(",", fusions$gene1),c("gene1", "strand1", "contig1", "breakpoint1")], c("gene", "strand", "contig", "breakpoint")),
		setNames(fusions[grepl(",", fusions$gene2),c("gene2", "strand2", "contig2", "breakpoint2")], c("gene", "strand", "contig", "breakpoint"))
	)
	exons <- rbind(exons, data.frame(
		contig=intergenicBreakpoints$contig,
		type="exon",
		start=intergenicBreakpoints$breakpoint-5000,
		end=intergenicBreakpoints$breakpoint+5000,
		strand=sub(".*/", "", intergenicBreakpoints$strand, perl=T),
		attributes="",
		geneName=intergenicBreakpoints$gene,
		transcript=intergenicBreakpoints$gene,
		exonNumber="intergenic"
	))
}

drawVerticalGradient <- function(left, right, y, color, selection=NULL) {
	# check if gradient should only be drawn in part of the region
	if (!is.null(selection)) {
		y <- y[selection]
		left <- left[selection]
		right <- right[selection]
	}
	# draw gradient
	for (i in 1:length(y)) {
		polygon(
			c(left[1:i], right[1:i]),
			c(y[1:i], y[i:1]),
			border=NA,
			col=rgb(col2rgb(color)["red",], col2rgb(color)["green",], col2rgb(color)["blue",], col2rgb(color, alpha=T)["alpha",]*(1/length(y)), max=255)
		)
	}
}

drawCurlyBrace <- function(left, right, top, bottom, tip) {
	smoothness <- 20
	x <- cumsum(exp(-seq(-2.5, 2.5, len=smoothness)^2))
	x <- x/max(x)
	y <- seq(top, bottom, len=smoothness)
	lines(left+(tip-left)+x*(left-tip), y)
	lines(tip+x*(right-tip), y)
}

drawIdeogram <- function(adjust, left, right, y, ideograms, contig, breakpoint) {
	# define design of ideogram
	bandColors <- c(gneg="#ffffff", gpos25="#bbbbbb", gpos50="#888888", gpos75="#444444", gpos100="#000000", acen="#ec4f4f", stalk="#0000ff")
	ideograms$color <- bandColors[ideograms$giemsa]
	arcSteps <- 30 # defines roundness of arc
	curlyBraceHeight <- 0.03
	ideogramHeight <- 0.04
	ideogramWidth <- 0.4
	# extract bands of given contig
	bands <- ideograms[ideograms$contig==contig,]
	if (nrow(bands) == 0)
		stop(paste("Giemsa bands of contig", contig, "not found"))
	# scale width of ideogram to fit inside given region
	bands$left <- bands$start / max(ideograms$end) * ideogramWidth
	bands$right <- bands$end / max(ideograms$end) * ideogramWidth
	# left/right-align ideograms
	offset <- ifelse(adjust=="left", left, right - max(bands$right))
	bands$left <- bands$left + offset
	bands$right <- bands$right + offset
	# draw curly braces
	tip <- min(bands$left) + (max(bands$right)-min(bands$left)) / (max(bands$end)-min(bands$start)) * breakpoint
	drawCurlyBrace(left, right, y-0.05+curlyBraceHeight, y-0.05, tip)
	# draw title of chromosome
	text((max(bands$right)+min(bands$left))/2, y+0.07, paste("chromosome", contig), font=2)
	# draw name of band
	bandName <- bands[which(bands$start <= breakpoint & bands$end >= breakpoint), "name"]
	text(tip, y+0.03, bandName)
	# draw start of chromosome
	leftArcX <- bands[1,"left"] + (1+cos(seq(pi/2,1.5*pi,len=arcSteps))) * (bands[1,"right"]-bands[1,"left"])
	leftArcY <- y + sin(seq(pi/2,1.5*pi,len=arcSteps)) * (ideogramHeight/2)
	polygon(leftArcX, leftArcY, col=bands[1,"color"])
	# draw bands
	centromereStart <- NULL
	centromereEnd <- NULL
	for (band in 2:(nrow(bands)-1)) {
		if (bands[band,"giemsa"] != "acen") {
			rect(bands[band,"left"], y-ideogramHeight/2, bands[band,"right"], y+ideogramHeight/2, col=bands[band,"color"])
		} else { # draw centromere
			if (is.null(centromereStart)) {
				polygon(c(bands[band,"left"], bands[band,"right"], bands[band,"left"]), c(y-ideogramHeight/2, y, y+ideogramHeight/2), col=bands[band,"color"])
				centromereStart <- bands[band,"left"]
			} else {
				polygon(c(bands[band,"right"], bands[band,"left"], bands[band,"right"]), c(y-ideogramHeight/2, y, y+ideogramHeight/2), col=bands[band,"color"])
				centromereEnd <- bands[band,"right"]
			}
		}
	}
	# draw end of chromosome
	band <- nrow(bands)
	rightArcX <- bands[band,"right"] - (1+cos(seq(1.5*pi,pi/2,len=arcSteps))) * (bands[band,"right"]-bands[band,"left"])
	rightArcY <- y + sin(seq(pi/2,1.5*pi,len=arcSteps)) * ideogramHeight/2
	polygon(rightArcX, rightArcY, col=bands[band,"color"])
	# draw gradients for 3D effect
	drawVerticalGradient(leftArcX, rep(centromereStart, arcSteps), leftArcY, rgb(0,0,0,0.8), 1:round(arcSteps*0.4)) # black from top on p-arm
	drawVerticalGradient(leftArcX, rep(centromereStart, arcSteps), leftArcY, rgb(1,1,1,0.7), round(arcSteps*0.4):round(arcSteps*0.1)) # white to top on p-arm
	drawVerticalGradient(leftArcX, rep(centromereStart, arcSteps), leftArcY, rgb(1,1,1,0.7), round(arcSteps*0.4):round(arcSteps*0.6)) # white to bottom on p-arm
	drawVerticalGradient(leftArcX, rep(centromereStart, arcSteps), leftArcY, rgb(0,0,0,0.9), arcSteps:round(arcSteps*0.5)) # black from bottom on p-arm
	drawVerticalGradient(rightArcX, rep(centromereEnd, arcSteps), rightArcY, rgb(0,0,0,0.8), 1:round(arcSteps*0.4)) # black from top on q-arm
	drawVerticalGradient(rightArcX, rep(centromereEnd, arcSteps), rightArcY, rgb(1,1,1,0.7), round(arcSteps*0.4):round(arcSteps*0.1)) # white to top on q-arm
	drawVerticalGradient(rightArcX, rep(centromereEnd, arcSteps), rightArcY, rgb(1,1,1,0.7), round(arcSteps*0.4):round(arcSteps*0.6)) # white to bottom on q-arm
	drawVerticalGradient(rightArcX, rep(centromereEnd, arcSteps), rightArcY, rgb(0,0,0,0.9), arcSteps:round(arcSteps*0.5)) # black from bottom on q-arm
}

drawCoverage <- function(left, right, y, coverage, start, end, color) {
	# draw coverage as bars
	if (!is.null(coverage)) {
		coverageData <- as.numeric(coverage[IRanges(start, end)])
		for (position in 1:length(coverageData))
			rect(left+(position-1)/(end-start)*(right-left), y, left+position/(end-start)*(right-left), y+coverageData[position]*0.1, col=color, border=NA)
	}
}

drawStrand <- function(left, right, y, color, strand) {
	# draw strand
	lines(c(left+0.001, right-0.001), c(y, y), col=color, lwd=2)
	lines(c(left+0.001, right-0.001), c(y, y), col=rgb(1,1,1,0.1), lwd=1)
	# indicate orientation
	if (strand %in% c("+", "-")) {
		if (right - left > 0.01)
			for (i in seq(left+0.01, right-0.01, by=sign(right-left-2*0.01)*0.01)) {
				arrows(i, y, i+0.001*ifelse(strand=="+", 1, -1), y, col=color, length=0.05, lwd=2, angle=60)
				arrows(i, y, i+0.001*ifelse(strand=="+", 1, -1), y, col=rgb(1,1,1,0.1), length=0.05, lwd=1, angle=60)
			}
	}
}

drawExon <- function(left, right, y, color, title, type) {
	gradientSteps <- 10 # defines smoothenes of gradient
	exonHeight <- 0.03
	if (title != "dummy" && title != "intergenic") {
		if (type == "CDS") {
			# draw coding regions as thicker bars
			rect(left, y+exonHeight, right, y+exonHeight/2-0.001, col=color, border=NA)
			rect(left, y-exonHeight, right, y-exonHeight/2+0.001, col=color, border=NA)
			# draw border
			lines(c(left, left, right, right), c(y+exonHeight/2, y+exonHeight, y+exonHeight, y+exonHeight/2), col=getDarkColor(color))
			lines(c(left, left, right, right), c(y-exonHeight/2, y-exonHeight, y-exonHeight, y-exonHeight/2), col=getDarkColor(color))
			# draw gradients for 3D effect
			drawVerticalGradient(rep(left, gradientSteps), rep(right, gradientSteps), seq(y+0.03, y+0.015, len=gradientSteps), rgb(0,0,0,0.2))
			drawVerticalGradient(rep(left, gradientSteps), rep(right, gradientSteps), seq(y-0.03, y-0.015, len=gradientSteps), rgb(0,0,0,0.3))
		} else {
			rect(left, y+exonHeight/2, right, y-exonHeight/2, col=color, border=getDarkColor(color))
			# draw gradients for 3D effect
			drawVerticalGradient(rep(left, gradientSteps), rep(right, gradientSteps), seq(y, y+exonHeight/2, len=gradientSteps), rgb(1,1,1,0.6))
			drawVerticalGradient(rep(left, gradientSteps), rep(right, gradientSteps), seq(y, y-exonHeight/2, len=gradientSteps), rgb(1,1,1,0.6))
			# add exon label
			text((left+right)/2, y, title, cex=0.9)
		}
	}
}

drawCircos <- function(fusion, fusions, ideograms, minConfidenceForCircosPlot) {
	# initialize with empty circos plot
	circos.clear()
	circos.initializeWithIdeogram(cytoband=ideograms, plotType=NULL)
	# use gene names as labels or <contig>:<position> for intergenic breakpoints
	geneLabels <- data.frame(
		contig=c(fusions[fusion,"contig1"], fusions[fusion,"contig2"]),
		start=c(fusions[fusion,"breakpoint1"], fusions[fusion,"breakpoint2"])
	)
	geneLabels$end <- geneLabels$start + 1
	geneLabels$gene <- c(fusions[fusion,"gene1"], fusions[fusion,"gene2"])
	geneLabels$gene <- ifelse(grepl(",", geneLabels$gene), paste0(geneLabels$contig, ":", geneLabels$start), geneLabels$gene)
	# draw gene labels
	circos.genomicLabels(geneLabels, labels.column=4, side="outside")
	# draw chromosome labels in connector plot
	for (contig in unique(ideograms$contig)) {
		set.current.cell(track.index=2, sector.index=contig) # draw in gene label connector track (track.index=2)
		circos.text(CELL_META$xcenter, CELL_META$ycenter, contig, cex=0.75)
	}
	# draw ideograms
	circos.genomicIdeogram(cytoband=ideograms)
	# draw arcs
	confidenceRank <- c(low=0, medium=1, high=2)
	for (i in c(setdiff(1:nrow(fusions), fusion), fusion)) { # draw fusion of interest last, such that its arc is on top
		f <- fusions[i,]
		if (confidenceRank[f$confidence] >= confidenceRank[minConfidenceForCircosPlot])
			circos.link(
				f$contig1, f$breakpoint1,
				f$contig2, f$breakpoint2,
				lwd=2, col=ifelse(
					(f$gene1 != fusions[fusion,"gene1"] | f$gene2 != fusions[fusion,"gene2"]) && (f$gene1 != fusions[fusion,"gene2"] | f$gene2 != fusions[fusion,"gene1"]),
					rgb(1,0.7,0.7), # pale arcs for other fusions
					rgb(1,0,0) # solid arc for fusion of interest
				)
			)
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
	# sort coding exons last, such that they are drawn over the border of non-coding exons
	candidateExons <- unique(candidateExons[order(candidateExons$start, -rank(candidateExons$type)),])
	return(candidateExons)
}

for (fusion in 1:nrow(fusions)) {

	message(paste0("Drawing fusion #", fusion, ": ", fusions[fusion,"gene1"], ":", fusions[fusion,"gene2"]))

	exons1 <- findExons(exons, fusions[fusion,"gene1"], fusions[fusion,"direction1"], fusions[fusion,"contig1"], fusions[fusion,"breakpoint1"])
	if (nrow(exons1) == 0) {
		par(mfrow=c(1,1))
		plot(0, 0, type="l", xaxt="n", yaxt="n", xlab="", ylab="")
		text(0, 0, paste0("Error: exon coordinates of ", fusions[fusion,"gene1"], " not found in\n", exonsFile))
		next
	}
	exons2 <- findExons(exons, fusions[fusion,"gene2"], fusions[fusion,"direction2"], fusions[fusion,"contig2"], fusions[fusion,"breakpoint2"])
	if (nrow(exons2) == 0) {
		par(mfrow=c(1,1))
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
		coverageNormalization <- ifelse(
			squishIntrons, # => ignore intronic coverage
			max(
				apply(exons1[,c("start", "end")], 1, function(e) {max(as.numeric(coverage1[IRanges(e["start"],e["end"])]))}),
				apply(exons2[,c("start", "end")], 1, function(e) {max(as.numeric(coverage2[IRanges(e["start"],e["end"])]))})
			),
			max(max(coverage1), max(coverage2))
		)
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

	# layout: fusion on top, circos plot on bottom left, statistics on bottom right
	layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))
	par(mar=c(0, 0, 0, 0))
	plot(0, 0, type="l", xlim=c(-0.12, 1.12), ylim=c(0.4, 1), bty="n", xaxt="n", yaxt="n")

	# vertical coordinates of layers
	yIdeograms <- 0.94
	yGeneNames <- 0.61
	yBreakpointLabels <- 0.86
	yCoverage <- 0.72
	yExons <- 0.67
	yFusion <- 0.5
	yTranscript <- 0.43

	# draw ideograms
	if (printIdeograms) {
		drawIdeogram("left", min(exons1$left), max(exons1$right), yIdeograms, ideograms, fusions[fusion,"contig1"], fusions[fusion,"breakpoint1"])
		drawIdeogram("right", gene2Offset, gene2Offset+max(exons2$right), yIdeograms, ideograms, fusions[fusion,"contig2"], fusions[fusion,"breakpoint2"])
	}

	# draw gene & transcript names
	text(max(exons1$right)/2, yGeneNames, fusions[fusion,"gene1"], font=2)
	if (!grepl(",", head(exons1$transcript,1)))
		text(max(exons1$right)/2, yGeneNames-0.03, head(exons1$transcript,1), cex=0.9)
	text(gene2Offset+max(exons2$right)/2, yGeneNames, fusions[fusion,"gene2"], font=2)
	if (!grepl(",", head(exons2$transcript,1)))
		text(gene2Offset+max(exons2$right)/2, yGeneNames-0.03, head(exons2$transcript,1), cex=0.9)

	# label breakpoints
	text(breakpoint1+0.01, yBreakpointLabels, paste0("breakpoint\n", fusions[fusion,"contig1"], ":", fusions[fusion,"breakpoint1"]), adj=c(1,0.5))
	text(gene2Offset+breakpoint2-0.01, yBreakpointLabels, paste0("breakpoint\n", fusions[fusion,"contig2"], ":", fusions[fusion,"breakpoint2"]), adj=c(0,0.5))

	# draw coverage axis
	lines(c(-0.02, -0.01, -0.01, -0.02), c(yCoverage, yCoverage, yCoverage+0.1, yCoverage+0.1))
	text(-0.025, yCoverage, "0", adj=c(1,0.5), cex=0.9)
	text(-0.025, yCoverage+0.1, coverageNormalization, adj=c(1,0.5), cex=0.9)
	text(-0.05, yCoverage+0.04, "Coverage", srt=90, cex=0.9)
	rect(min(exons1$left), yCoverage, max(exons1$right), yCoverage+0.1, col="#eeeeee", border=NA)
	rect(gene2Offset+min(exons2$left), yCoverage, gene2Offset+max(exons2$right), yCoverage+0.1, col="#eeeeee", border=NA)

	# plot coverage 1
	if (squishIntrons) {
		for (exon in 1:nrow(exons1))
			drawCoverage(exons1[exon,"left"], exons1[exon,"right"], yCoverage, coverage1, exons1[exon,"start"], exons1[exon,"end"], color1)
	} else {
		drawCoverage(min(exons1$left), max(exons1$right), yCoverage, coverage1, min(exons1$start), max(exons1$end), color1)
	}

	# plot coverage 2
	if (squishIntrons) {
		for (exon in 1:nrow(exons2))
			drawCoverage(gene2Offset+exons2[exon,"left"], gene2Offset+exons2[exon,"right"], yCoverage, coverage2, exons2[exon,"start"], exons2[exon,"end"], color2)
	} else {
		drawCoverage(gene2Offset+min(exons2$left), gene2Offset+max(exons2$right), yCoverage, coverage2, min(exons2$start), max(exons2$end), color2)
	}

	# plot gene 1
	drawStrand(0, max(exons1$right), yExons, darkColor1, exons1[1,"strand"])
	for (exon in 1:nrow(exons1))
		drawExon(exons1[exon,"left"], exons1[exon,"right"], yExons, color1, exons1[exon,"exonNumber"], exons1[exon,"type"])

	# plot gene 2
	drawStrand(gene2Offset, gene2Offset+max(exons2$right), yExons, col=darkColor2, exons2[1,"strand"])
	for (exon in 1:nrow(exons2))
		drawExon(gene2Offset+exons2[exon,"left"], gene2Offset+exons2[exon,"right"], yExons, color2, exons2[exon,"exonNumber"], exons2[exon,"type"])

	# plot gene1 of fusion
	if (fusions[fusion,"direction1"] == "downstream") {
		# plot strands
		drawStrand(fusionOffset1, fusionOffset1+breakpoint1, yFusion, col=darkColor1, exons1[1,"strand"])
		# plot exons
		for (exon in 1:nrow(exons1))
			if (exons1[exon,"start"] <= fusions[fusion,"breakpoint1"])
				drawExon(fusionOffset1+exons1[exon,"left"], fusionOffset1+min(breakpoint1, exons1[exon,"right"]), yFusion, color1, exons1[exon,"exonNumber"], exons1[exon,"type"])
		# plot trajectories
		lines(c(0, 0, fusionOffset1), c(yExons+0.03, yExons-0.05, yFusion+0.03), col="red", lty=2)
		lines(c(breakpoint1, breakpoint1, fusionOffset1+breakpoint1), c(yBreakpointLabels-0.03, yExons-0.05, yFusion+0.03), col="red", lty=2)
	} else if (fusions[fusion,"direction1"] == "upstream") {
		# plot strands
		drawStrand(fusionOffset1, fusionOffset2, yFusion, col=darkColor1, chartr("+-", "-+", exons1[1,"strand"]))
		# plot exons
		for (exon in 1:nrow(exons1))
			if (exons1[exon,"end"]+1 >= fusions[fusion,"breakpoint1"])
				drawExon(fusionOffset1+max(exons1$right)-exons1[exon,"right"], min(fusionOffset2, fusionOffset1+max(exons1$right)-exons1[exon,"left"]), yFusion, color1, exons1[exon,"exonNumber"], exons1[exon,"type"])
		# plot trajectories
		lines(c(max(exons1$right), max(exons1$right), fusionOffset1), c(yExons+0.03, yExons-0.05, yFusion+0.03), col="red", lty=2)
		lines(c(breakpoint1, breakpoint1, fusionOffset1+max(exons1$right)-breakpoint1), c(yBreakpointLabels-0.03, yExons-0.05, yFusion+0.03), col="red", lty=2)
	}
	
	# plot gene2 of fusion
	if (fusions[fusion,"direction2"] == "downstream") {
		# plot strands
		drawStrand(fusionOffset2, fusionOffset2+breakpoint2, yFusion, col=darkColor2, chartr("+-", "-+", exons2[1,"strand"]))
		# plot exons
		for (exon in 1:nrow(exons2))
			if (exons2[exon,"start"] <= fusions[fusion,"breakpoint2"])
				drawExon(max(fusionOffset2, fusionOffset2+breakpoint2-exons2[exon,"right"]), fusionOffset2+breakpoint2-exons2[exon,"left"], yFusion, color2, exons2[exon,"exonNumber"], exons2[exon,"type"])
		# plot trajectories
		lines(c(gene2Offset, gene2Offset, fusionOffset2+breakpoint2), c(yExons+0.03, yExons-0.05, yFusion+0.03), col="red", lty=2)
		lines(c(gene2Offset+breakpoint2, gene2Offset+breakpoint2, fusionOffset2), c(yBreakpointLabels-0.03, yExons-0.05, yFusion+0.03), col="red", lty=2)
	} else if (fusions[fusion,"direction2"] == "upstream") {
		# plot strands
		drawStrand(fusionOffset2, fusionOffset2+max(exons2$right)-breakpoint2, yFusion, col=darkColor2, exons2[1,"strand"])
		# plot exons
		for (exon in 1:nrow(exons2))
			if (exons2[exon,"end"]+1 >= fusions[fusion,"breakpoint2"])
				drawExon(max(fusionOffset2, fusionOffset2+exons2[exon,"left"]-breakpoint2), fusionOffset2+exons2[exon,"right"]-breakpoint2, yFusion, color2, exons2[exon,"exonNumber"], exons2[exon,"type"])
		# plot trajectories
		lines(c(gene2Offset+max(exons2$right), gene2Offset+max(exons2$right), fusionOffset2+max(exons2$right)-breakpoint2), c(yExons+0.03, yExons-0.05, yFusion+0.03), col="red", lty=2)
		lines(c(gene2Offset+breakpoint2, gene2Offset+breakpoint2, fusionOffset2), c(yBreakpointLabels-0.03, yExons-0.05, yFusion+0.03), col="red", lty=2)
	}
	
	if (printFusionTranscript && fusions[fusion,"fusion_transcript"] != ".") {
		# print fusion transcript colored by gene of origin
		fusion_transcript1 <- gsub("\\|.*", "", fusions[fusion,"fusion_transcript"], perl=T)
		fusion_transcript1 <- substr(fusion_transcript1, max(1, nchar(fusion_transcript1)-30), nchar(fusion_transcript1))
		fusion_transcript2 <- gsub(".*\\|", "", fusions[fusion,"fusion_transcript"], perl=T)
		fusion_transcript2 <- substr(fusion_transcript2, 1, min(nchar(fusion_transcript2), 30))
		text(fusionOffset2, yTranscript, fusion_transcript1, col=darkColor1, adj=c(1,0.5))
		text(fusionOffset2, yTranscript, fusion_transcript1, col=darkColor2, adj=c(0,0.5))
	}

	if (!printCircos) {
		plot(0, 0, type="l", xlim=c(0, 1), ylim=c(0, 1), bty="n", xaxt="n", yaxt="n")
	} else {
		drawCircos(fusion, fusions, ideograms, minConfidenceForCircosPlot)
	}

	plot(0, 0, type="l", xlim=c(0, 1), ylim=c(0, 1), bty="n", xaxt="n", yaxt="n")
	if (printStats) {
		# print statistics about supporting alignments
		text(0, 0.5, paste("Split reads in", fusions[fusion,"gene1"], "=", fusions[fusion,"split_reads1"]), adj=c(0,0.5))
		text(0, 0.45, paste("Split reads in", fusions[fusion,"gene2"], "=", fusions[fusion,"split_reads2"]), adj=c(0,0.5))
		text(0, 0.40, paste("Discordant mates =", fusions[fusion,"discordant_mates"]), adj=c(0,0.5))
	}

}

devNull <- dev.off()
