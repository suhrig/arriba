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
parseFileParameter <- function(parameter, args, mandatory=FALSE) {
	fileName <- sub(paste0("^--", parameter, "="), "", args[tail(which(grepl(paste0("^--", parameter, "="), args)), 1)], perl=T)
	if (length(fileName) > 0) {
		if (file.access(fileName) == -1)
			stop(paste("Cannot read file:", fileName))
	} else {
		if (mandatory)
			stop(paste0("Missing mandatory argument: --", parameter))
		fileName = ""
	}
	return(fileName)
}

if (any(grepl("^--help", args)) || length(args) == 0)
	stop("Usage: draw_fusions.R --annotation=annotation.gtf --fusions=fusions.tsv --output=output.pdf [--alignments=Aligned.out.bam] [--cytobands=cytobands.tsv] [--minConfidenceForCircosPlot=medium] [--proteinDomains=protein_domains.gff3] [--squishIntrons=TRUE] [--printExonLabels=TRUE] [--pdfWidth=11.692] [--pdfHeight=8.267] [--color1=#e5a5a5] [--color2=#a7c4e5] [--mergeDomainsOverlappingBy=0.9]")
exonsFile <- parseFileParameter("annotation", args, T)
fusionsFile <- parseFileParameter("fusions", args, T)
outputFile <- parseStringParameter("output", args)
if (outputFile == "")
	stop("Missing mandatory argument: --output")
alignmentsFile <- parseFileParameter("alignments", args)
cytobandsFile <- parseFileParameter("cytobands", args)
if (cytobandsFile == "")
	warning("Missing parameter '--cytobands'. No ideograms and circos plots will be drawn.")
minConfidenceForCircosPlot <- parseStringParameter("minConfidenceForCircosPlot", args, "medium")
if (!(minConfidenceForCircosPlot %in% c("low", "medium", "high")))
	stop("Invalid argument to --minConfidenceForCircosPlot")
proteinDomainsFile <- parseFileParameter("proteinDomains", args)
squishIntrons <- parseBooleanParameter("squishIntrons", args, T)
printExonLabels <- parseBooleanParameter("printExonLabels", args, T)
pdfWidth <- as.numeric(parseStringParameter("pdfWidth", args, "11.692"))
pdfHeight <- as.numeric(parseStringParameter("pdfHeight", args, "8.267"))
color1 <- parseStringParameter("color1", args, "#e5a5a5")
color2 <- parseStringParameter("color2", args, "#a7c4e5")
mergeDomainsOverlappingBy <- as.numeric(parseStringParameter("mergeDomainsOverlappingBy", args, 0.9))

# check if required packages are installed
if (!suppressPackageStartupMessages(require(GenomicRanges)))
	warning("Package 'GenomicRanges' is not installed. No circos plots will be drawn.")
if (!suppressPackageStartupMessages(require(circlize)))
	warning("Package 'circlize' is not installed. No circos plots will be drawn.")
if (alignmentsFile != "")
	if (!suppressPackageStartupMessages(require(GenomicAlignments)))
		stop("Package 'GenomicAlignments' must be installed when '--alignments' is used")

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

pdf(outputFile, onefile=T, width=pdfWidth, height=pdfHeight, title=fusionsFile)
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
	sub("^chr", "", sub("^chrM", "MT", contig), perl=T)
}

# read cytoband annotation
cytobands <- NULL
if (cytobandsFile != "") {
	cytobands <- read.table(cytobandsFile, header=T, colClasses=c("character", "numeric", "numeric", "character", "character"))
	cytobands <- cytobands[order(cytobands$contig, cytobands$start, cytobands$end),]
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

# read protein domain annotation
proteinDomains <- NULL
if (proteinDomainsFile != "") {
	message("Loading protein domains")
	proteinDomains <- read.table(proteinDomainsFile, header=F, sep="\t", comment.char="", quote="", stringsAsFactors=F)[,c(1,4,5,7,9)]
	colnames(proteinDomains) <- c("contig", "start", "end", "strand", "attributes")
	proteinDomains$color <- sub(";.*", "", sub(".*color=", "", proteinDomains$attributes, perl=T), perl=T)
	proteinDomains$proteinDomainName <- sapply(sub(";.*", "", sub(".*Name=", "", proteinDomains$attributes, perl=T), perl=T), URLdecode)
	proteinDomains$proteinDomainID <- sub(";.*", "", sub(".*protein_domain_id=", "", proteinDomains$attributes, perl=T), perl=T)
	proteinDomains <- proteinDomains[,colnames(proteinDomains) != "attributes"]
}

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

drawIdeogram <- function(adjust, left, right, y, cytobands, contig, breakpoint) {
	# define design of ideogram
	bandColors <- c(gneg="#ffffff", gpos25="#bbbbbb", gpos50="#888888", gpos75="#444444", gpos100="#000000", acen="#ec4f4f", stalk="#0000ff")
	cytobands$color <- bandColors[cytobands$giemsa]
	arcSteps <- 30 # defines roundness of arc
	curlyBraceHeight <- 0.03
	ideogramHeight <- 0.04
	ideogramWidth <- 0.4
	# extract bands of given contig
	bands <- cytobands[cytobands$contig==contig,]
	if (nrow(bands) == 0)
		stop(paste("Giemsa bands of contig", contig, "not found"))
	# scale width of ideogram to fit inside given region
	bands$left <- bands$start / max(cytobands$end) * ideogramWidth
	bands$right <- bands$end / max(cytobands$end) * ideogramWidth
	# left/right-align cytobands
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
			rect(left, y+exonHeight, right, y+exonHeight/2-0.0015, col=color, border=NA)
			rect(left, y-exonHeight, right, y-exonHeight/2+0.0015, col=color, border=NA)
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

drawCircos <- function(fusion, fusions, cytobands, minConfidenceForCircosPlot) {
	# initialize with empty circos plot
	circos.clear()
	circos.initializeWithIdeogram(cytoband=cytobands, plotType=NULL)
	# use gene names as labels or <contig>:<position> for intergenic breakpoints
	geneLabels <- data.frame(
		contig=c(fusions[fusion,"contig1"], fusions[fusion,"contig2"]),
		start=c(fusions[fusion,"breakpoint1"], fusions[fusion,"breakpoint2"])
	)
	geneLabels$end <- geneLabels$start + 1
	geneLabels$gene <- c(fusions[fusion,"gene1"], fusions[fusion,"gene2"])
	geneLabels$gene <- ifelse(grepl(",", geneLabels$gene), paste0(geneLabels$contig, ":", geneLabels$start), geneLabels$gene)
	# draw gene labels
	circos.genomicLabels(geneLabels, labels.column=4, side="outside", cex=1)
	# draw chromosome labels in connector plot
	for (contig in unique(cytobands$contig)) {
		set.current.cell(track.index=2, sector.index=contig) # draw in gene label connector track (track.index=2)
		circos.text(CELL_META$xcenter, CELL_META$ycenter, contig, cex=0.85)
	}
	# draw ideograms
	circos.genomicIdeogram(cytoband=cytobands)
	# draw arcs
	confidenceRank <- c(low=0, medium=1, high=2)
	for (i in c(setdiff(1:nrow(fusions), fusion), fusion)) { # draw fusion of interest last, such that its arc is on top
		f <- fusions[i,]
		if (confidenceRank[f$confidence] >= confidenceRank[minConfidenceForCircosPlot] || i==fusion)
			circos.link(
				f$contig1, f$breakpoint1,
				f$contig2, f$breakpoint2,
				lwd=2, col=ifelse(i==fusion, rgb(1,0,0), rgb(1,0.7,0.7))
			)
	}
}

drawProteinDomains <- function(fusion, exons1, exons2, proteinDomains, color1, color2, mergeDomainsOverlappingBy) {

	exonHeight <- 0.2
	exonsY <- 0.5
	geneNamesY <- exonsY - exonHeight/2 - 0.05

	# find coding exons
	codingExons1 <- exons1[exons1$type == "CDS",]
	codingExons2 <- exons2[exons2$type == "CDS",]

	# cut off coding regions beyond breakpoint
	if (fusion$direction1 == "upstream") {
		codingExons1 <- codingExons1[codingExons1$end >= fusion$breakpoint1,]
		codingExons1$start <- ifelse(codingExons1$start < fusion$breakpoint1, fusion$breakpoint1, codingExons1$start)
	} else {
		codingExons1 <- codingExons1[codingExons1$start <= fusion$breakpoint1,]
		codingExons1$end <- ifelse(codingExons1$end > fusion$breakpoint1, fusion$breakpoint1, codingExons1$end)
	}
	if (fusion$direction2 == "upstream") {
		codingExons2 <- codingExons2[codingExons2$end >= fusion$breakpoint2,]
		codingExons2$start <- ifelse(codingExons2$start < fusion$breakpoint2, fusion$breakpoint2, codingExons2$start)
	} else {
		codingExons2 <- codingExons2[codingExons2$start <= fusion$breakpoint2,]
		codingExons2$end <- ifelse(codingExons2$end > fusion$breakpoint2, fusion$breakpoint2, codingExons2$end)
	}

	# find overlapping domains
	exonsGRanges1 <- GRanges(codingExons1$contig, IRanges(codingExons1$start, codingExons1$end), strand=codingExons1$strand)
	exonsGRanges2 <- GRanges(codingExons2$contig, IRanges(codingExons2$start, codingExons2$end), strand=codingExons2$strand)
	domainsGRanges <- GRanges(proteinDomains$contig, IRanges(proteinDomains$start, proteinDomains$end), strand=proteinDomains$strand)
	domainsGRanges$proteinDomainName <- proteinDomains$proteinDomainName
	domainsGRanges$proteinDomainID <- proteinDomains$proteinDomainID
	domainsGRanges$color <- proteinDomains$color
	domainsGRanges <- domainsGRanges[suppressWarnings(unique(queryHits(findOverlaps(domainsGRanges, union(exonsGRanges1, exonsGRanges2)))))]

	# group overlapping domains by domain ID
	domainsGRangesList <- GRangesList(lapply(unique(domainsGRanges$proteinDomainID), function(x) { domainsGRanges[domainsGRanges$proteinDomainID == x] }))

	# trim protein domains to exon boundaries
	trimDomains <- function(domainsGRangesList, exonsGRanges) {
		do.call(
			"rbind",
			lapply(
				domainsGRangesList,
				function(x) {
					intersected <- as.data.frame(reduce(suppressWarnings(intersect(x, exonsGRanges))))
					if (nrow(intersected) > 0) {
						intersected$proteinDomainName <- head(x$proteinDomainName, 1)
						intersected$proteinDomainID <- head(x$proteinDomainID, 1)
						intersected$color <- head(x$color, 1)
					} else {
						intersected$proteinDomainName <- character()
						intersected$proteinDomainID <- character()
						intersected$color <- character()
					}
					return(intersected)
				}
			)
		)
	}
	retainedDomains1 <- trimDomains(domainsGRangesList, exonsGRanges1)
	retainedDomains2 <- trimDomains(domainsGRangesList, exonsGRanges2)

	# calculate length of coding exons
	codingExons1$length <- codingExons1$end - codingExons1$start + 1
	codingExons2$length <- codingExons2$end - codingExons2$start + 1

	# abort, if there are no coding regions
	if (sum(exons1$type == "CDS") + sum(exons2$type == "CDS") == 0) {
		text(0.5, 0.5, "Genes are not protein-coding.")
		return(NULL)
	}
	codingLength1 <- sum(codingExons1$length)
	codingLength2 <- sum(codingExons2$length)
	if (codingLength1 + codingLength2 == 0) {
		text(0.5, 0.5, "No coding regions retained in fusion transcript.")
		return(NULL)
	}
	antisenseTranscription1 <- sub("/.*", "", fusion$strand1) != sub(".*/", "", fusion$strand1)
	antisenseTranscription2 <- sub("/.*", "", fusion$strand2) != sub(".*/", "", fusion$strand2)
	if ((codingLength1 == 0 || antisenseTranscription1) && (codingLength2 == 0 || antisenseTranscription2)) {
		text(0.5, 0.5, "No coding regions due to antisense transcription.")
		return(NULL)
	}

	# remove introns from protein domains
	removeIntronsFromProteinDomains <- function(codingExons, retainedDomains) {
		if (nrow(codingExons) == 0) return(NULL)
		cumulativeIntronLength <- 0
		previousExonEnd <- 0
		for (exon in 1:nrow(codingExons)) {
			if (codingExons[exon,"start"] > previousExonEnd)
				cumulativeIntronLength <- cumulativeIntronLength + codingExons[exon,"start"] - previousExonEnd
			domainsInExon <- which(retainedDomains$start >= codingExons[exon,"start"] & retainedDomains$start <= codingExons[exon,"end"])
			retainedDomains[domainsInExon,"start"] <- retainedDomains[domainsInExon,"start"] - cumulativeIntronLength
			domainsInExon <- which(retainedDomains$end >= codingExons[exon,"start"] & retainedDomains$end <= codingExons[exon,"end"])
			retainedDomains[domainsInExon,"end"] <- retainedDomains[domainsInExon,"end"] - cumulativeIntronLength
			previousExonEnd <- codingExons[exon,"end"]
		}
		# merge adjacent domains
		retainedDomains <- do.call(
			"rbind",
			lapply(
				unique(retainedDomains$proteinDomainID),
				function(x) {
					domain <- retainedDomains[retainedDomains$proteinDomainID == x,]
					merged <- reduce(GRanges(domain$seqnames, IRanges(domain$start, domain$end), strand=domain$strand))
					merged$proteinDomainName <- head(domain$proteinDomainName, 1)
					merged$proteinDomainID <- head(domain$proteinDomainID, 1)
					merged$color <- head(domain$color, 1)
					return(as.data.frame(merged))
				}
			)
		)
		return(retainedDomains)
	}
	retainedDomains1 <- removeIntronsFromProteinDomains(codingExons1, retainedDomains1)
	retainedDomains2 <- removeIntronsFromProteinDomains(codingExons2, retainedDomains2)

	# abort, if no domains are retained
	if (is.null(retainedDomains1) && is.null(retainedDomains2)) {
		text(0.5, 0.5, "No protein domains retained in fusion.")
		return(NULL)
	}

	# merge domains with similar coordinates
	mergeSimilarDomains <- function(domains, mergeDomainsOverlappingBy) {
		if (is.null(domains)) return(domains)
		merged <- domains[F,] # create empty data frame
		domains <- domains[order(domains$end - domains$start, decreasing=T),] # start with bigger domains => bigger domains are retained
		for (domain in rownames(domains)) {
			if (!any((abs(merged$start - domains[domain,"start"]) + abs(merged$end - domains[domain,"end"])) / (domains[domain,"end"] - domains[domain,"start"]) <= 1-mergeDomainsOverlappingBy))
				merged <- rbind(merged, domains[domain,])
		}
		return(merged)
	}
	retainedDomains1 <- mergeSimilarDomains(retainedDomains1, mergeDomainsOverlappingBy)
	retainedDomains2 <- mergeSimilarDomains(retainedDomains2, mergeDomainsOverlappingBy)

	# normalize length to 1
	codingExons1$length <- codingExons1$length / (codingLength1 + codingLength2)
	codingExons2$length <- codingExons2$length / (codingLength1 + codingLength2)
	retainedDomains1$start <- retainedDomains1$start / (codingLength1 + codingLength2)
	retainedDomains1$end <- retainedDomains1$end / (codingLength1 + codingLength2)
	retainedDomains2$start <- retainedDomains2$start / (codingLength1 + codingLength2)
	retainedDomains2$end <- retainedDomains2$end / (codingLength1 + codingLength2)

	# draw coding regions
	rect(0, exonsY-exonHeight/2, sum(codingExons1$length), exonsY+exonHeight/2, col=color1, border=NA)
	rect(sum(codingExons1$length), exonsY-exonHeight/2, sum(codingExons1$length) + sum(codingExons2$length), exonsY+exonHeight/2, col=color2, border=NA)

	# indicate exon boundaries as dotted lines
	exonBoundaries <- cumsum(c(codingExons1$length, codingExons2$length))
	if (length(exonBoundaries) > 1) {
		exonBoundaries <- exonBoundaries[1:(length(exonBoundaries)-1)]
		for (exonBoundary in exonBoundaries)
			lines(c(exonBoundary, exonBoundary), c(exonsY-exonHeight, exonsY+exonHeight), col="white", lty=3)
	}

	# find overlapping domains
	# nest if one is contained in another
	# stack if they overlap partially
	nestDomains <- function(domains) {
		if (length(unlist(domains)) == 0) return(domains)
		domains <- domains[order(domains$end - domains$start, decreasing=T),]
		rownames(domains) <- 1:nrow(domains)
		# find nested domains and make tree structure
		domains$parent <- 0
		for (domain in rownames(domains))
			domains[domains$start >= domains[domain,"start"] & domains$end <= domains[domain,"end"] & rownames(domains) != domain,"parent"] <- domain
		# find partially overlapping domains
		maxOverlappingDomains <- max(coverage(IRanges(domains$start*10e6, domains$end*10e6)))
		padding <- 1 / maxOverlappingDomains * 0.4
		domains$y <- 0
		domains$height <- 0
		adjustPositionAndHeight <- function(parentDomain, y, height, padding, e) {
			for (domain in which(e$domains$parent == parentDomain)) {
				overlappingDomains <- which((e$domains$start >= e$domains[domain,"start"] & e$domains$start <= e$domains[domain,"end"] |
				                             e$domains$end   >= e$domains[domain,"start"] & e$domains$end   <= e$domains[domain,"end"]) &
				                             e$domains$parent == parentDomain)
				e$domains[domain,"height"] <- height/length(overlappingDomains) - padding * (length(overlappingDomains)-1) / length(overlappingDomains)
				e$domains[domain,"y"] <- y + (which(domain==overlappingDomains)-1) * (e$domains[domain,"height"] + padding)
				adjustPositionAndHeight(domain, e$domains[domain,"y"]+padding, e$domains[domain,"height"]-2*padding, padding, e)
			}
		}
		adjustPositionAndHeight(0, 0, 1, padding, environment())
		domains <- domains[order(domains$height, decreasing=T),] # draw nested domains last
		return(domains)
	}
	retainedDomains1 <- nestDomains(retainedDomains1)
	retainedDomains2 <- nestDomains(retainedDomains2)
	retainedDomains1$y <- exonsY - exonHeight/2 + 0.025 + (exonHeight-2*0.025) * retainedDomains1$y
	retainedDomains2$y <- exonsY - exonHeight/2 + 0.025 + (exonHeight-2*0.025) * retainedDomains2$y
	retainedDomains1$height <- retainedDomains1$height * (exonHeight-2*0.025)
	retainedDomains2$height <- retainedDomains2$height * (exonHeight-2*0.025)

	# draw domains
	drawProteinDomainRect <- function(left, bottom, right, top, color) {
		rect(left, bottom, right, top, col=color, border=getDarkColor(color))
		# draw gradients for 3D effect
		gradientSteps <- 20
		drawVerticalGradient(rep(left, gradientSteps), rep(right, gradientSteps), seq(top, bottom, len=gradientSteps), rgb(1,1,1,0.7))
		drawVerticalGradient(rep(left, gradientSteps), rep(right, gradientSteps), seq(bottom, bottom+(top-bottom)*0.4, len=gradientSteps), rgb(0,0,0,0.1))
	}
	if (length(unlist(retainedDomains1)) > 0)
		for (domain in 1:nrow(retainedDomains1))
			drawProteinDomainRect(retainedDomains1[domain,"start"], retainedDomains1[domain,"y"], retainedDomains1[domain,"end"], retainedDomains1[domain,"y"]+retainedDomains1[domain,"height"], retainedDomains1[domain,"color"])
	if (length(unlist(retainedDomains2)) > 0)
		for (domain in 1:nrow(retainedDomains2))
			drawProteinDomainRect(sum(codingExons1$length)+retainedDomains2[domain,"start"], retainedDomains2[domain,"y"], sum(codingExons1$length)+retainedDomains2[domain,"end"], retainedDomains2[domain,"y"]+retainedDomains2[domain,"height"], retainedDomains2[domain,"color"])

	# draw gene names, if there are coding exons
	if (codingLength1 > 0)
		text(sum(codingExons1$length)/2, geneNamesY, fusion$gene1, font=2)
	if (codingLength2 > 0)
		text(sum(codingExons1$length)+sum(codingExons2$length)/2, geneNamesY, fusion$gene2, font=2)

	# calculate how many non-adjacent unique domains there are
	# we need this info to know where to place labels vertically
	countUniqueDomains <- function(domains) {
		uniqueDomains <- 0
		if (length(unlist(domains)) > 0) {
			uniqueDomains <- 1
			if (nrow(domains) > 1) {
				previousDomain <- domains[1,"proteinDomainID"]
				for (domain in 2:nrow(domains)) {
					if (previousDomain != domains[domain,"proteinDomainID"])
						uniqueDomains <- uniqueDomains + 1
					previousDomain <- domains[domain,"proteinDomainID"]
				}
			}
		}
		return(uniqueDomains)
	}
	if (length(unlist(retainedDomains1)) > 0)
		retainedDomains1 <- retainedDomains1[order(retainedDomains1$start),]
	uniqueDomains1 <- countUniqueDomains(retainedDomains1)
	if (length(unlist(retainedDomains2)) > 0)
		retainedDomains2 <- retainedDomains2[order(retainedDomains2$end, decreasing=T),]
	uniqueDomains2 <- countUniqueDomains(retainedDomains2)

	# draw title of plot
	titleY <- exonsY + exonHeight/2 + (uniqueDomains1 + 1) * 0.05
	text(0.5, titleY, "RETAINED PROTEIN DOMAINS", font=2)

	# draw domain labels for gene1
	if (length(unlist(retainedDomains1)) > 0) {
		previousConnectorX <- -1
		previousLabelX <- -1
		labelY <- exonsY + exonHeight/2 + uniqueDomains1 * 0.05
		for (domain in 1:nrow(retainedDomains1)) {
			# if possible avoid overlapping lines of labels
			connectorX <- min(retainedDomains1[domain,"start"] + 0.01, (retainedDomains1[domain,"start"] + retainedDomains1[domain,"end"])/2)
			if (connectorX - previousConnectorX < 0.01 && retainedDomains1[domain,"end"] > previousConnectorX + 0.01)
				connectorX <- previousConnectorX + 0.01
			labelX <- max(connectorX, previousLabelX) + 0.02
			# use a signle label for adjacent domains of same type
			adjacentDomainsOfSameType <- domain + 1 <= nrow(retainedDomains1) && retainedDomains1[domain+1,"proteinDomainID"] == retainedDomains1[domain,"proteinDomainID"]
			if (adjacentDomainsOfSameType) {
				labelX <- retainedDomains1[domain+1,"start"] + 0.015
			} else {
				text(labelX, labelY, retainedDomains1[domain,"proteinDomainName"], adj=c(0,0.5), col=getDarkColor(retainedDomains1[domain,"color"]))
			}
			lines(c(labelX-0.005, connectorX, connectorX), c(labelY, labelY, retainedDomains1[domain,"y"]+retainedDomains1[domain,"height"]), col=getDarkColor(retainedDomains1[domain,"color"]))
			if (!adjacentDomainsOfSameType)
				labelY <- labelY - 0.05
			previousConnectorX <- connectorX
			previousLabelX <- labelX
		}
	}

	# draw domain labels for gene2
	if (length(unlist(retainedDomains2)) > 0) {
		previousConnectorX <- 100
		previousLabelX <- 100
		labelY <- exonsY - exonHeight/2 - (uniqueDomains2+1) * 0.05
		for (domain in 1:nrow(retainedDomains2)) {
			# if possible avoid overlapping connector lines of labels
			connectorX <- sum(codingExons1$length) + max(retainedDomains2[domain,"end"] - 0.01, (retainedDomains2[domain,"start"] + retainedDomains2[domain,"end"])/2)
			if (previousConnectorX - connectorX < 0.01 && sum(codingExons1$length) + retainedDomains2[domain,"start"] < previousConnectorX - 0.01)
				connectorX <- previousConnectorX - 0.01
			labelX <- min(connectorX, previousLabelX) - 0.02
			# use a signle label for adjacent domains of same type
			adjacentDomainsOfSameType <- domain + 1 <= nrow(retainedDomains2) && retainedDomains2[domain+1,"proteinDomainID"] == retainedDomains2[domain,"proteinDomainID"]
			if (adjacentDomainsOfSameType) {
				labelX <- sum(codingExons1$length) + retainedDomains2[domain+1,"end"] - 0.015
			} else {
				text(labelX, labelY, retainedDomains2[domain,"proteinDomainName"], adj=c(1,0.5), col=getDarkColor(retainedDomains2[domain,"color"]))
			}
			lines(c(labelX+0.005, connectorX, connectorX), c(labelY, labelY, retainedDomains2[domain,"y"]), col=getDarkColor(retainedDomains2[domain,"color"]))
			if (!adjacentDomainsOfSameType)
				labelY <- labelY + 0.05
			previousConnectorX <- connectorX
			previousLabelX <- labelX
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

	# layout: fusion on top, circos plot on bottom left, protein domains on bottom center, statistics on bottom right
	layout(matrix(c(1,1,1,2,3,4), 2, 3, byrow=TRUE), widths=c(0.9, 1.2, 0.9))
	par(mar=c(0, 0, 0, 0))
	plot(0, 0, type="l", xlim=c(-0.12, 1.12), ylim=c(0.4, 1.1), bty="n", xaxt="n", yaxt="n")

	# vertical coordinates of layers
	yIdeograms <- 0.94
	yGeneNames <- 0.61
	yBreakpointLabels <- 0.86
	yCoverage <- 0.72
	yExons <- 0.67
	yFusion <- 0.5
	yTranscript <- 0.43

	# draw ideograms
	if (!is.null(cytobands)) {
		drawIdeogram("left", min(exons1$left), max(exons1$right), yIdeograms, cytobands, fusions[fusion,"contig1"], fusions[fusion,"breakpoint1"])
		drawIdeogram("right", gene2Offset, gene2Offset+max(exons2$right), yIdeograms, cytobands, fusions[fusion,"contig2"], fusions[fusion,"breakpoint2"])
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
	if (alignmentsFile != "") {
		lines(c(-0.02, -0.01, -0.01, -0.02), c(yCoverage, yCoverage, yCoverage+0.1, yCoverage+0.1))
		text(-0.025, yCoverage, "0", adj=c(1,0.5), cex=0.9)
		text(-0.025, yCoverage+0.1, coverageNormalization, adj=c(1,0.5), cex=0.9)
		text(-0.05, yCoverage+0.04, "Coverage", srt=90, cex=0.9)
		rect(min(exons1$left), yCoverage, max(exons1$right), yCoverage+0.1, col="#eeeeee", border=NA)
		rect(gene2Offset+min(exons2$left), yCoverage, gene2Offset+max(exons2$right), yCoverage+0.1, col="#eeeeee", border=NA)
	}

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
	
	if (fusions[fusion,"fusion_transcript"] != ".") {
		# print fusion transcript colored by gene of origin
		fusion_transcript1 <- gsub("\\|.*", "", fusions[fusion,"fusion_transcript"], perl=T)
		fusion_transcript1 <- substr(fusion_transcript1, max(1, nchar(fusion_transcript1)-30), nchar(fusion_transcript1))
		fusion_transcript2 <- gsub(".*\\|", "", fusions[fusion,"fusion_transcript"], perl=T)
		fusion_transcript2 <- substr(fusion_transcript2, 1, min(nchar(fusion_transcript2), 30))
		# check for non-template bases
		non_template_bases <- gsub(".*\\|([^|]*)\\|.*", "\\1", fusions[fusion,"fusion_transcript"])
		if (non_template_bases == fusions[fusion,"fusion_transcript"]) # no non-template bases found
			non_template_bases <- ""
		# divide non-template bases half-and-half for centered alignment
		non_template_bases1 <- substr(non_template_bases, 1, floor(nchar(non_template_bases)/2))
		non_template_bases2 <- substr(non_template_bases, ceiling(nchar(non_template_bases)/2), nchar(non_template_bases))
		# transcript 1
		text(fusionOffset2, yTranscript, bquote(.(fusion_transcript1) * phantom(.(non_template_bases1))), col=darkColor1, adj=c(1,0.5))
		# transcript 2
		text(fusionOffset2, yTranscript, bquote(phantom(.(non_template_bases2)) * .(fusion_transcript2)), col=darkColor2, adj=c(0,0.5))
		# non-template bases
		text(fusionOffset2, yTranscript, non_template_bases1, adj=c(1,0.5))
		text(fusionOffset2, yTranscript, non_template_bases2, adj=c(0,0.5))
	}

	if (is.null(cytobands) || !("circlize" %in% names(sessionInfo()$otherPkgs)) || !("GenomicRanges" %in% names(sessionInfo()$otherPkgs))) {
		plot(0, 0, type="l", xlim=c(0, 1), ylim=c(0, 1), bty="n", xaxt="n", yaxt="n")
	} else {
		par(mar=c(0, 4, 0, 0))
		drawCircos(fusion, fusions, cytobands, minConfidenceForCircosPlot)
		par(mar=c(0, 0, 0, 0))
	}

	plot(0, 0, type="l", xlim=c(-0.1, 1.1), ylim=c(0, 1), bty="n", xaxt="n", yaxt="n")
	if (!is.null(proteinDomains))
		drawProteinDomains(fusions[fusion,], exons1, exons2, proteinDomains, color1, color2, mergeDomainsOverlappingBy)

	# print statistics about supporting alignments
	plot(0, 0, type="l", xlim=c(0, 1), ylim=c(0, 1), bty="n", xaxt="n", yaxt="n")
	text(0, 0.575, "SUPPORTING READ COUNT", font=2, adj=c(0,0.5))
	text(0, 0.525, paste("Split reads in", fusions[fusion,"gene1"], "=", fusions[fusion,"split_reads1"]), adj=c(0,0.5))
	text(0, 0.475, paste("Split reads in", fusions[fusion,"gene2"], "=", fusions[fusion,"split_reads2"]), adj=c(0,0.5))
	text(0, 0.425, paste("Discordant mates =", fusions[fusion,"discordant_mates"]), adj=c(0,0.5))

}

devNull <- dev.off()
