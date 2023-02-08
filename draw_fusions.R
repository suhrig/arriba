#!/usr/bin/env Rscript

# print warnings as they happen instead of collecting them for after a loop ends
options(warn=1)

# define valid parameters
parameters <- list(
	fusions=list("fusionsFile", "file", "fusions.tsv", T),
	annotation=list("exonsFile", "file", "annotation.gtf", T),
	output=list("outputFile", "string", "output.pdf", T),
	alignments=list("alignmentsFile", "file", "Aligned.sortedByCoord.out.bam"),
	cytobands=list("cytobandsFile", "file", "cytobands.tsv"),
	minConfidenceForCircosPlot=list("minConfidenceForCircosPlot", "string", "medium"),
	proteinDomains=list("proteinDomainsFile", "file", "protein_domains.gff3"),
	sampleName=list("sampleName", "string", ""),
	squishIntrons=list("squishIntrons", "bool", T),
	printExonLabels=list("printExonLabels", "bool", T),
	render3dEffect=list("render3dEffect", "bool", T),
	pdfWidth=list("pdfWidth", "numeric", 11.692),
	pdfHeight=list("pdfHeight", "numeric", 8.267),
	color1=list("color1", "string", "#e5a5a5"),
	color2=list("color2", "string", "#a7c4e5"),
	mergeDomainsOverlappingBy=list("mergeDomainsOverlappingBy", "numeric", 0.9),
	optimizeDomainColors=list("optimizeDomainColors", "bool", F),
	fontSize=list("fontSize", "numeric", 1),
	fontFamily=list("fontFamily", "string", "Helvetica"),
	showIntergenicVicinity=list("showIntergenicVicinity", "string", "0"),
	transcriptSelection=list("transcriptSelection", "string", "provided"),
	fixedScale=list("fixedScale", "numeric", 0),
	coverageRange=list("coverageRange", "string", "0")
)

# print help if necessary
args <- commandArgs(trailingOnly=T)
if (any(grepl("^--help", args)) || length(args) == 0) {
	usage <- "Usage: draw_fusions.R"
	for (parameter in names(parameters)) {
		usage <- paste0(usage, " ")
		if (length(parameters[[parameter]]) <= 3 || !parameters[[parameter]][[4]])
			usage <- paste0(usage, "[")
		usage <- paste0(usage, "--", parameter, "=", parameters[[parameter]][[3]])
		if (length(parameters[[parameter]]) <= 3 || !parameters[[parameter]][[4]])
			usage <- paste0(usage, "]")
	}
	message(usage)
	quit("no", ifelse(length(args) == 0, 1, 0))
}

# make sure mandatory arguments are present
for (parameter in names(parameters))
	if (length(parameters[[parameter]]) > 3 && parameters[[parameter]][[4]])
		if (!any(grepl(paste0("^--", parameter, "="), args), perl=T))
			stop(paste0("Missing mandatory argument: --", parameter))

# set default values
for (parameter in names(parameters))
	assign(parameters[[parameter]][[1]], ifelse(parameters[[parameter]][[2]] == "file", "", parameters[[parameter]][[3]]))

# parse command-line parameters
for (arg in args) {
	argName <- sub("=.*", "", sub("^--", "", arg, perl=T), perl=T)
	argValue <- sub("^[^=]*=", "", arg, perl=T)
	if (!(argName %in% names(parameters)) || !grepl("^--", arg, perl=T))
		stop(paste("Unknown parameter:", arg))
	if (parameters[[argName]][[2]] == "bool") {
		if (argValue %in% c("TRUE", "T", "FALSE", "F")) {
			assign(parameters[[argName]][[1]], as.logical(argValue))
		} else {
			stop(paste0("Invalid argument to --", argName))
		}
	} else if (parameters[[argName]][[2]] == "string") {
		assign(parameters[[argName]][[1]], argValue)
	} else if (parameters[[argName]][[2]] == "numeric") {
		if (is.na(suppressWarnings(as.numeric(argValue))))
			stop(paste0("Invalid argument to --", argName))
		assign(parameters[[argName]][[1]], as.numeric(argValue))
	} else if (parameters[[argName]][[2]] == "file") {
		if (file.access(argValue) == -1)
			stop(paste("Cannot read file:", argValue))
		assign(parameters[[argName]][[1]], argValue)
	}
}

# validate values of parameters
if (cytobandsFile == "")
	warning("Missing parameter '--cytobands'. No ideograms and circos plots will be drawn.")
if (!(minConfidenceForCircosPlot %in% c("none", "low", "medium", "high")))
	stop("Invalid argument to --minConfidenceForCircosPlot")
showIntergenicVicinity <- as.list(unlist(strsplit(showIntergenicVicinity, ",", fixed=T)))
if (!(length(showIntergenicVicinity) %in% c(1,4)))
	stop(paste0("Invalid argument to --showIntergenicVicinity"))
showIntergenicVicinity <- lapply(showIntergenicVicinity, function(x) {
	if (x == "closestGene") {
		return("exon")
	} else if (x == "closestProteinCodingGene") {
		return("CDS")
	} else if (is.na(suppressWarnings(as.numeric(x))) || as.numeric(x) < 0) {
		stop(paste0("Invalid argument to --showIntergenicVicinity"))
	} else {
		return(as.numeric(x))
	}
})
if (length(showIntergenicVicinity) == 1)
	showIntergenicVicinity <- rep(showIntergenicVicinity, 4)
if (squishIntrons)
	if (any(!is.numeric(unlist(showIntergenicVicinity))) || any(showIntergenicVicinity > 0))
		stop("--squishIntrons must be disabled, when --showIntergenicVicinity is > 0")
if (!(transcriptSelection %in% c("coverage", "provided", "canonical")))
	stop("Invalid argument to --transcriptSelection")
if (fixedScale < 0)
	stop("Invalid argument to --fixedScale")
if (!(fontFamily %in% names(pdfFonts())))
	stop(paste0("Unknown font: ", fontFamily, ". Available fonts: ", paste(names(pdfFonts()), collapse=", ")))
coverageRange <- suppressWarnings(as.numeric(unlist(strsplit(coverageRange, ",", fixed=T))))
if (!(length(coverageRange) %in% 1:2) || any(is.na(coverageRange)) || any(coverageRange < 0))
	stop("Invalid argument to --coverageRange")

# check if required packages are installed
if (!suppressPackageStartupMessages(require(GenomicRanges)))
	warning("Package 'GenomicRanges' is not installed. No protein domains and circos plots will be drawn.")
if (!suppressPackageStartupMessages(require(circlize)))
	warning("Package 'circlize' is not installed. No circos plots will be drawn.")
if (alignmentsFile != "")
	if (!suppressPackageStartupMessages(require(GenomicAlignments)))
		stop("Package 'GenomicAlignments' must be installed when '--alignments' is used")

# define colors
changeColorBrightness <- function(color, delta) {
	rgb(
		min(255,max(0,col2rgb(color)["red",]+delta)),
		min(255,max(0,col2rgb(color)["green",]+delta)),
		min(255,max(0,col2rgb(color)["blue",]+delta)),
		maxColorValue=255
	)
}
getDarkColor <- function(color) { changeColorBrightness(color, -100) }
getBrightColor <- function(color) { changeColorBrightness(color, +190) }
darkColor1 <- getDarkColor(color1)
darkColor2 <- getDarkColor(color2)
circosColors <- c(translocation="#000000", duplication="#00bb00", deletion="#ff0000", inversion="#0000ff")

# convenience functions to add/remove "chr" prefix
addChr <- function(contig) {
	ifelse(contig == "MT", "chrM", paste0("chr", contig))
}
removeChr <- function(contig) {
	sub("^chr", "", sub("^chrM", "MT", contig, perl=T), perl=T)
}

# convenience function to check if a value is between two others
between <- function(value, start, end) {
	value >= start & value <= end
}

# read fusions
fusions <- read.table(fusionsFile, stringsAsFactors=F, sep="\t", header=T, comment.char="", quote="")
if (colnames(fusions)[1] == "X.gene1") { # Arriba output
	colnames(fusions)[colnames(fusions) %in% c("X.gene1", "strand1.gene.fusion.", "strand2.gene.fusion.")] <- c("gene1", "strand1", "strand2")
	fusions$display_contig1 <- sub(":[^:]*$", "", fusions$breakpoint1, perl=T)
	fusions$display_contig2 <- sub(":[^:]*$", "", fusions$breakpoint2, perl=T)
	fusions$contig1 <- removeChr(fusions$display_contig1)
	fusions$contig2 <- removeChr(fusions$display_contig2)
	fusions$breakpoint1 <- as.numeric(sub(".*:", "", fusions$breakpoint1, perl=T))
	fusions$breakpoint2 <- as.numeric(sub(".*:", "", fusions$breakpoint2, perl=T))
	fusions$split_reads1 <- fusions$split_reads1
	fusions$split_reads2 <- fusions$split_reads2
	fusions$type <- sub(".*(translocation|duplication|deletion|inversion).*", "\\1", fusions$type, perl=T)
	fusions$fusion_transcript <- gsub("[()^$]", "", fusions$fusion_transcript)
} else if (colnames(fusions)[1] == "X.FusionName") { # STAR-Fusion
	fusions$gene1 <- sub("\\^.*", "", fusions$LeftGene, perl=T)
	fusions$gene2 <- sub("\\^.*", "", fusions$RightGene, perl=T)
	fusions$strand1 <- sub(".*:(.)$", "\\1/\\1", fusions$LeftBreakpoint, perl=T)
	fusions$strand2 <- sub(".*:(.)$", "\\1/\\1", fusions$RightBreakpoint, perl=T)
	fusions$display_contig1 <- sub(":[^:]*:[^:]*$", "", fusions$LeftBreakpoint, perl=T)
	fusions$display_contig2 <- sub(":[^:]*:[^:]*$", "", fusions$RightBreakpoint, perl=T)
	fusions$contig1 <- removeChr(fusions$display_contig1)
	fusions$contig2 <- removeChr(fusions$display_contig2)
	fusions$breakpoint1 <- as.numeric(sub(".*:([^:]*):[^:]*$", "\\1", fusions$LeftBreakpoint, perl=T))
	fusions$breakpoint2 <- as.numeric(sub(".*:([^:]*):[^:]*$", "\\1", fusions$RightBreakpoint, perl=T))
	fusions$direction1 <- ifelse(grepl(":\\+$", fusions$LeftBreakpoint, perl=T), "downstream", "upstream")
	fusions$direction2 <- ifelse(grepl(":\\+$", fusions$RightBreakpoint, perl=T), "upstream", "downstream")
	fusions$gene_id1 <- sub(".*\\^", "", fusions$LeftGene, perl=T)
	fusions$gene_id2 <- sub(".*\\^", "", fusions$RightGene, perl=T)
	fusions$transcript_id1 <- ifelse(rep(!("CDS_LEFT_ID" %in% colnames(fusions)), nrow(fusions)), ".", fusions$CDS_LEFT_ID)
	fusions$transcript_id2 <- ifelse(rep(!("CDS_RIGHT_ID" %in% colnames(fusions)), nrow(fusions)), ".", fusions$CDS_RIGHT_ID)
	fusions$fusion_transcript <- ifelse(rep(!("FUSION_CDS" %in% colnames(fusions)), nrow(fusions)), ".", toupper(sub("([a-z]*)", "\\1|", fusions$FUSION_CDS, perl=T)))
	fusions$reading_frame <- ifelse(rep(!("PROT_FUSION_TYPE" %in% colnames(fusions)), nrow(fusions)), ".", ifelse(fusions$PROT_FUSION_TYPE == "INFRAME", "in-frame", ifelse(fusions$PROT_FUSION_TYPE == "FRAMESHIFT", "out-of-frame", ".")))
	fusions$split_reads <- fusions$JunctionReadCount
	fusions$discordant_mates <- fusions$SpanningFragCount
	fusions$site1 <- rep("exon", nrow(fusions))
	fusions$site2 <- rep("exon", nrow(fusions))
	fusions$confidence <- rep("high", nrow(fusions))
	fusions$type <- ifelse(fusions$contig1 != fusions$contig2, "translocation", ifelse(fusions$direction1 == fusions$direction2, "inversion", ifelse((fusions$direction1 == "downstream") == (fusions$breakpoint1 < fusions$breakpoint2), "deletion", "duplication")))
} else {
	stop("Unrecognized fusion file format")
}

pdf(outputFile, onefile=T, width=pdfWidth, height=pdfHeight, title=ifelse(sampleName != "", sampleName, fusionsFile))
par(family=fontFamily)

if (nrow(fusions) == 0) {
	plot(0, 0, type="l", xaxt="n", yaxt="n", xlab="", ylab="")
	text(0, 0, "empty input file")
	warning("empty input file")
	dev.off()
	quit("no")
}

# read cytoband annotation
cytobands <- NULL
if (cytobandsFile != "") {
	cytobands <- read.table(cytobandsFile, header=T, colClasses=c("character", "numeric", "numeric", "character", "character"))
	cytobands <- cytobands[order(cytobands$contig, cytobands$start, cytobands$end),]
}

# read exon annotation
message("Loading annotation")
exons <- scan(exonsFile, what=list(contig="",src="",type="",start=0,end=0,score="",strand="",frame="",attributes=""), sep="\t", comment.char="#", quote='"', multi.line=F)
attr(exons, "row.names") <- .set_row_names(length(exons[[1]]))
class(exons) <- "data.frame"
exons <- exons[exons$type %in% c("exon","CDS"),c("contig","type","start","end","strand","attributes")]
exons$contig <- removeChr(exons$contig)
parseGtfAttribute <- function(attribute, gtf) {
	parsed <- sub(paste0(".*", attribute, "[ =]([^;]+).*"), "\\1", gtf$attributes, perl=T)
	failedToParse <- parsed == gtf$attributes
	if (any(failedToParse)) {
		warning(paste0("Failed to parse '", attribute, "' attribute of ", sum(failedToParse), " record(s)."))
		parsed <- ifelse(failedToParse, "", parsed)
	}
	return(parsed)
}
exons$geneID <- parseGtfAttribute("gene_id", exons)
exons$geneName <- parseGtfAttribute("gene_name", exons)
exons$geneName <- ifelse(exons$geneName == "", exons$geneID, exons$geneName)
exons$transcript <- parseGtfAttribute("transcript_id", exons)
exons$exonNumber <- ifelse(rep(printExonLabels, nrow(exons)), parseGtfAttribute("exon_number", exons), "")

# read protein domain annotation
proteinDomains <- NULL
if (proteinDomainsFile != "") {
	message("Loading protein domains")
	proteinDomains <- scan(proteinDomainsFile, what=list(contig="",src="",type="",start=0,end=0,score="",strand="",frame="",attributes=""), sep="\t", comment.char="", quote="", multi.line=F)
	attr(proteinDomains, "row.names") <- .set_row_names(length(proteinDomains[[1]]))
	class(proteinDomains) <- "data.frame"
	proteinDomains$color <- parseGtfAttribute("color", proteinDomains)
	proteinDomains$proteinDomainName <- sapply(parseGtfAttribute("Name", proteinDomains), URLdecode)
	proteinDomains$proteinDomainID <- parseGtfAttribute("protein_domain_id", proteinDomains)
}

# insert dummy annotations for intergenic breakpoints
if (any(fusions$site1 == "intergenic" | fusions$site2 == "intergenic")) {
	intergenicBreakpoints <- rbind(
		setNames(fusions[fusions$site1 == "intergenic",c("gene1", "strand1", "contig1", "breakpoint1")], c("gene", "strand", "contig", "breakpoint")),
		setNames(fusions[fusions$site2 == "intergenic",c("gene2", "strand2", "contig2", "breakpoint2")], c("gene", "strand", "contig", "breakpoint"))
	)
	exons <- rbind(exons, data.frame(
		contig=intergenicBreakpoints$contig,
		type="intergenic",
		start=sapply(intergenicBreakpoints$breakpoint-1000, max, 1),
		end=intergenicBreakpoints$breakpoint+1000,
		strand=".",
		attributes="",
		geneName=intergenicBreakpoints$gene,
		geneID=paste0(intergenicBreakpoints$contig, ":", intergenicBreakpoints$breakpoint),
		transcript=paste0(intergenicBreakpoints$contig, ":", intergenicBreakpoints$breakpoint),
		exonNumber="intergenic"
	))
	fusions[fusions$site1 == "intergenic","gene_id1"] <- paste0(fusions[fusions$site1 == "intergenic","contig1"], ":", fusions[fusions$site1 == "intergenic","breakpoint1"])
	fusions[fusions$site2 == "intergenic","gene_id2"] <- paste0(fusions[fusions$site2 == "intergenic","contig2"], ":", fusions[fusions$site2 == "intergenic","breakpoint2"])
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
if (!render3dEffect) # nullify function, if no 3D effect should be drawn
	drawVerticalGradient <- function(left, right, y, color, selection=NULL) { }

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
	bandColors <- setNames(rgb(100:0, 100:0, 100:0, maxColorValue=100), paste0("gpos", 0:100))
	bandColors <- c(bandColors, gneg="#ffffff", acen="#ec4f4f", stalk="#0000ff")
	cytobands$color <- bandColors[cytobands$giemsa]
	arcSteps <- 30 # defines roundness of arc
	curlyBraceHeight <- 0.03
	ideogramHeight <- 0.04
	ideogramWidth <- 0.4
	# extract bands of given contig
	bands <- cytobands[cytobands$contig==contig,]
	if (nrow(bands) == 0) {
		warning(paste("Ideogram of contig", contig, "cannot be drawn, because no Giemsa staining information is available."))
		return(NULL)
	}
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
	text((max(bands$right)+min(bands$left))/2, y+0.07, paste("chromosome", contig), font=2, cex=fontSize, adj=c(0.5,0))
	# draw name of band
	bandName <- bands[which(between(breakpoint, bands$start, bands$end)), "name"]
	text(tip, y+0.03, bandName, cex=fontSize, adj=c(0.5,0))
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
	# if there is no centromere, make an artificial one with length zero
	if (is.null(centromereStart) || is.null(centromereEnd)) {
		centromereStart <- bands[1,"right"]
		centromereEnd <- bands[1,"right"]
	}
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
	maxResolution <- 5000 # max number of data points to draw coverage
	# draw coverage as bars
	if (!is.null(coverage)) {
		coverageData <- as.numeric(coverage[IRanges(sapply(start, max, min(start(coverage))), sapply(end, min, max(end(coverage))))])
		# downsample to maxResolution, if there are too many data points
		coverageData <- aggregate(coverageData, by=list(round(1:length(coverageData) * (right-left) * maxResolution/length(coverageData))), mean)$x
		polygon(c(left, seq(left, right, length.out=length(coverageData)), right), c(y, y+coverageData*0.1, y), col=color, border=NA)
	}
}

drawStrand <- function(left, right, y, color, strand) {
	if (strand %in% c("+", "-")) {
		# draw strand
		lines(c(left+0.001, right-0.001), c(y, y), col=color, lwd=2)
		lines(c(left+0.001, right-0.001), c(y, y), col=rgb(1,1,1,0.1), lwd=1)
		# indicate orientation
		if (right - left > 0.01)
			for (i in seq(left+0.005, right-0.005, by=sign(right-left-2*0.005)*0.01)) {
				arrows(i, y, i+0.001*ifelse(strand=="+", 1, -1), y, col=color, length=0.05, lwd=2, angle=60)
				arrows(i, y, i+0.001*ifelse(strand=="+", 1, -1), y, col=rgb(1,1,1,0.1), length=0.05, lwd=1, angle=60)
			}
	}
}

drawExon <- function(left, right, y, color, title, type) {
	gradientSteps <- 10 # defines smoothness of gradient
	exonHeight <- 0.03
	if (type == "CDS") {
		# draw coding regions as thicker bars
		rect(left, y+exonHeight, right, y+exonHeight/2-0.001, col=color, border=NA)
		rect(left, y-exonHeight, right, y-exonHeight/2+0.001, col=color, border=NA)
		# draw border
		lines(c(left, left, right, right), c(y+exonHeight/2, y+exonHeight, y+exonHeight, y+exonHeight/2), col=getDarkColor(color), lend=2)
		lines(c(left, left, right, right), c(y-exonHeight/2, y-exonHeight, y-exonHeight, y-exonHeight/2), col=getDarkColor(color), lend=2)
		# draw gradients for 3D effect
		drawVerticalGradient(rep(left, gradientSteps), rep(right, gradientSteps), seq(y+0.03, y+0.015, len=gradientSteps), rgb(0,0,0,0.2))
		drawVerticalGradient(rep(left, gradientSteps), rep(right, gradientSteps), seq(y-0.03, y-0.015, len=gradientSteps), rgb(0,0,0,0.3))
	} else if (type == "exon") {
		rect(left, y+exonHeight/2, right, y-exonHeight/2, col=color, border=getDarkColor(color))
		# draw gradients for 3D effect
		drawVerticalGradient(rep(left, gradientSteps), rep(right, gradientSteps), seq(y, y+exonHeight/2, len=gradientSteps), rgb(1,1,1,0.6))
		drawVerticalGradient(rep(left, gradientSteps), rep(right, gradientSteps), seq(y, y-exonHeight/2, len=gradientSteps), rgb(1,1,1,0.6))
		# add exon label
		text((left+right)/2, y, title, cex=0.9*fontSize)
	}
}

drawCircos <- function(fusion, fusions, cytobands, minConfidenceForCircosPlot, circosColors) {
	# check if Giemsa staining information is available
	for (contig in unlist(fusions[fusion,c("contig1", "contig2")])) {
		if (!any(cytobands$contig==contig)) {
			warning(paste0("Circos plot cannot be drawn, because no Giemsa staining information is available for contig ", contig, "."))
			# draw empty plots as placeholder
			plot(0, 0, type="l", xlim=c(0, 1), ylim=c(0, 1), bty="n", xaxt="n", yaxt="n", xlab="", ylab="")
			plot(0, 0, type="l", xlim=c(0, 1), ylim=c(0, 1), bty="n", xaxt="n", yaxt="n", xlab="", ylab="")
			return(NULL)
		}
	}
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
	geneLabels$gene <- ifelse(c(fusions[fusion,"site1"], fusions[fusion,"site2"]) == "intergenic", paste0(c(fusions[fusion,"display_contig1"], fusions[fusion,"display_contig2"]), ":", geneLabels$start), geneLabels$gene)
	# draw gene labels
	circos.genomicLabels(geneLabels, labels.column=4, side="outside", cex=fontSize, labels_height=0.27)
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
		if (any(cytobands$contig==f$contig1) && any(cytobands$contig==f$contig2)) # ignore viral contigs, because we have no cytoband information for them
			if (minConfidenceForCircosPlot != "none" && confidenceRank[f$confidence] >= confidenceRank[minConfidenceForCircosPlot] || i==fusion)
				circos.link(
					f$contig1, f$breakpoint1,
					f$contig2, f$breakpoint2,
					lwd=2, col=ifelse(i==fusion, circosColors[f$type], getBrightColor(circosColors[f$type]))
				)
	}
	# draw legend
	plot(0, 0, type="l", xlim=c(0, 1), ylim=c(0, 1), bty="n", xaxt="n", yaxt="n", ylab="", xlab="")
	legend(x="top", legend=names(circosColors), col=sapply(circosColors, getBrightColor), lwd=3, ncol=2, box.lty=0)
}

drawProteinDomains <- function(fusion, exons1, exons2, proteinDomains, color1, color2, mergeDomainsOverlappingBy, optimizeDomainColors) {

	exonHeight <- 0.2
	exonsY <- 0.5
	geneNamesY <- exonsY - exonHeight/2 - 0.05

	# find coding exons
	codingExons1 <- exons1[exons1$type == "CDS" & fusion$site1 != "intergenic",]
	codingExons2 <- exons2[exons2$type == "CDS" & fusion$site2 != "intergenic",]

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
	if ((codingLength1 == 0 || grepl("\\.$", fusion$strand1)) && (codingLength2 == 0 || grepl("\\.$", fusion$strand2))) {
		text(0.5, 0.5, "Failed to determine retained protein domains due to lack of strand information.")
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
			domainsInExon <- which(between(retainedDomains$start, codingExons[exon,"start"], codingExons[exon,"end"]))
			retainedDomains[domainsInExon,"start"] <- retainedDomains[domainsInExon,"start"] - cumulativeIntronLength
			domainsInExon <- which(between(retainedDomains$end, codingExons[exon,"start"], codingExons[exon,"end"]))
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
			if (!any((abs(merged$start - domains[domain,"start"]) + abs(merged$end - domains[domain,"end"])) / (domains[domain,"end"] - domains[domain,"start"] + 1) <= 1-mergeDomainsOverlappingBy))
				merged <- rbind(merged, domains[domain,])
		}
		return(merged)
	}
	retainedDomains1 <- mergeSimilarDomains(retainedDomains1, mergeDomainsOverlappingBy)
	retainedDomains2 <- mergeSimilarDomains(retainedDomains2, mergeDomainsOverlappingBy)

	# if desired, reassign colors to protein domains to maximize contrast
	if (optimizeDomainColors) {
		uniqueDomains <- unique(c(retainedDomains1$proteinDomainID, retainedDomains2$proteinDomainID))
		# make rainbow of pretty pastell colors
		colors <- rainbow(length(uniqueDomains))
		colors <- apply(col2rgb(colors), 2, function(x) { 0.3 + y/255 * 0.7 }) # make pastell colors
		colors <- apply(colors, 2, function(x) {rgb(x["red"], x["green"], x["blue"])}) # convert back to rgb
		# reassign colors
		names(colors) <- uniqueDomains
		retainedDomains1$color <- colors[retainedDomains1$proteinDomainID]
		retainedDomains2$color <- colors[retainedDomains2$proteinDomainID]
	}

	# reverse exons and protein domains, if on the reverse strand
	if (any(codingExons1$strand == "-")) {
		codingExons1$length <- rev(codingExons1$length)
		temp <- retainedDomains1$end
		retainedDomains1$end <- codingLength1 - retainedDomains1$start
		retainedDomains1$start <- codingLength1 - temp
	}
	if (any(codingExons2$strand == "-")) {
		codingExons2$length <- rev(codingExons2$length)
		temp <- retainedDomains2$end
		retainedDomains2$end <- codingLength2 - retainedDomains2$start
		retainedDomains2$start <- codingLength2 - temp
	}

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
		maxOverlappingDomains <- max(1, as.integer(coverage(IRanges(domains$start*10e6, domains$end*10e6))))
		padding <- 1 / maxOverlappingDomains * 0.4
		domains$y <- 0
		domains$height <- 0
		adjustPositionAndHeight <- function(parentDomain, y, height, padding, e) {
			for (domain in which(e$domains$parent == parentDomain)) {
				overlappingDomains <- which((between(e$domains$start, e$domains[domain,"start"], e$domains[domain,"end"]) |
				                             between(e$domains$end  , e$domains[domain,"start"], e$domains[domain,"end"])) &
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
		text(sum(codingExons1$length)/2, geneNamesY, fusion$gene1, font=2, cex=fontSize)
	if (codingLength2 > 0)
		text(sum(codingExons1$length)+sum(codingExons2$length)/2, geneNamesY, fusion$gene2, font=2, cex=fontSize)

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
	titleY <- exonsY + exonHeight/2 + (uniqueDomains1 + 2) * 0.05
	text(0.5, titleY+0.01, "RETAINED PROTEIN DOMAINS", adj=c(0.5, 0), font=2, cex=fontSize)
	text(0.5, titleY, ifelse(fusion$reading_frame %in% c("in-frame", "out-of-frame"), paste(fusion$reading_frame, "fusion"), ifelse(fusion$reading_frame == "stop-codon", "stop codon before fusion junction", "reading frame unclear")), adj=c(0.5, 1), cex=fontSize)

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
				text(labelX, labelY, retainedDomains1[domain,"proteinDomainName"], adj=c(0,0.5), col=getDarkColor(retainedDomains1[domain,"color"]), cex=fontSize)
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
				text(labelX, labelY, retainedDomains2[domain,"proteinDomainName"], adj=c(1,0.5), col=getDarkColor(retainedDomains2[domain,"color"]), cex=fontSize)
			}
			lines(c(labelX+0.005, connectorX, connectorX), c(labelY, labelY, retainedDomains2[domain,"y"]), col=getDarkColor(retainedDomains2[domain,"color"]))
			if (!adjacentDomainsOfSameType)
				labelY <- labelY + 0.05
			previousConnectorX <- connectorX
			previousLabelX <- labelX
		}
	}

}

findExons <- function(exons, contig, geneID, direction, breakpoint, coverage, transcriptId, transcriptSelection) {
	# use the provided transcript if desired
	if (transcriptSelection == "provided" && transcriptId != "." && transcriptId != "") {
		candidateExons <- exons[exons$transcript == transcriptId,]
		if (nrow(candidateExons) == 0) {
			warning(paste0("Unknown transcript given in fusions file (", transcriptId, "), selecting a different one"))
		} else {
			return(candidateExons)
		}
	}

	if (transcriptSelection == "canonical") {
		candidateExons <- exons[exons$geneID == geneID & exons$contig == contig,]
	} else {
		# look for exon with breakpoint as splice site
		transcripts <- exons[exons$geneID == geneID & exons$contig == contig & exons$type == "exon" & (direction == "downstream" & abs(exons$end - breakpoint) <= 2 | direction == "upstream" & abs(exons$start - breakpoint) <= 2),"transcript"]
		candidateExons <- exons[exons$transcript %in% transcripts,]
		# if none was found, use all exons of the gene closest to the breakpoint
		if (nrow(candidateExons) == 0)
			candidateExons <- exons[exons$geneID == geneID & exons$contig == contig,]
		# if we have coverage information, use the transcript with the highest coverage if there are multiple hits
		if (!is.null(coverage)) {
			highestCoverage <- -1
			transcriptWithHighestCoverage <- NULL
			lengthOfTranscriptWithHighestCoverage <- 0
			for (transcript in unique(candidateExons$transcript)) {
				exonsOfTranscript <- candidateExons[candidateExons$transcript==transcript,]
				exonsOfTranscript$start <- sapply(exonsOfTranscript$start, max, min(start(coverage)))
				exonsOfTranscript$end <- sapply(exonsOfTranscript$end, min, max(end(coverage)))
				lengthOfTranscript <- sum(exonsOfTranscript$end - exonsOfTranscript$start + 1)
				coverageSum <- sum(as.numeric(coverage[IRanges(exonsOfTranscript$start, exonsOfTranscript$end)]))
				# we prefer shorter transcripts over longer ones, because otherwise there is a bias towards transcripts with long UTRs
				# => a longer transcript must have substantially higher coverage to replace a shorter one
				substantialDifference <- (1 - min(lengthOfTranscript, lengthOfTranscriptWithHighestCoverage) / max(lengthOfTranscript, lengthOfTranscriptWithHighestCoverage)) / 10
				if (lengthOfTranscript > lengthOfTranscriptWithHighestCoverage && coverageSum * (1-substantialDifference) > highestCoverage ||
				    lengthOfTranscript < lengthOfTranscriptWithHighestCoverage && coverageSum > highestCoverage * (1-substantialDifference)) {
					highestCoverage <- coverageSum
					transcriptWithHighestCoverage <- transcript
					lengthOfTranscriptWithHighestCoverage <- lengthOfTranscript
				}
			}
			if (highestCoverage > 0)
				candidateExons <- candidateExons[candidateExons$transcript==transcriptWithHighestCoverage,]
		}
		# if the gene has multiple transcripts, search for transcripts which encompass the breakpoint
		if (length(unique(candidateExons$transcript)) > 1) {
			transcriptStart <- aggregate(candidateExons$start, by=list(candidateExons$transcript), min)
			rownames(transcriptStart) <- transcriptStart[,1]
			transcriptEnd <- aggregate(candidateExons$end, by=list(candidateExons$transcript), max)
			rownames(transcriptEnd) <- transcriptEnd[,1]
			encompassingExons <- between(breakpoint, transcriptStart[candidateExons$transcript,2], transcriptEnd[candidateExons$transcript,2])
			if (any(encompassingExons))
				candidateExons <- candidateExons[encompassingExons,]
		}
	}

	# find the consensus transcript, if there are multiple hits
	if (length(unique(candidateExons$transcript)) > 1) {
		consensusTranscript <-
			ifelse(grepl("appris_principal_1", candidateExons$attributes), 12,
			ifelse(grepl("appris_principal_2", candidateExons$attributes), 11,
			ifelse(grepl("appris_principal_3", candidateExons$attributes), 10,
			ifelse(grepl("appris_principal_4", candidateExons$attributes), 9,
			ifelse(grepl("appris_principal_5", candidateExons$attributes), 8,
			ifelse(grepl("appris_principal", candidateExons$attributes), 7,
			ifelse(grepl("appris_candidate_longest", candidateExons$attributes), 6,
			ifelse(grepl("appris_candidate", candidateExons$attributes), 5,
			ifelse(grepl("appris_alternative_1", candidateExons$attributes), 4,
			ifelse(grepl("appris_alternative_2", candidateExons$attributes), 3,
			ifelse(grepl("appris_alternative", candidateExons$attributes), 2,
			ifelse(grepl("CCDS", candidateExons$attributes), 1,
			0
		))))))))))))
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
	candidateExons <- unique(candidateExons[candidateExons$transcript == head(unique(candidateExons$transcript), 1),])
	return(candidateExons)
}

findClosestGene <- function(exons, contig, breakpoint, extraConditions) {

	# find exons near breakpoint (extraConditions must define what is considered "near")
	closestExons <- exons[exons$contig == contig & extraConditions,] # find closest exon
	closestExons <- exons[exons$contig == contig & exons$geneID %in% closestExons$geneID,] # select all exons of closest gene

	# when more than one gene found with the given name, use the closest one
	if (length(unique(closestExons$geneID)) > 1) { # more than one gene found with the given name => use the closest one
		distanceToBreakpoint <- aggregate(1:nrow(closestExons), by=list(closestExons$geneID), function(x) { min(abs(closestExons[x,"start"]-breakpoint), abs(closestExons[x,"end"]-breakpoint)) })
		closestGene <- head(distanceToBreakpoint[distanceToBreakpoint[,2] == min(distanceToBreakpoint[,2]),1], 1)
		closestExons <- closestExons[closestExons$geneID == closestGene,]
	}

	# when no gene was found, return default values
	if (nrow(closestExons) == 0) {
		return(IRanges(max(1, breakpoint-1000), breakpoint+1000))
	} else {
		return(IRanges(min(closestExons$start), max(closestExons$end)))
	}
}

# main loop starts here
for (fusion in 1:nrow(fusions)) {

	message(paste0("Drawing fusion #", fusion, ": ", fusions[fusion,"gene1"], ":", fusions[fusion,"gene2"]))

	# if showIntergenicVicinity is a number, take it as is
	# if it is a keyword (closestGene/closestProteinCodingGene), determine the range dynamically
	showVicinity <- rep(0, 4)
	if (fusions[fusion,"site1"] == "intergenic") {
		showVicinity[1] <- ifelse(
			is.numeric(showIntergenicVicinity[[1]]),
			showIntergenicVicinity[[1]],
			fusions[fusion,"breakpoint1"] - start(findClosestGene(exons, fusions[fusion,"contig1"], fusions[fusion,"breakpoint1"], exons$end < fusions[fusion,"breakpoint1"] & exons$type == showIntergenicVicinity[[1]]))
		)
		showVicinity[2] <- ifelse(
			is.numeric(showIntergenicVicinity[[2]]),
			showIntergenicVicinity[[2]],
			end(findClosestGene(exons, fusions[fusion,"contig1"], fusions[fusion,"breakpoint1"], exons$start > fusions[fusion,"breakpoint1"] & exons$type == showIntergenicVicinity[[2]])) - fusions[fusion,"breakpoint1"]
		)
	}
	if (fusions[fusion,"site2"] == "intergenic") {
		showVicinity[3] <- ifelse(
			is.numeric(showIntergenicVicinity[[3]]),
			showIntergenicVicinity[[3]],
			fusions[fusion,"breakpoint2"] - start(findClosestGene(exons, fusions[fusion,"contig2"], fusions[fusion,"breakpoint2"], exons$end < fusions[fusion,"breakpoint2"] & exons$type == showIntergenicVicinity[[3]]))
		)
		showVicinity[4] <- ifelse(
			is.numeric(showIntergenicVicinity[[4]]),
			showIntergenicVicinity[[4]],
			end(findClosestGene(exons, fusions[fusion,"contig2"], fusions[fusion,"breakpoint2"], exons$start > fusions[fusion,"breakpoint2"] & exons$type == showIntergenicVicinity[[4]])) - fusions[fusion,"breakpoint2"]
		)
	}

	# compute coverage from alignments file
	coverage1 <- NULL
	coverage2 <- NULL
	if (alignmentsFile != "") {
		# determine range in which we need to compute the coverage
		determineCoverageRegion <- function(exons, geneID, contig, breakpoint, showVicinityLeft, showVicinityRight) {
			closestGene <- findClosestGene(exons, contig, breakpoint, exons$geneID == geneID)
			return(IRanges(min(start(closestGene), breakpoint-showVicinityLeft), max(end(closestGene), breakpoint+showVicinityRight)))
		}
		coverageRegion1 <- determineCoverageRegion(exons, fusions[fusion,"gene_id1"], fusions[fusion,"contig1"], fusions[fusion,"breakpoint1"], showVicinity[1], showVicinity[2])
		coverageRegion2 <- determineCoverageRegion(exons, fusions[fusion,"gene_id2"], fusions[fusion,"contig2"], fusions[fusion,"breakpoint2"], showVicinity[3], showVicinity[4])
		# function which reads alignments from BAM file with & without "chr" prefix
		readCoverage <- function(alignmentsFile, contig, coverageRegion) {
			coverageData <- tryCatch(
				{
					alignments <- readGAlignments(alignmentsFile, param=ScanBamParam(which=GRanges(contig, coverageRegion)))
					coverage(alignments)[[contig]]
				},
				error=function(e) {
					alignments <- readGAlignments(alignmentsFile, param=ScanBamParam(which=GRanges(addChr(contig), coverageRegion)))
					coverage(alignments)[[addChr(contig)]]
				}
			)
			if (exists("alignments")) rm(alignments)
			return(coverageData)
		}
		# get coverage track
		coverage1 <- readCoverage(alignmentsFile, fusions[fusion,"contig1"], coverageRegion1)
		coverage2 <- readCoverage(alignmentsFile, fusions[fusion,"contig2"], coverageRegion2)
		# shrink coverage range to chromosome boundaries to avoid subscript out of bounds errors
		coverageRegion1 <- IRanges(max(start(coverageRegion1), min(start(coverage1))), min(end(coverageRegion1), max(end(coverage1))))
		coverageRegion2 <- IRanges(max(start(coverageRegion2), min(start(coverage2))), min(end(coverageRegion2), max(end(coverage2))))
	}

	# find all exons belonging to the fused genes
	exons1 <- findExons(exons, fusions[fusion,"contig1"], fusions[fusion,"gene_id1"], fusions[fusion,"direction1"], fusions[fusion,"breakpoint1"], coverage1, fusions[fusion,"transcript_id1"], transcriptSelection)
	if (nrow(exons1) == 0) {
		par(mfrow=c(1,1))
		plot(0, 0, type="l", xaxt="n", yaxt="n", xlab="", ylab="")
		text(0, 0, paste0("exon coordinates of ", fusions[fusion,"gene1"], " not found in\n", exonsFile))
		warning(paste("exon coordinates of", fusions[fusion,"gene1"], "not found"))
		next
	}
	exons2 <- findExons(exons, fusions[fusion,"contig2"], fusions[fusion,"gene_id2"], fusions[fusion,"direction2"], fusions[fusion,"breakpoint2"], coverage2, fusions[fusion,"transcript_id2"], transcriptSelection)
	if (nrow(exons2) == 0) {
		par(mfrow=c(1,1))
		plot(0, 0, type="l", xaxt="n", yaxt="n", xlab="", ylab="")
		text(0, 0, paste0("exon coordinates of ", fusions[fusion,"gene2"], " not found in\n", exonsFile))
		warning(paste("exon coordinates of", fusions[fusion,"gene2"], "not found"))
		next
	}

	# in case of intergenic breakpoints, show the vicinity
	if (sum(showVicinity) > 0) {
		if (fusions[fusion,"site1"] == "intergenic") {
			for (geneID in unique(exons[exons$contig == fusions[fusion,"contig1"] & exons$exonNumber != "intergenic" &
			                            (between(exons$end, fusions[fusion,"breakpoint1"]-showVicinity[1], fusions[fusion,"breakpoint1"]+showVicinity[2]) |
			                             between(exons$start, fusions[fusion,"breakpoint1"]-showVicinity[1], fusions[fusion,"breakpoint1"]+showVicinity[2])),"geneID"]))
				exons1 <- rbind(exons1, findExons(exons, fusions[fusion,"contig1"], geneID, fusions[fusion,"direction1"], fusions[fusion,"breakpoint1"], coverage1, fusions[fusion,"transcript_id1"], transcriptSelection))
			# crop genes that are only partially within user-defined vicinity, because the coverage data is incomplete for those
			exons1 <- exons1[exons1$start >= fusions[fusion,"breakpoint1"]-showVicinity[1] & exons1$end <= fusions[fusion,"breakpoint1"]+showVicinity[2] | exons1$exonNumber == "intergenic",]
		}
		if (fusions[fusion,"site2"] == "intergenic") {
			for (geneID in unique(exons[exons$contig == fusions[fusion,"contig2"] & exons$exonNumber != "intergenic" &
			                            (between(exons$end, fusions[fusion,"breakpoint2"]-showVicinity[3], fusions[fusion,"breakpoint2"]+showVicinity[4]) |
			                             between(exons$start, fusions[fusion,"breakpoint2"]-showVicinity[3], fusions[fusion,"breakpoint2"]+showVicinity[4])),"geneID"]))
				exons2 <- rbind(exons2, findExons(exons, fusions[fusion,"contig2"], geneID, fusions[fusion,"direction2"], fusions[fusion,"breakpoint2"], coverage2, fusions[fusion,"transcript_id2"], transcriptSelection))
			# crop genes that are only partially within user-defined vicinity, because the coverage data is incomplete for those
			exons2 <- exons2[exons2$start >= fusions[fusion,"breakpoint2"]-showVicinity[3] & exons2$end <= fusions[fusion,"breakpoint2"]+showVicinity[4] | exons2$exonNumber == "intergenic",]
		}
	}

	# normalize coverage
	if (alignmentsFile != "") {
		coverageNormalization <- function(coverage, coverageRegion, exons) {
			max(1, ifelse(
				squishIntrons, # => ignore intronic coverage
				max(as.numeric(coverage[IRanges(sapply(exons$start,max,min(start(coverage))),sapply(exons$end,min,max(end(coverage))))])),
				round(quantile(coverage[coverageRegion], 0.9999)) # ignore coverage spikes from read-attracting regions
			))
		}
		coverageNormalization1 <- ifelse(head(coverageRange,1) == 0, coverageNormalization(coverage1, coverageRegion1, exons1), head(coverageRange,1))
		coverageNormalization2 <- ifelse(tail(coverageRange,1) == 0, coverageNormalization(coverage2, coverageRegion2, exons2), tail(coverageRange,1))
		if (length(coverageRange) == 1 && coverageRange == 0) { # harmonize scales of gene1 and gene2
			coverageNormalization1 <- max(coverageNormalization1, coverageNormalization2)
			coverageNormalization2 <- max(coverageNormalization1, coverageNormalization2)
		}
		coverage1 <- coverage1/coverageNormalization1
		coverage2 <- coverage2/coverageNormalization2
		coverage1[coverage1 > 1] <- 1
		coverage2[coverage2 > 1] <- 1
	}

	# sort coding exons last, such that they are drawn over the border of non-coding exons
	exons1 <- exons1[order(exons1$start, -rank(exons1$type)),]
	exons2 <- exons2[order(exons2$start, -rank(exons2$type)),]

	# insert dummy exons, if breakpoints are outside the gene (e.g., in UTRs)
	# this avoids plotting artifacts
	breakpoint1 <- fusions[fusion,"breakpoint1"]
	breakpoint2 <- fusions[fusion,"breakpoint2"]
	if (breakpoint1 < min(exons1$start)) {
		exons1 <- rbind(c(exons1[1,"contig"], "dummy", max(1,breakpoint1-1000), max(1,breakpoint1-1000), exons1[1,"strand"], "", "dummy", exons1[1,"geneID"], exons1[1,"transcript"], ""), exons1)
	} else if (breakpoint1 > max(exons1$end)) {
		exons1 <- rbind(exons1, c(exons1[1,"contig"], "dummy", breakpoint1+1000, breakpoint1+1000, exons1[1,"strand"], "", "dummy", exons1[1,"geneID"], exons1[1,"transcript"], ""))
	}
	if (breakpoint2 < min(exons2$start)) {
		exons2 <- rbind(c(exons2[1,"contig"], "dummy", max(1,breakpoint2-1000), max(1,breakpoint2-1000), exons2[1,"strand"], "", "dummy", exons2[1,"geneID"], exons2[1,"transcript"], ""), exons2)
	} else if (breakpoint2 > max(exons2$end)) {
		exons2 <- rbind(exons2, c(exons2[1,"contig"], "dummy", breakpoint2+1000, breakpoint2+1000, exons2[1,"strand"], "", "dummy", exons2[1,"geneID"], exons2[1,"transcript"], ""))
	}
	exons1$start <- as.integer(exons1$start)
	exons1$end <- as.integer(exons1$end)
	exons2$start <- as.integer(exons2$start)
	exons2$end <- as.integer(exons2$end)

	exons1$left <- exons1$start
	exons1$right <- exons1$end
	exons2$left <- exons2$start
	exons2$right <- exons2$end

	squishedIntronSize <- 200
	if (squishIntrons) {
		# hide introns in gene1
		cumulativeIntronLength <- 0
		previousExonEnd <- -squishedIntronSize
		for (exon in 1:nrow(exons1)) {
			if (breakpoint1 > previousExonEnd+1 && breakpoint1 < exons1[exon,"left"])
				breakpoint1 <- (breakpoint1-previousExonEnd) / (exons1[exon,"left"]-previousExonEnd) * squishedIntronSize + previousExonEnd - cumulativeIntronLength
			if (exons1[exon,"left"] > previousExonEnd) {
				cumulativeIntronLength <- cumulativeIntronLength + exons1[exon,"left"] - previousExonEnd - squishedIntronSize
				previousExonEnd <- exons1[exon,"right"]
			}
			if (breakpoint1 >= exons1[exon,"left"] && breakpoint1 <= exons1[exon,"right"]+1)
				breakpoint1 <- breakpoint1 - cumulativeIntronLength
			exons1[exon,"left"] <- exons1[exon,"left"] - cumulativeIntronLength
			exons1[exon,"right"] <- exons1[exon,"right"] - cumulativeIntronLength
		}

		# hide introns in gene2
		cumulativeIntronLength <- 0
		previousExonEnd <- -squishedIntronSize
		for (exon in 1:nrow(exons2)) {
			if (breakpoint2 > previousExonEnd+1 && breakpoint2 < exons2[exon,"left"])
				breakpoint2 <- (breakpoint2-previousExonEnd) / (exons2[exon,"left"]-previousExonEnd) * squishedIntronSize + previousExonEnd - cumulativeIntronLength
			if (exons2[exon,"left"] > previousExonEnd) {
				cumulativeIntronLength <- cumulativeIntronLength + exons2[exon,"left"] - previousExonEnd - squishedIntronSize
				previousExonEnd <- exons2[exon,"right"]
			}
			if (breakpoint2 >= exons2[exon,"left"] && breakpoint2 <= exons2[exon,"right"]+1)
				breakpoint2 <- breakpoint2 - cumulativeIntronLength
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

	# scale exon sizes to fit on page
	scalingFactor <- max(exons1$right) + max(exons2$right)
	if (fixedScale > 0) {
		if (fixedScale >= scalingFactor) {
			scalingFactor <- fixedScale
		} else {
			warning(paste("fallback to automatic scaling, because value for --fixedScale is too small to fit transcripts on canvas (increase it to", scalingFactor, "to avoid this)"))
		}
	}
	exons1$left <- exons1$left / scalingFactor
	exons1$right <- exons1$right / scalingFactor
	exons2$left <- exons2$left / scalingFactor
	exons2$right <- exons2$right / scalingFactor
	breakpoint1 <- breakpoint1 / scalingFactor
	breakpoint2 <- breakpoint2 / scalingFactor

	# shift gene2 to the right border of the page
	gene2Offset <- 1 + 0.05 - max(exons2$right)

	# center fusion horizontally
	fusionOffset1 <- (max(exons1$right)+gene2Offset)/2 - ifelse(fusions[fusion,"direction1"] == "downstream", breakpoint1, max(exons1$right)-breakpoint1)
	fusionOffset2 <- fusionOffset1 + ifelse(fusions[fusion,"direction1"] == "downstream", breakpoint1, max(exons1$right)-breakpoint1)

	# layout: fusion on top, circos plot on bottom left, protein domains on bottom center, statistics on bottom right
	layout(matrix(c(1,1,1,2,4,5,3,4,5), 3, 3, byrow=TRUE), widths=c(1.1, 1.2, 0.7), heights=c(1.55, 1.2, 0.25))
	par(mar=c(0, 0, 0, 0))
	plot(0, 0, type="l", xlim=c(-0.12, 1.12), ylim=c(0.4, 1.1), bty="n", xaxt="n", yaxt="n", xlab="", ylab="")

	# vertical coordinates of layers
	ySampleName <- 1.04
	yIdeograms <- ifelse(alignmentsFile != "", 0.94, 0.84)
	yBreakpointLabels <- ifelse(alignmentsFile != "", 0.86, 0.76)
	yCoverage <- 0.72
	yExons <- 0.67
	yGeneNames <- 0.58
	yFusion <- 0.5
	yTranscript <- 0.45
	yScale <- 0.407
	yTrajectoryBreakpointLabels <- yBreakpointLabels - 0.035
	yTrajectoryExonTop <- yExons + 0.03
	yTrajectoryExonBottom <- yExons - 0.055
	yTrajectoryFusion <- yFusion + 0.03

	# print sample name (title of page)
	text(0.5, ySampleName, sampleName, font=2, cex=fontSize*1.5, adj=c(0.5,0))

	# draw ideograms
	if (!is.null(cytobands)) {
		drawIdeogram("left", min(exons1$left), max(exons1$right), yIdeograms, cytobands, fusions[fusion,"contig1"], fusions[fusion,"breakpoint1"])
		drawIdeogram("right", gene2Offset, gene2Offset+max(exons2$right), yIdeograms, cytobands, fusions[fusion,"contig2"], fusions[fusion,"breakpoint2"])
	}

	# draw gene & transcript names
	if (fusions[fusion,"gene1"] != ".")
		text(max(exons1$right)/2, yGeneNames, fusions[fusion,"gene1"], font=2, cex=fontSize, adj=c(0.5,0))
	if (fusions[fusion,"site1"] != "intergenic")
		text(max(exons1$right)/2, yGeneNames-0.01, head(exons1$transcript,1), cex=0.9*fontSize, adj=c(0.5,1))
	if (fusions[fusion,"gene2"] != ".")
		text(gene2Offset+max(exons2$right)/2, yGeneNames, fusions[fusion,"gene2"], font=2, cex=fontSize, adj=c(0.5,0))
	if (fusions[fusion,"site2"] != "intergenic")
		text(gene2Offset+max(exons2$right)/2, yGeneNames-0.01, head(exons2$transcript,1), cex=0.9*fontSize, adj=c(0.5,1))

	# if multiple genes in the vicinity are shown, label them
	if (fusions[fusion,"site1"] == "intergenic")
		for (gene in unique(exons1$geneName)) {
			exonsOfGene <- exons1[exons1$geneName == gene & exons1$type != "dummy",]
			if (any(exonsOfGene$type == "exon"))
				text(mean(c(min(exonsOfGene$left), max(exonsOfGene$right))), yExons-0.04, gene, cex=0.9*fontSize, adj=c(0.5,1))
		}
	if (fusions[fusion,"site2"] == "intergenic")
		for (gene in unique(exons2$geneName)) {
			exonsOfGene <- exons2[exons2$geneName == gene & exons2$type != "dummy",]
			if (any(exonsOfGene$type == "exon"))
				text(gene2Offset+mean(c(min(exonsOfGene$left), max(exonsOfGene$right))), yExons-0.04, gene, cex=0.9*fontSize, adj=c(0.5,1))
		}

	# label breakpoints
	text(breakpoint1+0.01, yBreakpointLabels-0.03, paste0("breakpoint1\n", fusions[fusion,"display_contig1"], ":", fusions[fusion,"breakpoint1"]), adj=c(1,0), cex=fontSize)
	text(gene2Offset+breakpoint2-0.01, yBreakpointLabels-0.03, paste0("breakpoint2\n", fusions[fusion,"display_contig2"], ":", fusions[fusion,"breakpoint2"]), adj=c(0,0), cex=fontSize)

	# draw coverage axis
	if (alignmentsFile != "") {
		# left axis (gene1)
		lines(c(-0.02, -0.01, -0.01, -0.02), c(yCoverage, yCoverage, yCoverage+0.1, yCoverage+0.1))
		text(-0.025, yCoverage, "0", adj=c(1,0.5), cex=0.9*fontSize)
		text(-0.025, yCoverage+0.1, coverageNormalization1, adj=c(1,0.5), cex=0.9*fontSize)
		text(-0.05, yCoverage+0.08, "Coverage", srt=90, cex=0.9*fontSize, adj=c(1,0.5))

		# right axis (gene2)
		if (length(coverageRange) == 2) { # separate axes for gene1 and gene2
			rightCoverageAxisX <- gene2Offset+max(exons2$right)
			lines(c(rightCoverageAxisX+0.02, rightCoverageAxisX+0.01, rightCoverageAxisX+0.01, rightCoverageAxisX+0.02), c(yCoverage, yCoverage, yCoverage+0.1, yCoverage+0.1))
			text(rightCoverageAxisX+0.025, yCoverage, "0", adj=c(0,0.5), cex=0.9*fontSize)
			text(rightCoverageAxisX+0.025, yCoverage+0.1, coverageNormalization2, adj=c(0,0.5), cex=0.9*fontSize)
			text(rightCoverageAxisX+0.05, yCoverage+0.08, "Coverage", srt=90, cex=0.9*fontSize, adj=c(1,0.5))
		}

		# plot coverage 1
		rect(min(exons1$left), yCoverage, max(exons1$right), yCoverage+0.1, col="#eeeeee", border=NA)
		if (squishIntrons) {
			for (exon in 1:nrow(exons1))
				if (exons1[exon,"type"] != "CDS") # don't draw coverage twice for coding regions
					drawCoverage(exons1[exon,"left"], exons1[exon,"right"], yCoverage, coverage1, exons1[exon,"start"], exons1[exon,"end"], color1)
		} else {
			drawCoverage(min(exons1$left), max(exons1$right), yCoverage, coverage1, min(exons1$start), max(exons1$end), color1)
		}

		# plot coverage 2
		rect(gene2Offset+min(exons2$left), yCoverage, gene2Offset+max(exons2$right), yCoverage+0.1, col="#eeeeee", border=NA)
		if (squishIntrons) {
			for (exon in 1:nrow(exons2))
				if (exons2[exon,"type"] != "CDS") # don't draw coverage twice for coding regions
					drawCoverage(gene2Offset+exons2[exon,"left"], gene2Offset+exons2[exon,"right"], yCoverage, coverage2, exons2[exon,"start"], exons2[exon,"end"], color2)
		} else {
			drawCoverage(gene2Offset+min(exons2$left), gene2Offset+max(exons2$right), yCoverage, coverage2, min(exons2$start), max(exons2$end), color2)
		}
	}

	# plot gene 1
	lines(c(min(exons1$left), max(exons1$right)), c(yExons, yExons), col=darkColor1)
	for (gene in unique(exons1$geneName))
		drawStrand(min(exons1[exons1$geneName == gene,"left"]), max(exons1[exons1$geneName == gene,"right"]), yExons, darkColor1, head(exons1[exons1$geneName == gene,"strand"],1))
	for (exon in 1:nrow(exons1))
		drawExon(exons1[exon,"left"], exons1[exon,"right"], yExons, color1, exons1[exon,"exonNumber"], exons1[exon,"type"])

	# plot gene 2
	lines(c(gene2Offset, gene2Offset+max(exons2$right)), c(yExons, yExons), col=darkColor2)
	for (gene in unique(exons2$geneName))
		drawStrand(gene2Offset+min(exons2[exons2$geneName == gene,"left"]), gene2Offset+max(exons2[exons2$geneName == gene,"right"]), yExons, darkColor2, head(exons2[exons2$geneName == gene,"strand"],1))
	for (exon in 1:nrow(exons2))
		drawExon(gene2Offset+exons2[exon,"left"], gene2Offset+exons2[exon,"right"], yExons, color2, exons2[exon,"exonNumber"], exons2[exon,"type"])

	# plot gene1 of fusion
	if (fusions[fusion,"direction1"] == "downstream") {
		# plot strands
		lines(c(fusionOffset1, fusionOffset1+breakpoint1), c(yFusion, yFusion), col=darkColor1)
		for (gene in unique(exons1$geneName)) {
			exonsOfGene <- exons1[exons1$geneName == gene,]
			if (min(exonsOfGene$start) <= fusions[fusion,"breakpoint1"])
				drawStrand(fusionOffset1+min(exonsOfGene$left), fusionOffset1+min(breakpoint1, max(exonsOfGene$right)), yFusion, col=darkColor1, exonsOfGene$strand[1])
		}
		# plot exons
		for (exon in 1:nrow(exons1))
			if (exons1[exon,"start"] <= fusions[fusion,"breakpoint1"])
				drawExon(fusionOffset1+exons1[exon,"left"], fusionOffset1+min(breakpoint1, exons1[exon,"right"]), yFusion, color1, exons1[exon,"exonNumber"], exons1[exon,"type"])
		# plot trajectories
		lines(c(0, 0, fusionOffset1), c(yTrajectoryExonTop, yTrajectoryExonBottom, yTrajectoryFusion), col="red", lty=2)
		lines(c(breakpoint1, breakpoint1, fusionOffset1+breakpoint1), c(yTrajectoryBreakpointLabels, yTrajectoryExonBottom, yTrajectoryFusion), col="red", lty=2)
	} else if (fusions[fusion,"direction1"] == "upstream") {
		# plot strands
		lines(c(fusionOffset1, fusionOffset2), c(yFusion, yFusion), col=darkColor1)
		for (gene in unique(exons1$geneName)) {
			exonsOfGene <- exons1[exons1$geneName == gene,]
			if (max(exonsOfGene$end+1) >= fusions[fusion,"breakpoint1"])
				drawStrand(fusionOffset2-max(exonsOfGene$right)+breakpoint1, min(fusionOffset2, fusionOffset2-min(exonsOfGene$left)+breakpoint1), yFusion, col=darkColor1, chartr("+-", "-+", exonsOfGene$strand[1]))
		}
		# plot exons
		for (exon in 1:nrow(exons1))
			if (exons1[exon,"end"]+1 >= fusions[fusion,"breakpoint1"])
				drawExon(fusionOffset1+max(exons1$right)-exons1[exon,"right"], min(fusionOffset2, fusionOffset1+max(exons1$right)-exons1[exon,"left"]), yFusion, color1, exons1[exon,"exonNumber"], exons1[exon,"type"])
		# plot trajectories
		lines(c(max(exons1$right), max(exons1$right), fusionOffset1), c(yTrajectoryExonTop, yTrajectoryExonBottom, yTrajectoryFusion), col="red", lty=2)
		lines(c(breakpoint1, breakpoint1, fusionOffset1+max(exons1$right)-breakpoint1), c(yTrajectoryBreakpointLabels, yTrajectoryExonBottom, yTrajectoryFusion), col="red", lty=2)
	}
	
	# plot gene2 of fusion
	if (fusions[fusion,"direction2"] == "downstream") {
		# plot strands
		lines(c(fusionOffset2, fusionOffset2+breakpoint2), c(yFusion, yFusion), col=darkColor2)
		for (gene in unique(exons2$geneName)) {
			exonsOfGene <- exons2[exons2$geneName == gene,]
			if (min(exonsOfGene$start) <= fusions[fusion,"breakpoint2"])
				drawStrand(max(fusionOffset2, fusionOffset2+breakpoint2-max(exonsOfGene$right)), fusionOffset2+breakpoint2-min(exonsOfGene$left), yFusion, col=darkColor2, chartr("+-", "-+", exonsOfGene$strand[1]))
		}
		# plot exons
		for (exon in 1:nrow(exons2))
			if (exons2[exon,"start"] <= fusions[fusion,"breakpoint2"])
				drawExon(max(fusionOffset2, fusionOffset2+breakpoint2-exons2[exon,"right"]), fusionOffset2+breakpoint2-exons2[exon,"left"], yFusion, color2, exons2[exon,"exonNumber"], exons2[exon,"type"])
		# plot trajectories
		lines(c(gene2Offset, gene2Offset, fusionOffset2+breakpoint2), c(yTrajectoryExonTop, yTrajectoryExonBottom, yTrajectoryFusion), col="red", lty=2)
		lines(c(gene2Offset+breakpoint2, gene2Offset+breakpoint2, fusionOffset2), c(yTrajectoryBreakpointLabels, yTrajectoryExonBottom, yTrajectoryFusion), col="red", lty=2)
	} else if (fusions[fusion,"direction2"] == "upstream") {
		# plot strands
		lines(c(fusionOffset2, fusionOffset2+max(exons2$right)-breakpoint2), c(yFusion, yFusion), col=darkColor2)
		for (gene in unique(exons2$geneName)) {
			exonsOfGene <- exons2[exons2$geneName == gene,]
			if (max(exonsOfGene$end+1) >= fusions[fusion,"breakpoint2"])
				drawStrand(max(fusionOffset2, fusionOffset2+min(exonsOfGene$left)-breakpoint2), fusionOffset2+max(exonsOfGene$right)-breakpoint2, yFusion, col=darkColor2, exonsOfGene$strand[1])
		}
		# plot exons
		for (exon in 1:nrow(exons2))
			if (exons2[exon,"end"]+1 >= fusions[fusion,"breakpoint2"])
				drawExon(max(fusionOffset2, fusionOffset2+exons2[exon,"left"]-breakpoint2), fusionOffset2+exons2[exon,"right"]-breakpoint2, yFusion, color2, exons2[exon,"exonNumber"], exons2[exon,"type"])
		# plot trajectories
		lines(c(gene2Offset+max(exons2$right), gene2Offset+max(exons2$right), fusionOffset2+max(exons2$right)-breakpoint2), c(yTrajectoryExonTop, yTrajectoryExonBottom, yTrajectoryFusion), col="red", lty=2)
		lines(c(gene2Offset+breakpoint2, gene2Offset+breakpoint2, fusionOffset2), c(yTrajectoryBreakpointLabels, yTrajectoryExonBottom, yTrajectoryFusion), col="red", lty=2)
	}
	
	if (fusions[fusion,"fusion_transcript"] != ".") {
		# print fusion transcript colored by gene of origin
		fusion_transcript1 <- gsub("\\|.*", "", fusions[fusion,"fusion_transcript"], perl=T)
		fusion_transcript1 <- substr(fusion_transcript1, max(1, nchar(fusion_transcript1)-30), nchar(fusion_transcript1))
		fusion_transcript2 <- gsub(".*\\|", "", fusions[fusion,"fusion_transcript"], perl=T)
		fusion_transcript2 <- substr(fusion_transcript2, 1, min(nchar(fusion_transcript2), 30))
		# check for non-template bases
		non_template_bases <- gsub(".*\\|([^|]*)\\|.*", "\\1", fusions[fusion,"fusion_transcript"], perl=T)
		if (non_template_bases == fusions[fusion,"fusion_transcript"]) # no non-template bases found
			non_template_bases <- ""
		# divide non-template bases half-and-half for centered alignment
		non_template_bases1 <- substr(non_template_bases, 1, floor(nchar(non_template_bases)/2))
		non_template_bases2 <- substr(non_template_bases, ceiling(nchar(non_template_bases)/2+0.5), nchar(non_template_bases))
		# transcript 1
		text(fusionOffset2, yTranscript, bquote(.(fusion_transcript1) * phantom(.(non_template_bases1))), col=darkColor1, adj=c(1,0.5), cex=fontSize)
		# transcript 2
		text(fusionOffset2, yTranscript, bquote(phantom(.(non_template_bases2)) * .(fusion_transcript2)), col=darkColor2, adj=c(0,0.5), cex=fontSize)
		# non-template bases
		text(fusionOffset2, yTranscript, non_template_bases1, adj=c(1,0.5), cex=fontSize)
		text(fusionOffset2, yTranscript, non_template_bases2, adj=c(0,0.5), cex=fontSize)
	}

	# draw scale
	realScale <- max(exons1$end - exons1$start, exons2$end - exons2$start)
	mapScale <- max(exons1$right - exons1$left, exons2$right - exons2$left)
	# choose scale which is closest to desired scale length
	desiredScaleSize <- 0.2
	realScale <- desiredScaleSize / mapScale * realScale
	mapScale <- desiredScaleSize
	realScaleOptimalFit <- signif(realScale, 1) # round to most significant digit
	mapScaleOptimalFit <- realScaleOptimalFit / realScale * mapScale
	# draw scale line
	lines(c(1-mapScaleOptimalFit, 1), c(yScale, yScale)) # scale line
	lines(c(1-mapScaleOptimalFit, 1-mapScaleOptimalFit), c(yScale-0.007, yScale+0.007)) # left whisker
	lines(c(1, 1), c(yScale-0.007, yScale+0.007)) # right whisker
	# draw units above scale line
	realScaleThousands <- max(0, min(3, floor(log10(realScaleOptimalFit)/3)))
	scaleUnits <- c("bp", "kbp", "Mbp", "Gbp")
	scaleLabel <- paste(realScaleOptimalFit/max(1,1000^realScaleThousands), scaleUnits[realScaleThousands+1])
	text(1-mapScaleOptimalFit/2, yScale+0.005, scaleLabel, adj=c(0.5,0), cex=fontSize*0.9)
	if (squishIntrons)
		text(1-mapScaleOptimalFit/2, yScale-0.005, "introns not to scale", adj=c(0.5,1), cex=fontSize*0.9, font=3)

	# draw circos plot
	if (is.null(cytobands) || !("circlize" %in% names(sessionInfo()$otherPkgs)) || !("GenomicRanges" %in% names(sessionInfo()$otherPkgs))) {
		plot(0, 0, type="l", xlim=c(0, 1), ylim=c(0, 1), bty="n", xaxt="n", yaxt="n", xlab="", ylab="")
		plot(0, 0, type="l", xlim=c(0, 1), ylim=c(0, 1), bty="n", xaxt="n", yaxt="n", xlab="", ylab="")
	} else {
		par(mar=c(2,4,0,0), xpd=NA)
		drawCircos(fusion, fusions, cytobands, minConfidenceForCircosPlot, circosColors)
		par(mar=c(0,0,0,0), xpd=F)
	}

	# draw protein domains
	plot(0, 0, type="l", xlim=c(-0.1, 1.1), ylim=c(0, 1), bty="n", xaxt="n", yaxt="n", xlab="", ylab="")
	par(xpd=NA)
	if (!is.null(proteinDomains) && "GenomicRanges" %in% names(sessionInfo()$otherPkgs))
		drawProteinDomains(fusions[fusion,], exons1, exons2, proteinDomains, color1, color2, mergeDomainsOverlappingBy, optimizeDomainColors)
	par(xpd=F)

	# print statistics about supporting alignments
	plot(0, 0, type="l", xlim=c(0, 1), ylim=c(0, 1), bty="n", xaxt="n", yaxt="n", xlab="", ylab="")
	text(0, 0.575, "SUPPORTING READ COUNT", font=2, adj=c(0,0), cex=fontSize)
	if ("split_reads" %in% colnames(fusions)) { # STAR-Fusion reports split reads from both breakpoints combined
		text(0, 0.525, paste0("Split reads = ", fusions[fusion,"split_reads"], "\n", "Discordant mates = ", fusions[fusion,"discordant_mates"]), adj=c(0,1), cex=fontSize)
	} else { # Arriba reports split reads separately for the two breakpoints
		text(
			0, 0.525,
			paste0(
				"Split reads at breakpoint1 = ", fusions[fusion,"split_reads1"], "\n",
				"Split reads at breakpoint2 = ", fusions[fusion,"split_reads2"], "\n",
				"Discordant mates = ", fusions[fusion,"discordant_mates"]
			),
			adj=c(0,1), cex=fontSize
		)
	}

}

devNull <- dev.off()
message("Done")
