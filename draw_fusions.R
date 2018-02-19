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
	stop("usage: draw_fusions.R --annotation=annotation.gtf --fusions=fusions.tsv --output=output.pdf [--squishIntrons=TRUE] [--printStats=TRUE] [--printFusionTranscript=TRUE] [--pdfPaper=a4r] [--pdfWidth=11] [--pdfHeight=7]")
exonsFile <- parseStringParameter("annotation", args)
if (file.access(exonsFile) == -1)
	stop(sprintf("Exon annotation file (%s) does not exist", exonsFile))
fusionsFile <- parseStringParameter("fusions", args)
if (file.access(fusionsFile) == -1)
	stop(sprintf("Fusions file (%s) does not exist", fusionsFile))
outputFile <- parseStringParameter("output", args)
if (outputFile == "")
	stop("Output file not specified")
squishIntrons <- parseBooleanParameter("squishIntrons", args, T)
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

# read exon annotation
exons <- read.table(exonsFile, header=F, sep="\t", comment.char="#", quote="", stringsAsFactors=F)[,c(1, 3, 4, 5, 7, 9)]
colnames(exons) <- c("contig", "type", "start", "end", "strand", "attributes")
exons <- exons[exons$type %in% c("exon", "CDS") & grepl("exon_number ", exons$attributes),]
exons$geneName <- gsub(".*gene_name \"?([^;\"]+)\"?;.*", "\\1", exons$attributes)
exons$transcript <- gsub(".*transcript_id \"?([^;\"]+)\"?;.*", "\\1", exons$attributes)
exons$exonNumber <- gsub(".*exon_number \"?([^;\"]+)\"?;.*", "\\1", exons$attributes)
exons$contig <- sub("chr", "", sub("chrM", "MT", exons$contig))

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

drawStrands <- function(start, end, y, col, strand, scalingFactor, drawLabels=T, gene="") {
	# draw strand
	lines(c(start, end), c(y, y), col=col, lwd=2)
	# optionally draw gene label
	if (gene != "")
		text(start-0.1*scalingFactor, y, gene, srt=90)
	# indicate orientation
	if (strand %in% c("+", "-")) {
		for (i in seq(start+scalingFactor*0.01, end-scalingFactor*0.01, by=sign(end-start-2*scalingFactor*0.01)*scalingFactor*0.01))
			arrows(i, y, i+scalingFactor*0.001*ifelse(strand=="+", 1, -1), y, col=col, length=0.08, lwd=2, angle=60)
	}
}

drawExon <- function(start, end, y, color, title, type) {
	if (title != "dummy") {
		if (type == "CDS") {
			# draw coding regions as thicker bars
			rect(start, y+0.05, end, y+0.025, col=color, border=NA)
			rect(start, y-0.05, end, y-0.025, col=color, border=NA)
		} else {
			rect(start, y+0.025, end, y-0.025, col=color, border=NA)
			text((start+end)/2, y, title, cex=0.75)
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
	exons1 <- findExons(exons, fusions[fusion,"gene1"], fusions[fusion,"direction1"], fusions[fusion,"contig1"], fusions[fusion,"breakpoint1"])
	if (nrow(exons1) == 0) { plot(0, 0, type="l", xaxt="n", yaxt="n", xlab="", ylab=""); text(0, 0, paste0("Error: exon coordinates of ", fusions[fusion,"gene1"], " not found in\n", exonsFile)); next }
	exons2 <- findExons(exons, fusions[fusion,"gene2"], fusions[fusion,"direction2"], fusions[fusion,"contig2"], fusions[fusion,"breakpoint2"])
	if (nrow(exons2) == 0) { plot(0, 0, type="l", xaxt="n", yaxt="n", xlab="", ylab=""); text(0, 0, paste0("Error: exon coordinates of ", fusions[fusion,"gene2"], " not found in\n", exonsFile)); next }

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

	if (squishIntrons) {
		# hide introns in gene1
		cumulativeIntronLength <- 0
		previousExonEnd <- -200
		for (exon in 1:nrow(exons1)) {
			if (breakpoint1 > previousExonEnd+1 && breakpoint1 < exons1[exon,"start"])
				breakpoint1 <- (breakpoint1-previousExonEnd) / (exons1[exon,"start"]-previousExonEnd) * 200 + previousExonEnd - cumulativeIntronLength
			if (exons1[exon,"start"] > previousExonEnd)
				cumulativeIntronLength <- cumulativeIntronLength + exons1[exon,"start"] - previousExonEnd - 200
			if (breakpoint1 >= exons1[exon,"start"] && breakpoint1 <= exons1[exon,"end"]+1)
				breakpoint1 <- breakpoint1 - cumulativeIntronLength
			previousExonEnd <- exons1[exon,"end"]
			exons1[exon,"start"] <- exons1[exon,"start"] - cumulativeIntronLength
			exons1[exon,"end"] <- exons1[exon,"end"] - cumulativeIntronLength
		}

		# hide introns in gene2
		cumulativeIntronLength <- 0
		previousExonEnd <- -200
		for (exon in 1:nrow(exons2)) {
			if (breakpoint2 > previousExonEnd+1 && breakpoint2 < exons2[exon,"start"])
				breakpoint2 <- (breakpoint2-previousExonEnd) / (exons2[exon,"start"]-previousExonEnd) * 200 + previousExonEnd - cumulativeIntronLength
			if (exons2[exon,"start"] > previousExonEnd)
				cumulativeIntronLength <- cumulativeIntronLength + exons2[exon,"start"] - previousExonEnd - 200
			if (breakpoint2 >= exons2[exon,"start"] && breakpoint2 <= exons2[exon,"end"]+1)
				breakpoint2 <- breakpoint2 - cumulativeIntronLength
			previousExonEnd <- exons2[exon,"end"]
			exons2[exon,"start"] <- exons2[exon,"start"] - cumulativeIntronLength
			exons2[exon,"end"] <- exons2[exon,"end"] - cumulativeIntronLength
		}
	} else { # don't squish introns
		# shift exon coordinates to align the gene to the left border of the plot
		exons1$end <- exons1$end - min(exons1$start)
		breakpoint1 <- breakpoint1 - min(exons1$start)
		exons1$start <- exons1$start - min(exons1$start)
		exons2$end <- exons2$end - min(exons2$start)
		breakpoint2 <- breakpoint2 - min(exons2$start)
		exons2$start <- exons2$start - min(exons2$start)
	}

	scalingFactor <- max(exons1$end, exons2$end, ifelse(fusions[fusion,"direction1"] == "downstream", breakpoint1, max(exons1$end)-breakpoint1)+ifelse(fusions[fusion,"direction2"] == "downstream", breakpoint2, max(exons2$end)-breakpoint2))

	par(mar=c(0, 0, 0, 0))
	plot(0, 0, type="l", xlim=c(-0.12*scalingFactor, scalingFactor*1.12), ylim=c(-0.1, 1), bty="n", xaxt="n", yaxt="n")

	# plot gene 1 (at top of the page)
	drawStrands(-0.02*scalingFactor, max(exons1$end)+0.02*scalingFactor, 0.85, "darkolivegreen4", exons1[1,"strand"], scalingFactor, gene=fusions[fusion,"gene1"])
	for (exon in 1:nrow(exons1))
		drawExon(exons1[exon,"start"], exons1[exon,"end"], 0.85, "darkolivegreen2", exons1[exon,"exonNumber"], exons1[exon,"type"])

	# plot gene 2 (at the bottom of the page)
	drawStrands(-0.02*scalingFactor, max(exons2$end)+0.02*scalingFactor, 0.15, col="deepskyblue4", exons2[1,"strand"], scalingFactor, gene=fusions[fusion,"gene2"])
	for (exon in 1:nrow(exons2))
		drawExon(exons2[exon,"start"], exons2[exon,"end"], 0.15, "deepskyblue2", exons2[exon,"exonNumber"], exons2[exon,"type"])
	
	fusionOffset <- ifelse(fusions[fusion,"direction1"] == "downstream", breakpoint1, max(exons1$end)-breakpoint1)
	
	# plot gene1 of fusion (middle left of the page)
	if (fusions[fusion,"direction1"] == "downstream") {
		# plot strands
		drawStrands(-0.02*scalingFactor, breakpoint1, 0.5, col="darkolivegreen4", exons1[1,"strand"], scalingFactor, drawLabels=F, gene="fusion")
		# plot exons
		for (exon in 1:nrow(exons1))
			if (exons1[exon,"start"] <= breakpoint1)
				drawExon(exons1[exon,"start"], min(breakpoint1, exons1[exon,"end"]), 0.5, "darkolivegreen2", exons1[exon,"exonNumber"], exons1[exon,"type"])
		# plot trajectories
		lines(c(breakpoint1, breakpoint1), c(0.8, 0.55), col="red", lty=2)
		lines(c(0, 0), c(0.8, 0.55), col="red", lty=2)
		lines(c(0, 0), c(0.8, 0.90), col="red", lty=2)
	} else if (fusions[fusion,"direction1"] == "upstream") {
		# plot strands
		drawStrands(-0.02*scalingFactor, fusionOffset, 0.5, col="darkolivegreen4", chartr("+-", "-+", exons1[1,"strand"]), scalingFactor, gene="fusion", drawLabels=F)
		# plot exons
		for (exon in 1:nrow(exons1))
			if (exons1[exon,"end"]+1 >= breakpoint1)
				drawExon(max(exons1$end)-exons1[exon,"end"], min(fusionOffset, max(exons1$end)-exons1[exon,"start"]), 0.5, "darkolivegreen2", exons1[exon,"exonNumber"], exons1[exon,"type"])
		# plot trajectories
		lines(c(breakpoint1, fusionOffset), c(0.8, 0.55), col="red", lty=2)
		lines(c(max(exons1$end), 0), c(0.8, 0.55), col="red", lty=2)
		lines(c(max(exons1$end), max(exons1$end)), c(0.8, 0.90), col="red", lty=2)
	}
	# indicate breakpoint1 with dashed line
	lines(c(breakpoint1, breakpoint1), c(0.8, 0.92), col="red", lty=2)
	text(breakpoint1, 0.94, paste0("breakpoint\n", fusions[fusion,"contig1"], ":", fusions[fusion,"breakpoint1"]), cex=0.75)
	
	# plot gene2 of fusion (middle right of the page)
	if (fusions[fusion,"direction2"] == "downstream") {
		# plot strands
		drawStrands(fusionOffset, fusionOffset+breakpoint2+0.02*scalingFactor, 0.5, col="deepskyblue4", chartr("+-", "-+", exons2[1,"strand"]), scalingFactor, drawLabels=F)
		# plot exons
		for (exon in 1:nrow(exons2))
			if (exons2[exon,"start"] <= breakpoint2)
				drawExon(max(fusionOffset, fusionOffset+breakpoint2-exons2[exon,"end"]), fusionOffset+breakpoint2-exons2[exon,"start"], 0.5, "deepskyblue2", exons2[exon,"exonNumber"], exons2[exon,"type"])
		# plot trajectories
		lines(c(breakpoint2, fusionOffset), c(0.2, 0.45), col="red", lty=2)
		lines(c(0, fusionOffset+breakpoint2), c(0.2, 0.45), col="red", lty=2)
		lines(c(0, 0), c(0.2, 0.1), col="red", lty=2)
	} else if (fusions[fusion,"direction2"] == "upstream") {
		# plot strands
		drawStrands(fusionOffset, fusionOffset+max(exons2$end)-breakpoint2+0.02*scalingFactor, 0.5, col="deepskyblue4", exons2[1,"strand"], scalingFactor, drawLabels=F)
		# plot exons
		for (exon in 1:nrow(exons2))
			if (exons2[exon,"end"]+1 >= breakpoint2)
				drawExon(max(fusionOffset, fusionOffset+exons2[exon,"start"]-breakpoint2), fusionOffset+exons2[exon,"end"]-breakpoint2, 0.5, "deepskyblue2", exons2[exon,"exonNumber"], exons2[exon,"type"])
		# plot trajectories
		lines(c(breakpoint2, fusionOffset), c(0.2, 0.45), col="red", lty=2)
		lines(c(max(exons2$end), fusionOffset+max(exons2$end)-breakpoint2), c(0.2, 0.45), col="red", lty=2)
		lines(c(max(exons2$end), max(exons2$end)), c(0.2, 0.1), col="red", lty=2)
	}
	# indicate breakpoint2 with dashed line
	lines(c(breakpoint2, breakpoint2), c(0.2, 0.08), col="red", lty=2)
	text(breakpoint2, 0.06, paste0("breakpoint\n", fusions[fusion,"contig2"], ":", fusions[fusion,"breakpoint2"]), cex=0.75)
	
	if (printStats) {
		# print statistics about supporting alignments
		text(scalingFactor/2, 0, paste("Split reads in", fusions[fusion,"gene1"], "=", fusions[fusion,"split_reads1"]), cex=0.75)
		text(scalingFactor/2, -0.03, paste("Split reads in", fusions[fusion,"gene2"], "=", fusions[fusion,"split_reads2"]), cex=0.75)
		text(scalingFactor/2, -0.06, paste("Discordant mates =", fusions[fusion,"discordant_mates"]), cex=0.75)
	}

	if (printFusionTranscript && fusions[fusion,"fusion_transcript"] != ".") {
		# print fusion transcript colored by gene of origin
		fusion_transcript1 <- gsub("\\|.*", "", fusions[fusion,"fusion_transcript"], perl=T)
		fusion_transcript1 <- substr(fusion_transcript1, max(1, nchar(fusion_transcript1)-30), nchar(fusion_transcript1))
		fusion_transcript2 <- gsub(".*\\|", "", fusions[fusion,"fusion_transcript"], perl=T)
		fusion_transcript2 <- substr(fusion_transcript2, 1, min(nchar(fusion_transcript2), 30))
		text(scalingFactor/2, -0.09, bquote("Fusion transcript = " * phantom(.(fusion_transcript1)) * phantom(.(fusion_transcript2))), cex=0.75)
		text(scalingFactor/2, -0.09, bquote(phantom("Fusion transcript = ") * .(fusion_transcript1) * phantom(.(fusion_transcript2))), col="darkolivegreen4", cex=0.75)
		text(scalingFactor/2, -0.09, bquote(phantom("Fusion transcript = ") * phantom(.(fusion_transcript1)) * .(fusion_transcript2)), col="deepskyblue4", cex=0.75)
	}
}

devNull <- dev.off()
