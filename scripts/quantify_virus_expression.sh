#!/bin/bash

# parse command-line arguments
if [ $# -ne 2 ]; then
	echo Usage: $(basename "$0") alignments.bam virus_expression.tsv
	echo
	echo "Description: This script takes alignments in BAM format as input and counts the number of high-quality alignments mapping to viral contigs. It assumes that viral contigs start with the prefix 'AC_' or 'NC_'. An alignment is considered high-quality if it does not contain repeats, if it spans the complete read from start to end, and - in case of paired-end data - if both mates map in proper pair orientation. When reads map equally well to multiple related viral strains, only the highest expressed viral strain is reported (as measured by RPKM). Only viruses having at least 5% and at least 100bp of their genome expressed are reported."
	echo "Note: When a BAM file index exists, only the relevant contigs are read from the BAM file and processing is substantially faster."
	exit 1
fi 1>&2
INPUT="$1"
OUTPUT="$2"
VIRAL_CONTIGS="${VIRAL_CONTIGS-^[AN]C_}" # regular expression to identify viral genomes
KMER_LENGTH="${MAX_LENGTH-12}" # identify related viral strains by looking for common kmers of reads mapping to multiple viral genomes
MAX_SHARED_KMERS_PCT="${MAX_SHARED_KMERS-10}" # consider two viral strains to be related when they share this fraction (%) of kmers
MIN_COVERED_GENOME_PCT="${MIN_COVERED_GENOME_PCT-5}" # ignore viruses with only a small fraction (%) of their genome covered
MIN_COVERED_GENOME_BASES="${MIN_COVERED_GENOME_BASES-100}" # ignore viruses with only few bases of their genome covered

# tell bash to abort on error
set -e -u -o pipefail

# make sure required software is installed
if ! [[ $(samtools --version-only 2> /dev/null) =~ ^1\. ]]; then
	echo "samtools >= 1.0 must be installed" 1>&2
	exit 1
fi

# if there is a BAM index, make use of it by counting reads on non-viral contigs and collecting the names of expressed viruses
if [ -e "$INPUT.bai" ]; then
	TOTAL_MAPPED_READS=$(samtools idxstats "$INPUT" | awk -v viral_contigs="$VIRAL_CONTIGS" '!match($1,viral_contigs){sum+=$3} END{print sum}')
	EXPRESSED_VIRAL_CONTIGS=$(samtools idxstats "$INPUT" | awk -v viral_contigs="$VIRAL_CONTIGS" 'match($1,viral_contigs){print $1}')
fi

samtools view -F 4 -h "$INPUT" ${EXPRESSED_VIRAL_CONTIGS-} |

# quantify expression of viral contigs
awk -F '\t' -v OFS='\t' -v kmer_length="$KMER_LENGTH" -v max_shared_kmers_pct="$MAX_SHARED_KMERS_PCT" -v viral_contigs="$VIRAL_CONTIGS" -v min_covered_genome_pct="$MIN_COVERED_GENOME_PCT" -v min_covered_genome_bases="$MIN_COVERED_GENOME_BASES" -v total_mapped_reads="${TOTAL_MAPPED_READS-0}" '
{
	if (/^@SQ\t/) { # this is a header line => remember sizes of viral genomes

		sub(/SN:/, "", $2)
		sub(/LN:/, "", $3)
		ids[$2] = id++
		names[ids[$2]] = $2
		size[ids[$2]] = $3

	} else { # this is an alignment record => quantify expression

		total_mapped_reads++

		if (($2 % 4 >= 2 || $2 % 2 == 0) && # reads must be aligned in proper pair, unless they are single-end
		    match($3, viral_contigs) && # read must map to a viral contig
		    $6 ~ /^[0-9NMX]+$/ && # only count fully aligned reads (i.e., CIGAR string contains only N, M, or X)
		    $10 !~ /(AA.?){8}|(AC.?){8}|(AG.?){8}|(AT.?){8}|(CC.?){8}|(CG.?){8}|(CT.?){8}|(GG.?){8}|(GT.?){8}|(TT.?){8}/) { # ignore tandem repeat regions

			id = ids[$3]
			viral_mapped_reads[id]++

			# split read into kmers and make a list of which kmer belongs to which virus
			# later on, we find related viruses by searching for shared kmers
			for (i = 1; i + kmer_length <= length($10); i++) {
				kmer = substr($10, i, kmer_length)
				if (!already_seen[id,kmer]++) {
					viruses_by_kmer[kmer] = viruses_by_kmer[kmer] "," id
					kmer_count[id]++
				}
			}

			# compute covered bases of virus genome
			split($6, cigar_ops, "[0-9]+")
			split($6, cigar_lens, "[A-Z]")
			reference_pos = $4
			for (cigar = 1; cigar < length(cigar_lens); cigar++) {
				if (cigar_ops[cigar+1] ~ /[MX]/) {
					for (i = 1; i <= cigar_lens[cigar]; i++) {
						if (!already_counted[id,reference_pos]++)
							covered_genome[id]++
						reference_pos++
					}
				} else if (cigar_ops[cigar+1] ~ /[NI]/) {
					reference_pos += cigar_lens[cigar]
				}
			}

		}
	}
}

END {
	# for each virus, compute expression in RPKM
	for (virus in viral_mapped_reads)
		if (size[virus] > 0 && total_mapped_reads > 0)
			rpkm[virus] = 1000000000 * viral_mapped_reads[virus] / size[virus] / total_mapped_reads

	# remove viruses with sequence similary to viruses with higher expression
	for (kmer in viruses_by_kmer) {
		split(viruses_by_kmer[kmer], viruses, ",")
		for (i in viruses)
			for (j in viruses)
				if (viruses[i] in rpkm && viruses[j] in rpkm)
					if (rpkm[viruses[i]] > rpkm[viruses[j]] || rpkm[viruses[i]] == rpkm[viruses[j]] && viruses[i] < viruses[j]) # use virus name as tie-breaker
						if (++shared_kmers[viruses[i],viruses[j]] > kmer_count[viruses[j]] * max_shared_kmers_pct/100)
							remove[viruses[j]] = 1
	}

	# print output table sorted by RPKM
	print "VIRUS", "GENOME_SIZE", "COVERED_BASES", "COVERED_GENOME_FRACTION", "HIGH_QUALITY_ALIGNMENTS", "RPKM"
	for (virus in rpkm)
		if (!(virus in remove))
			if (covered_genome[virus] >= min_covered_genome_bases && covered_genome[virus]/size[virus] > min_covered_genome_pct/100)
				print names[virus], size[virus], covered_genome[virus], covered_genome[virus]/size[virus], viral_mapped_reads[virus], rpkm[virus] | "sort -k6,6gr"
}
' > "$OUTPUT"

