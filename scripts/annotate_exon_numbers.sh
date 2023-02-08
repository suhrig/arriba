#!/bin/bash

# parse command-line arguments
if [ $# -ne 3 ]; then
	echo Usage: $(basename "$0") fusions.tsv annotation.gtf output.tsv
	echo
	echo "Description: For each breakpoint, this script annotates the exon numbers in reference to the transcripts given in the columns 'transcript_id1' and 'transcript_id2'. It appends two columns 'exon_number1' and 'exon_number2'."
	exit 1
fi 1>&2
FUSIONS_TSV="$1"
ANNOTATION_GTF="$2"
OUTPUT_TSV="$3"

# extract only relevant GTF lines for speedup
cut -f23,24 "$FUSIONS_TSV" | tr '\t' '\n' | grep -v '^\.$' | sort -u |
grep -w -F -f /dev/stdin "$ANNOTATION_GTF" |

# annotate exon numbers
awk -F '\t' '

	# load exon coordinates from GTF file into memory
	FILENAME=="/dev/stdin" && $3=="exon" && /transcript_id / && /exon_number / {

		# extract transcript ID
		transcript=$9
		sub(/.*transcript_id[ "]*/, "", transcript)
		sub(/[;"].*/, "", transcript)

		# extract exon number
		exon=$9
		sub(/.*exon_number[ "]*/, "", exon)
		sub(/[;"].*/, "", exon)

		# cache coordinates
		transcript_id[++i]=transcript
		exon_start[i]=$4
		exon_end[i]=$5
		exon_number[i]=exon
	}

	# write header of fusion file
	FILENAME!="/dev/stdin" && FNR==1 {
		print $0"\texon_number1\texon_number2"
		for (i=1; i<=NF; i++) col[$i]=i
	}

	# annotate breakpoints with exon numbers
	FILENAME!="/dev/stdin" && FNR>1 {

		# get position
		position1=$col["breakpoint1"]; sub(/.*:/, "", position1)
		position2=$col["breakpoint2"]; sub(/.*:/, "", position2)

		# find exons overlapping breakpoints
		exon_number1=exon_number2="."
		for (i in transcript_id)
			if ($col["transcript_id1"]==transcript_id[i] && position1 >= exon_start[i]-2 && position1 <= exon_end[i]+2)
				exon_number1=exon_number[i]
		for (i in transcript_id)
			if ($col["transcript_id2"]==transcript_id[i] && position2 >= exon_start[i]-2 && position2 <= exon_end[i]+2)
				exon_number2=exon_number[i]

		print $0"\t"exon_number1"\t"exon_number2
	}
' /dev/stdin "$FUSIONS_TSV" > "$OUTPUT_TSV"

