#!/bin/bash

# parse command-line arguments
if [ $# -ne 3 ]; then
	echo Usage: $(basename "$0") fusions.tsv Aligned.sortedByCoord.out.bam output_prefix
	echo
	echo "Description: This script takes fusion predictions from Arriba (fusions.tsv) and extracts the fusion-supporting alignments listed in the column 'read_identifiers' from the given input BAM file (Aligned.sortedByCoord.out.bam). The input BAM file must be sorted and indexed. For each fusion, a separate BAM file is created containing only the fusion-supporting alignments. The created BAM files are named after the given output prefix and the rank of the fusion in Arriba's output file."
	exit 1
fi 1>&2
FUSIONS_FILE="$1"
ALIGNMENTS="$2"
OUTPUT_PREFIX="$3"

# tell bash to abort on error
set -e -u -o pipefail

# make sure required software is installed
if ! [[ $(samtools --version-only 2> /dev/null) =~ ^1\. ]]; then
	echo "samtools >= 1.0 must be installed" 1>&2
	exit 1
fi

ID=0
SEARCH_WINDOW=1000000 # search this many bases up- and downstream of breakpoint for alignments
tail -n+2 "$FUSIONS_FILE" | while read LINE; do
	ID=$((ID+1))
	echo "Extracting alignments of fusion $ID"

	# extract columns from Arriba output line
	BREAKPOINT1=$(cut -f5 <<<"$LINE")
	BREAKPOINT2=$(cut -f6 <<<"$LINE")
	CHROMOSOME1="${BREAKPOINT1%:*}"
	CHROMOSOME2="${BREAKPOINT2%:*}"
	POSITION1="${BREAKPOINT1##*:}"
	POSITION2="${BREAKPOINT2##*:}"
	if [ "$POSITION1" -lt $SEARCH_WINDOW ]; then POSITION1=$SEARCH_WINDOW; fi
	if [ "$POSITION2" -lt $SEARCH_WINDOW ]; then POSITION2=$SEARCH_WINDOW; fi

	# make a new BAM file containing only the alignments of the given fusion
	(
		samtools view -H "$ALIGNMENTS"
		(
			cut -f30 <<<"$LINE" | tr ',' '\n' # extract read identifiers from Arriba output
			samtools view "$ALIGNMENTS" $CHROMOSOME1:$((POSITION1-SEARCH_WINDOW))-$((POSITION1+SEARCH_WINDOW))
			samtools view "$ALIGNMENTS" $CHROMOSOME2:$((POSITION2-SEARCH_WINDOW))-$((POSITION2+SEARCH_WINDOW))
		) |
		awk -F '\t' '
			!duplicate_line[$0]++ && NF>1 && $1 in reads{print} # extract alignments
			NF==1{reads[$0]=1} # remember names of relevant reads
		'
	) |
	samtools view -Sbu - |
	samtools sort - | tee "${OUTPUT_PREFIX}_$ID.bam" |
	samtools index - "${OUTPUT_PREFIX}_$ID.bam.bai"
done

