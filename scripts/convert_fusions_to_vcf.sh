#!/bin/bash

# parse command-line arguments
if [ $# -ne 2 ]; then
	echo Usage: $(basename "$0") input_fusions.tsv output_fusions.vcf
	echo
	echo "Description: This script converts fusion predictions from Arriba's custom tab-separated format to the standards-compliant VCF 4.3 format."
	exit 1
fi 1>&2
INPUT="$1"
OUTPUT="$2"

# tell bash to abort on error
set -e -u -o pipefail

# print VCF header
echo '##fileformat=VCFv4.3
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=MATEID,Number=.,Type=String,Description="ID of mate breakends">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO' > "$OUTPUT"

# read TSV file and convert it to VCF line-by-line
FUSION=0
tail -n +2 "$INPUT" | while read LINE; do
	FUSION=$((FUSION+1))

	BREAKPOINT1=$(cut -f5 <<<"$LINE"); CHROMOSOME1="${BREAKPOINT1%:*}"; POSITION1="${BREAKPOINT1##*:}"
	BREAKPOINT2=$(cut -f6 <<<"$LINE"); CHROMOSOME2="${BREAKPOINT2%:*}"; POSITION2="${BREAKPOINT2##*:}"
	QUAL=$(cut -f15 <<<"$LINE" | sed -e 's/low/0.5/' -e 's/medium/2/' -e 's/high/5/')
	FUSION_TRANSCRIPT=$(cut -f28 <<<"$LINE")
	REF1=$(sed -e 's/|.*//' -e 's/[^ACGT]//g' -e 's/.*\(.\)$/\1/' -e 's/./\U&/g' <<<"$FUSION_TRANSCRIPT") # get last base before fusion junction
	REF2=$(sed -e 's/.*|//' -e 's/[^ACGT]//g' -e 's/^\(.\).*/\1/' -e 's/./\U&/g' <<<"$FUSION_TRANSCRIPT") # get first base after fusion junction
	if [ -z "$REF1" ]; then REF1="."; fi
	if [ -z "$REF2" ]; then REF2="."; fi
	NON_TEMPLATE_BASES=$(sed -n -e 's/./\U&/g' -e 's/.*|\([^|]*\)|.*/\1/p' <<<"$FUSION_TRANSCRIPT") # get bases between two pipes
	DIRECTION1=$(cut -f25 <<<"$LINE")
	DIRECTION2=$(cut -f26 <<<"$LINE")
	ALT1="$REF1$NON_TEMPLATE_BASES"
	ALT2="$NON_TEMPLATE_BASES$REF2"
	if [ "$DIRECTION1" =   "upstream" ]; then ALT1=$(rev <<<"$ALT1"); fi
	if [ "$DIRECTION2" = "downstream" ]; then ALT2=$(rev <<<"$ALT2"); fi
	if [ "$DIRECTION1" = "downstream" ]; then ALT2_BREAKPOINT="]$BREAKPOINT1]"; else ALT2_BREAKPOINT="[$BREAKPOINT1["; fi
	if [ "$DIRECTION2" = "downstream" ]; then ALT1_BREAKPOINT="]$BREAKPOINT2]"; else ALT1_BREAKPOINT="[$BREAKPOINT2["; fi
	if [ "$DIRECTION1" = "downstream" ]; then ALT1="$ALT1$ALT1_BREAKPOINT"; else ALT1="$ALT1_BREAKPOINT$ALT1"; fi
	if [ "$DIRECTION2" = "downstream" ]; then ALT2="$ALT2$ALT2_BREAKPOINT"; else ALT2="$ALT2_BREAKPOINT$ALT2"; fi
	FILTER="PASS"
	INFO="SVTYPE=BND"
	echo "$CHROMOSOME1	$POSITION1	${FUSION}a	$REF1	$ALT1	$QUAL	$FILTER	$INFO;MATEID=${FUSION}b" >> "$OUTPUT"
	echo "$CHROMOSOME2	$POSITION2	${FUSION}b	$REF2	$ALT2	$QUAL	$FILTER	$INFO;MATEID=${FUSION}a" >> "$OUTPUT"
done

