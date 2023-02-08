#!/bin/bash

# parse command-line arguments
if [ $# -ne 3 ]; then
	echo Usage: $(basename "$0") assembly.fa input_fusions.tsv output_fusions.vcf
	echo
	echo "Description: This script converts fusion predictions from Arriba's custom tab-separated format to the standards-compliant VCF 4.3 format."
	echo "Note: If a FastA index (.fai) does not exist for the given assembly file, it will be created on-the-fly."
	exit 1
fi 1>&2
ASSEMBLY="$1"
INPUT=$(cat "$2")
OUTPUT="$3"

# tell bash to abort on error
set -e -u -o pipefail

# make sure required software is installed
if ! [[ $(samtools --version-only 2> /dev/null) =~ ^1\. ]]; then
	echo "samtools >= 1.0 must be installed" 1>&2
	exit 1
fi

# create FastA index, if necessary
if [ ! -e "$ASSEMBLY.fai" ]; then
	echo "Indexing FastA file" 1>&2
	samtools faidx "$ASSEMBLY"
fi

# print VCF header
echo '##fileformat=VCFv4.3' > "$OUTPUT"
echo "$INPUT" | cut -f5-6 | tr '\t' '\n' | sed -e 's/:[0-9]*$//' | sort -u |
awk -F '\t' '
	FILENAME == "/dev/stdin" { contigs[$0] }
	FILENAME != "/dev/stdin" && $1 in contigs { print "##contig=<ID="$1",length="$2">" }
' /dev/stdin "$ASSEMBLY.fai" >> "$OUTPUT"
echo '##FILTER=<ID=PASS,Description="All filters passed">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=MATEID,Number=.,Type=String,Description="ID of mate breakends">
##INFO=<ID=GENE_NAME,Number=.,Type=String,Description="Name of gene hit by breakpoint">
##INFO=<ID=GENE_ID,Number=.,Type=String,Description="ID of gene hit by breakpoint">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO' >> "$OUTPUT"

# read TSV file and convert it to VCF line-by-line
FUSION=0
tail -n +2 <<<"$INPUT" | while read LINE; do
	FUSION=$((FUSION+1))
	SITE1=$(cut -f7 <<<"$LINE")
	SITE2=$(cut -f8 <<<"$LINE")
	GENE_NAME1=$([ "$SITE1" = "intergenic" ] || cut -f1 <<<"$LINE")
	GENE_NAME2=$([ "$SITE2" = "intergenic" ] || cut -f2 <<<"$LINE")
	GENE_ID1=$([ "$SITE1" = "intergenic" ] || cut -f21 <<<"$LINE")
	GENE_ID2=$([ "$SITE2" = "intergenic" ] || cut -f22 <<<"$LINE")
	BREAKPOINT1=$(cut -f5 <<<"$LINE"); CHROMOSOME1="${BREAKPOINT1%:*}"; POSITION1="${BREAKPOINT1##*:}"
	BREAKPOINT2=$(cut -f6 <<<"$LINE"); CHROMOSOME2="${BREAKPOINT2%:*}"; POSITION2="${BREAKPOINT2##*:}"
	QUAL=$(cut -f15 <<<"$LINE" | sed -e 's/low/0.5/' -e 's/medium/2/' -e 's/high/5/')
	REF1=$(samtools faidx "$ASSEMBLY" "$BREAKPOINT1-$POSITION1" | tail -n1 | tr '[:lower:]' '[:upper:]')
	REF2=$(samtools faidx "$ASSEMBLY" "$BREAKPOINT2-$POSITION2" | tail -n1 | tr '[:lower:]' '[:upper:]')
	NON_TEMPLATE_BASES=$(cut -f28 <<<"$LINE" | tr '[:lower:]' '[:upper:]' | sed -n -e 's/.*|\([^|]*\)|.*/\1/p') # get bases between two pipes in transcript sequence
	STRAND1=$(cut -f3 <<<"$LINE" | cut -f2 -d/)
	if [ "$STRAND1" = "-" ]; then NON_TEMPLATE_BASES=$(tr ATCG TAGC <<<"$NON_TEMPLATE_BASES"); fi # complement bases on reverse strand
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
	INFO1="SVTYPE=BND;MATEID=${FUSION}b;GENE_NAME=$GENE_NAME1;GENE_ID=$GENE_ID1"
	INFO2="SVTYPE=BND;MATEID=${FUSION}a;GENE_NAME=$GENE_NAME2;GENE_ID=$GENE_ID2"
	echo "$CHROMOSOME1	$POSITION1	${FUSION}a	$REF1	$ALT1	$QUAL	$FILTER	$INFO1" >> "$OUTPUT"
	echo "$CHROMOSOME2	$POSITION2	${FUSION}b	$REF2	$ALT2	$QUAL	$FILTER	$INFO2" >> "$OUTPUT"
done

