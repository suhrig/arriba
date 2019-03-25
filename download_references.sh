#!/bin/bash

declare -A ASSEMBLIES
ASSEMBLIES[hs37d5]="ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz"
ASSEMBLIES[hg19]="http://hgdownload.cse.ucsc.edu/goldenpath/hg19/bigZips/chromFa.tar.gz"
ASSEMBLIES[GRCh37]="ftp://ftp.ensembl.org/pub/grch37/release-87/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz"
ASSEMBLIES[hg38]="http://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.chromFa.tar.gz"
ASSEMBLIES[GRCh38]="ftp://ftp.ensembl.org/pub/release-93/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"

declare -A ANNOTATIONS
ANNOTATIONS[GENCODE19]="ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz"
ANNOTATIONS[RefSeq_hg19]="http://hgdownload.cse.ucsc.edu/goldenpath/hg19/database/refGene.txt.gz"
ANNOTATIONS[ENSEMBL87]="ftp://ftp.ensembl.org/pub/grch37/release-87/gtf/homo_sapiens/Homo_sapiens.GRCh37.87.chr.gtf.gz"
ANNOTATIONS[GENCODE28]="ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_28/gencode.v28.annotation.gtf.gz"
ANNOTATIONS[RefSeq_hg38]="http://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/refGene.txt.gz"
ANNOTATIONS[ENSEMBL93]="ftp://ftp.ensembl.org/pub/release-93/gtf/homo_sapiens/Homo_sapiens.GRCh38.93.chr.gtf.gz"

declare -A COMBINATIONS
COMBINATIONS["hs37d5+GENCODE19"]="hs37d5+GENCODE19"
COMBINATIONS["hs37d5+RefSeq"]="hs37d5+RefSeq_hg19"
COMBINATIONS["hs37d5+ENSEMBL87"]="hs37d5+ENSEMBL87"
COMBINATIONS["hg19+GENCODE19"]="hg19+GENCODE19"
COMBINATIONS["hg19+RefSeq"]="hg19+RefSeq_hg19"
COMBINATIONS["hg19+ENSEMBL87"]="hg19+ENSEMBL87"
COMBINATIONS["GRCh37+GENCODE19"]="GRCh37+GENCODE19"
COMBINATIONS["GRCh37+RefSeq"]="GRCh37+RefSeq_hg19"
COMBINATIONS["GRCh37+ENSEMBL87"]="GRCh37+ENSEMBL87"
COMBINATIONS["hg38+GENCODE28"]="hg38+GENCODE28"
COMBINATIONS["hg38+RefSeq"]="hg38+RefSeq_hg38"
COMBINATIONS["hg38+ENSEMBL93"]="hg38+ENSEMBL93"
COMBINATIONS["GRCh38+GENCODE28"]="GRCh38+GENCODE28"
COMBINATIONS["GRCh38+RefSeq"]="GRCh38+RefSeq_hg38"
COMBINATIONS["GRCh38+ENSEMBL93"]="GRCh38+ENSEMBL93"

if [ $# -gt 2 ] || [ -z "$1" ] || [ -z "${COMBINATIONS[$1]}" ] || ! [[ $2 =~ ^(|[1-9][0-9]*)$ ]]; then
	echo "Usage: $(basename $0) ASSEMBLY+ANNOTATION [THREADS]" 1>&2
	echo "Available assemblies and annotations: ${!COMBINATIONS[@]}" 1>&2
	exit 1
fi

ASSEMBLY="${COMBINATIONS[$1]%+*}"
ANNOTATION="${COMBINATIONS[$1]#*+}"
THREADS="${2-8}"

set -o pipefail
set -e -u

echo "Downloading assembly: ${ASSEMBLIES[$ASSEMBLY]}"
wget -q -O - "${ASSEMBLIES[$ASSEMBLY]}" |
if [[ ${ASSEMBLIES[$ASSEMBLY]} =~ \.tar\.gz$ ]]; then
	tar -x -O -z
elif [[ ${ASSEMBLIES[$ASSEMBLY]} =~ \.gz$ ]]; then
	gunzip -c
else
	cat
fi > "$ASSEMBLY.fa"

echo "Downloading annotation: ${ANNOTATIONS[$ANNOTATION]}"
wget -q -O - "${ANNOTATIONS[$ANNOTATION]}" |
if [[ ${ANNOTATIONS[$ANNOTATION]} =~ \.gz$ ]]; then
	gunzip -c
else
	cat
fi |
if [[ $ANNOTATION =~ RefSeq ]]; then
	# convert genePred to GTF
	awk -F '\t' -v OFS='\t' '
	function min(x, y) { return (x>y) ? y : x }
	function max(x, y) { return (x<y) ? y : x }
	{
		split($10, start, ",")
		split($11, end, ",")
		split($16, frame, ",")
		# remove stop codon from left end for coding genes on the minus strand
		if ($4=="-" && $14=="cmpl" && (start[1]!=$7 || (min(end[1],$8)-start[1]+frame[1])%3==0)) {
			$7+=3
			for (i in end)
				if ($7>=end[i] && $7<=end[i]+2)
					$7+=start[i+1]-end[i]
		}
		# remove stop codon from right end for coding genes on the plus strand
		if ($4=="+" && $15=="cmpl" && (end[$9]!=$8 || (end[$9]-max(start[$9],$7)+frame[$9])%3==0)) {
			$8-=3
			for (i in start)
				if ($8<=start[i] && $8>=start[i]-2)
					$8-=start[i]-end[i-1]
		}
		# append running number to duplicate uses of the same transcript ID
		if (transcripts[$2]++)
			$2=$2"_"transcripts[$2]
		# print one line for each exon
		for (i=1; i<=$9; i++) {
			exon=($4=="+") ? i : $9-i+1
			attributes="gene_id \""$13"\"; transcript_id \""$2"\"; exon_number \""exon"\"; exon_id \""$2"."exon"\"; gene_name \""$13"\";"
			print $3,"RefSeq","exon",start[i]+1,end[i],".",$4,".",attributes
			# print one line for each coding region
			if ($14~/cmpl/ && $7<=end[i] && $8>=start[i])
				print $3,"RefSeq","CDS",max($7,start[i])+1,min($8,end[i]),".",$4,frame[i],attributes
		}
	}' | sort -k1,1V -k4,4n -k5,5n -k3,3 -S4G
else
	cat
fi |
if ! grep -q -P '^>chr' "$ASSEMBLY.fa"; then
	sed -e 's/^chrM\t/MT\t/' -e 's/^chr//'
else
	sed -e 's/^MT\t/chrM\t/' -e 's/^\([1-9XY]\|[12][0-9]\)\t/\1\t/'
fi > "$ANNOTATION.gtf"

mkdir STAR_index_${ASSEMBLY}_${ANNOTATION}
STAR --runMode genomeGenerate --genomeDir STAR_index_${ASSEMBLY}_${ANNOTATION} --genomeFastaFiles "$ASSEMBLY.fa" --sjdbGTFfile "$ANNOTATION.gtf" --runThreadN "$THREADS" --sjdbOverhang 200

