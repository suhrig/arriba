#!/bin/bash

if [ $# -ne 7 ]; then
	echo "Usage: $(basename $0) STAR_genomeDir/ annotation.gtf assembly.fa blacklist.tsv read1.fastq.gz read2.fastq.gz threads" 1>&2
	exit 1
fi

# tell bash to be verbose and to abort on error
set -o pipefail
set -x -e -u

# get arguments
STAR_INDEX_DIR="$1"
ANNOTATION_GTF="$2"
ASSEMBLY_FA="$3"
BLACKLIST_TSV="$4"
READ1="$5"
READ2="$6"
THREADS="$7"

# find installation directory of arriba
BASE_DIR=$(dirname "$0")

# align FastQ files (STAR >=2.5.0a recommended)
# "--outSAMtype BAM Unsorted SortedByCoordinate" generates both, an unsorted and a coordinate-sorted output file
# the former is directly piped to extract_reads via "--outStd BAM_Unsorted"
# like so, read-through fusions are extracted while the alignment is running, instead of after
STAR \
	--runThreadN "$THREADS" \
	--genomeDir "$STAR_INDEX_DIR" --genomeLoad NoSharedMemory \
	--readFilesIn "$READ1" "$READ2" --readFilesCommand zcat \
	--outStd BAM_Unsorted --outSAMtype BAM Unsorted SortedByCoordinate \
	--outSAMunmapped Within \
	--outFilterMultimapNmax 1 --outFilterMismatchNmax 3 \
	--chimSegmentMin 10 --chimOutType SeparateSAMold --chimJunctionOverhangMin 10 --chimScoreMin 1 --chimScoreDropMax 30 --chimScoreJunctionNonGTAG 0 --chimScoreSeparation 1 --alignSJstitchMismatchNmax 5 -1 5 5 --chimSegmentReadGapMax 3 \
	--limitBAMsortRAM 50000000000 |
"$BASE_DIR/extract_reads" -g "$ANNOTATION_GTF" -r read_through.bam

# index normal alignments
samtools index Aligned.sortedByCoord.out.bam

# call arriba
"$BASE_DIR/arriba" \
	-c Chimeric.out.sam \
	-r read_through.bam \
	-x Aligned.sortedByCoord.out.bam \
	-o fusions.tsv \
	-O fusions.discarded.tsv \
	-a "$ASSEMBLY_FA" \
	-g "$ANNOTATION_GTF" \
	-b "$BLACKLIST_TSV" \
	-T -P \
#	-d structural_variants_from_WGS.tsv \
#	-k known_fusions_from_CancerGeneCensus.tsv # see section "Complete Fusion Export" at http://cancer.sanger.ac.uk/cosmic/download

