#!/bin/bash

if [ $# -lt 7 -o $# -gt 8 ]; then
	echo "Usage: $(basename $0) STAR_genomeDir/ annotation.gtf assembly.fa blacklist.tsv known_fusions.tsv threads read1.fastq.gz [read2.fastq.gz]" 1>&2
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
KNOWN_FUSIONS_TSV="$5"
THREADS="$6"
READ1="$7"
READ2="${8-}"

# find installation directory of arriba
BASE_DIR=$(dirname "$0")

# align FastQ files (STAR >=2.5.3a recommended)
STAR \
	--runThreadN "$THREADS" \
	--genomeDir "$STAR_INDEX_DIR" --genomeLoad NoSharedMemory \
	--readFilesIn "$READ1" "$READ2" --readFilesCommand zcat \
	--outStd BAM_Unsorted --outSAMtype BAM Unsorted --outSAMunmapped Within --outBAMcompression 0 \
	--outFilterMultimapNmax 1 --outFilterMismatchNmax 3 \
	--chimSegmentMin 10 --chimOutType WithinBAM SoftClip --chimJunctionOverhangMin 10 --chimScoreMin 1 --chimScoreDropMax 30 --chimScoreJunctionNonGTAG 0 --chimScoreSeparation 1 --alignSJstitchMismatchNmax 5 -1 5 5 --chimSegmentReadGapMax 3 |

tee Aligned.out.bam |

# call arriba
"$BASE_DIR/arriba" \
	-x /dev/stdin \
	-o fusions.tsv -O fusions.discarded.tsv \
	-a "$ASSEMBLY_FA" -g "$ANNOTATION_GTF" -b "$BLACKLIST_TSV" -k "$KNOWN_FUSIONS_TSV" -t "$KNOWN_FUSIONS_TSV" \
	-T -P \
#	-d structural_variants_from_WGS.tsv

# sorting and indexing is only required for visualization
if [[ $(samtools --version-only 2> /dev/null) =~ ^1\. ]]; then
	samtools sort -@ "$THREADS" -m 4G -T tmp -O bam Aligned.out.bam > Aligned.sortedByCoord.out.bam
	rm -f Aligned.out.bam
	samtools index Aligned.sortedByCoord.out.bam
else
	echo "samtools >= 1.0 required for sorting of alignments"
fi

