#!/bin/bash

# parse command-line arguments
if [ $# -ne 8 ]; then
	echo Usage: $(basename "$0") STAR_genomeDir/ annotation.gtf assembly.fa blacklist.tsv known_fusions.tsv protein_domains.gff3 alignments.bam
	echo
	echo "Description: This script has a similar function as the main script 'run_arriba.sh'. But instead of FastQ files, it takes existing alignments as input (alignments.bam) and realigns only those reads that are potentially relevant to fusion detection, namely unmapped and clipped reads. All other alignments are taken as they are. Since usually only a minor fraction of reads is unmapped or clipped, much of the computation time can be saved. The relevant reads are realigned using STAR, because Arriba needs SAM-compliant chimeric alignments and STAR is the only aligner that generates chimeric alignments in SAM-compliant format. This script is most useful if you want to reprocess existing alignments from an old version of STAR with a more recent version of STAR for the purpose of fusion detection in an efficient way, or if you have alignments from an entirely different aligner that does not support chimeric alignments in SAM-compliant format (e.g., HISAT2) and want to use Arriba in conjunction with it."
	echo "Note: It is not possible to mix assemblies with incompatible coordinates. For example, this script cannot be used if the prealigned BAM file (alignments.bam) is based on hg19 and the STAR index passed as argument 'STAR_genomeDir/' is based on hg38. However, it is safe to use assemblies with somewhat different contigs (as long as the main contigs are the same). For example, if the prealigned BAM file lacks viral contigs which are part of the STAR index, this script is applicable. Likewise, if the prealigned BAM file has some extra contigs which are not part of the STAR index, it is safe to use this script, because the script will simply treat reads mapping to such contigs as unmapped and will realign them."
	echo "Note: There can be slight differences in the fusion calls when using this script compared to using the main script 'run_arriba.sh'. For optimal fusion detection results, it is not recommended to use this script, but to realign all reads as is done by the main script. This script is provided as a convenience to users with a focus on computational efficiency."
	echo "Note: This script should not be used if you have generated a BAM file according to Arriba's recommended workflow. In this case, you can run Arriba directly on the BAM file. There is no need to realign any reads, because they are already aligned optimally."
	exit 1
fi 1>&2
STAR_INDEX_DIR="$1"
ANNOTATION_GTF="$2"
ASSEMBLY_FA="$3"
BLACKLIST_TSV="$4"
KNOWN_FUSIONS_TSV="$5"
TAGS_TSV="$KNOWN_FUSIONS_TSV" # different files can be used for filtering and tagging, but the provided one can be used for both
PROTEIN_DOMAINS_GFF3="$6"
THREADS=2 # there is no point in using more, because STAR never gets enough input to use more
ALIGNMENTS="$8"

# tell bash to abort on error
set -e -u -o pipefail

# make sure required software is installed
if ! [[ $(samtools --version-only 2> /dev/null) =~ ^1\. ]]; then
	echo "samtools >= 1.0 must be installed" 1>&2
	exit 1
fi
if ! [[ "$(STAR --version 2> /dev/null)" =~ ^2\.(7\.([6-9]|[1-9][0-9])|[8-9]|[1-9][0-9]) ]]; then
	echo "STAR >= 2.7.6a must be installed" 1>&2
	exit 1
fi
BASE_DIR=$(dirname "$0")
if [ -e "$BASE_DIR/../arriba" ]; then
	ARRIBA="$BASE_DIR/../arriba"
elif [ -n "$(which arriba 2> /dev/null)" ]; then
	ARRIBA=$(which arriba)
else
	echo "could not find arriba" 1>&2
	exit 1
fi

# auto-detect library layout (single-end vs. paired-end)
LAYOUT=$(samtools view "$ALIGNMENTS" | head -n1 | awk '{print ($2 % 2) ? "PE" : "SE"}' || exit 0)

(

mkfifo realign
if [ "$LAYOUT" = "SE" ]; then
	samtools view -F 2304 "$ALIGNMENTS"
elif [ "$LAYOUT" = "PE" ]; then
	samtools collate -u -f -r 1000000 -O "$ALIGNMENTS" | samtools view - |
	awk -F '\t' '$1==name { print prev; print $0 } { name=$1; prev=$0 }' # skip singletons
fi |

cut -f 1-11 | # drop SAM attributes

awk -F '\t' -v ASSEMBLY_FA="$ASSEMBLY_FA" -v LAYOUT="$LAYOUT" '
	BEGIN{
		# get list of contig names from assembly
		while (getline line < ASSEMBLY_FA)
			if (line~/^>/) {
				gsub(/^>|[ \t].*/, "", line); contig=line
			} else {
				contigs[contig]+=length($0)
			}
		# generate SAM header
		print "@HD\tVN:1.4\tSO:coordinate"
		for (contig in contigs)
			print "@SQ\tSN:"contig"\tLN:"contigs[contig]
	}

	# helper function to check if SAM flag is set
	function flag(f) { return ($2 % (2*f) >= f) }

	# find reads that need to be realigned
	function realign() {
		return (flag(4) || # unmapped
		        LAYOUT=="SE" && $6~/[0-9]{2,}S/ || # single-end and clipped
		        !flag(16) && $6~/^[0-9]{2,}S/ || # paired-end and forward strand and preclipped
		        flag(16) && $6~/[0-9]{2,}S$/ || # paired-end and reverse strand and postclipped
		        !($3 in contigs)) # contig not part of assembly
	}

	# extract paired-end reads for realignment
	LAYOUT=="PE" {
		if ($1==name1) { # we have encountered both mates
			if (realign1 || realign()) {
				print mate1 "\n" $0 > "realign"
			} else {
				print mate1 "\n" $0
			}
		} else {
			# cache mate1
			mate1=$0; name1=$1; realign1=realign()
		}
	}

	# extract single-end reads for realignment
	LAYOUT=="SE" {
		if (realign()) {
			print > "realign"
		} else {
			print
		}
	}

' & EXTRACTION=$!

# realign unmapped and clipped reads
STAR \
	--runThreadN "$THREADS" \
	--genomeDir "$STAR_INDEX_DIR" --genomeLoad NoSharedMemory \
	--readFilesIn realign --readFilesType SAM $LAYOUT \
	--outStd BAM_Unsorted --outSAMtype BAM Unsorted --outBAMcompression 0 \
	--outFilterMultimapNmax 50 --peOverlapNbasesMin 10 --alignSplicedMateMapLminOverLmate 0.5 --alignSJstitchMismatchNmax 5 -1 5 5 \
	--chimSegmentMin 10 --chimOutType WithinBAM HardClip --chimJunctionOverhangMin 10 --chimScoreDropMax 30 --chimScoreJunctionNonGTAG 0 --chimScoreSeparation 1 --chimSegmentReadGapMax 3 --chimMultimapNmax 50 > realigned.bam

wait $EXTRACTION
samtools view realigned.bam

) |

# call arriba
"$ARRIBA" \
        -x /dev/stdin \
        -o fusions.tsv -O fusions.discarded.tsv \
        -a "$ASSEMBLY_FA" -g "$ANNOTATION_GTF" -b "$BLACKLIST_TSV" -k "$KNOWN_FUSIONS_TSV" -t "$TAGS_TSV" -p "$PROTEIN_DOMAINS_GFF3"

rm -f realign realigned.bam SJ.out.tab

