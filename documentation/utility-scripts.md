Utility scripts
===============

The [folder `scripts`](https://github.com/suhrig/arriba/tree/master/scripts) contains small utility scripts for common tasks related to fusion detection.

Extract fusion-supporting alignments
------------------------------------

**Usage:**

```
extract_fusion-supporting_alignments.sh fusions.tsv Aligned.sortedByCoord.out.bam output_prefix
```

**Description:**

This script takes fusion predictions from Arriba (`fusions.tsv`) and extracts the fusion-supporting alignments listed in the column `read_identifiers` from the given input BAM file (`Aligned.sortedByCoord.out.bam`). The input BAM file must be sorted and indexed. For each fusion, a separate mini-BAM file is created containing only the fusion-supporting alignments. The created BAM files are named after the given output prefix and the rank of the fusion in Arriba's output file.

Convert fusions.tsv to VCF
--------------------------

**Usage:**

```
convert_fusions_to_vcf.sh assembly.fa input_fusions.tsv output_fusions.vcf
```

**Description:**

This script converts fusion predictions from Arriba's custom tab-separated format to the standards-compliant [Variant Call Format version 4.3](https://samtools.github.io/hts-specs/VCFv4.3.pdf).

If a FastA index (.fai) does not exist for the given assembly file, it will be created on-the-fly.

Run Arriba on prealigned BAM file
---------------------------------

**Usage:**

```
run_arriba_on_prealigned_bam.sh STAR_genomeDir/ annotation.gtf assembly.fa blacklist.tsv known_fusions.tsv protein_domains.gff3 threads alignments.bam
```

**Description:**

This script has a similar function as the main script `run_arriba.sh`. But instead of FastQ files, it takes existing alignments as input (`alignments.bam`) and realigns only those reads that are potentially relevant to fusion detection, namely unmapped and clipped reads. All other alignments are taken as they are. Since usually only a minor fraction of reads is unmapped or clipped, much of the computation time can be saved. The relevant reads are realigned using STAR, because Arriba needs SAM-compliant chimeric alignments and STAR is the only aligner that generates chimeric alignments in SAM-compliant format.

This script is most useful if you want to efficiently reprocess existing alignments from an old version of STAR with a more recent version of STAR for the purpose of fusion detection or if you have alignments from an entirely different aligner that does not support chimeric alignments in SAM-compliant format (e.g., HISAT2) and want to use Arriba in conjunction with it.

Note: It is not possible to mix assemblies with incompatible coordinates. For example, this script cannot be used if the prealigned BAM file (`alignments.bam`) is based on hg19 and the STAR index passed as argument `STAR_genomeDir/` is based on hg38. However, it is safe to use assemblies with somewhat different contigs (as long as the main contigs are the same). For example, if the prealigned BAM file lacks viral contigs which are part of the STAR index, this script is applicable. Likewise, if the prealigned BAM file has some extra contigs which are not part of the STAR index, it is safe to use this script, because the script will simply treat reads mapping to such contigs as unmapped and will realign them.

Note: There can be slight differences in the fusion calls when using this script compared to using the main script `run_arriba.sh`. For optimal fusion detection results, it is not recommended to use this script, but to realign all reads as is done by the main script. This script is provided as a convenience to users with a focus on computational efficiency.

Note: This script should not be used if you have generated a BAM file according to Arriba's [recommended workflow](workflow.md). In this case, you can run Arriba directly on the BAM file. There is no need to realign any reads, because they are already aligned optimally.

Quantify virus expression
-------------------------

**Usage:**

```
quantify_virus_expression.sh alignments.bam virus_expression.tsv
```

**Description:**

This script takes alignments in BAM format as input and counts the number of high-quality alignments mapping to viral contigs. It assumes that viral contigs start with the prefix `AC_` or `NC_`. An alignment is considered high-quality if it does not contain repeats, if it spans the complete read from start to end, and - in case of paired-end data - if both mates map in proper pair orientation. When reads map equally well to multiple related viral strains, only the highest expressed viral strain is reported (as measured by RPKM). Only viruses having at least 5% and at least 100bp of their genome expressed are reported.

Note: When a BAM file index exists, only the relevant contigs are read from the BAM file and processing is substantially faster.

Annotate exon numbers
---------------------

**Usage:**

```
annotate_exon_numbers.sh fusions.tsv annotation.gtf output.tsv
```

**Description:**

For each breakpoint, this script annotates the exon numbers in reference to the transcripts given in the columns `transcript_id1` and `transcript_id2`. It appends two columns `exon_number1` and `exon_number2`.

