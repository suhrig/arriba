The folder `scripts` contains small utility scripts for common tasks related to fusion detection.

Extracting fusion-supporting alignments
---------------------------------------

**Usage:**

```
extract_fusion-supporting_alignments.sh fusions.tsv Aligned.sortedByCoord.out.bam output_prefix
```

**Description:**

This script takes fusion predictions from Arriba (`fusions.tsv`) and extracts the fusion-supporting alignments listed in the column `read_identifiers` from the given input BAM file (`Aligned.sortedByCoord.out.bam`). The input BAM file must be sorted and indexed. For each fusion, a separate BAM file is created containing only the fusion-supporting alignments. The created BAM files are named after the given output prefix and the rank of the fusion in Arriba's output file.

Converting fusions.tsv to VCF format
------------------------------------

**Usage:**

```
convert_fusions_to_vcf.sh input_fusions.tsv output_fusions.vcf
```

**Description:**

This script converts fusion predictions from Arriba's custom tab-separated format to the standards-compliant [VCF 4.3 format](https://samtools.github.io/hts-specs/VCFv4.3.pdf).

Run Arriba on prealigned BAM file
---------------------------------

**Usage:**

```
run_arriba_on_prealigned_bam.sh STAR_genomeDir/ annotation.gtf assembly.fa blacklist.tsv known_fusions.tsv protein_domains.gff3 alignments.bam
```

**Description:**

This script has a similar function as the main script `run_arriba.sh`. But instead of FastQ files, it takes existing alignments as input (`alignments.bam`) and realigns only those reads that are potentially relevant to fusion detection, namely unmapped and clipped reads. All other alignments are taken as they are. Since usually only a minor fraction of reads is unmapped or clipped, much of the computation time can be saved. The relevant reads are realigned using STAR, because Arriba needs SAM-compliant chimeric alignments and STAR is the only aligner that generates chimeric alignments in SAM-compliant format.

This script is most useful if you want to reprocess existing alignments from an old version of STAR with a more recent version of STAR for the purpose of fusion detection in an efficient way, or if you have alignments from an entirely different aligner that does not support chimeric alignments in SAM-compliant format (e.g., HISAT2) and want to use Arriba in conjunction with it.

Note: It is not possible to mix assemblies with incompatible coordinates. For example, this script cannot be used if the prealigned BAM file (alignments.bam) is based on hg19 and the STAR index passed as argument `STAR_genomeDir/` is based on hg38. However, it is safe to use assemblies with somewhat different contigs (as long as the main contigs are the same). For example, if the prealigned BAM file lacks viral contigs which are part of the STAR index, this script is applicable. Likewise, if the prealigned BAM file has some extra contigs which are not part of the STAR index, it is safe to use this script, because the script will simply treat reads mapping to such contigs as unmapped and will realign them.

Note: There can be slight differences in the fusion calls when using this script compared to using the main script `run_arriba.sh`. For optimal fusion detection results, it is not recommended to use this script, but to realign all reads as is done by the main script. This script is provided as a convenience to users with a focus on computational efficiency.

Note: This script should not be used if you have generated a BAM file according to Arriba's [recommended workflow](workflow.md). In this case, you can run Arriba directly on the BAM file. There is no need to realign any reads, because they are already aligned optimally.

