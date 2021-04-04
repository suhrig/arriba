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

