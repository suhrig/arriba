Intragenic deletions
--------------------

Arriba can detect intragenic inversions and duplications, but not deletions. This is because deletions within a gene are difficult to distinguish from ordinary splicing in RNA-Seq data. Moreover, Arriba's statistical model to find significant events is not applicable to the identification of a significant lack of exon coverage. These questions are better answered by indel callers, whole-genome sequencing, or algorithms to identify differential exon expression. For these reasons, Arriba does not report any intragenic deletions.

RefSeq annotation
-----------------

It is recommended to use annotation from GENCODE or ENSEMBL to run Arriba. RefSeq annotation has less comprehensive annotation of splice sites. Moreover, RefSeq does not annotate the immunoglobulin/T-cell receptor loci. These shortcomings reduce the sensitivity of fusion detection. Users who want to use RefSeq nonetheless are advised to copy the immunoglobulin/T-cell receptor loci annotation from GENCODE/ENSEMBL as a workaround.

Memory consumption
------------------

Arriba usually consumes less than 10 GB of RAM. Samples with an extraordinary number of chimeric reads can require more memory. Approximately 1 GB of RAM is consumed per million chimeric read pairs, plus 4 GB of static overhead to load the assembly and gene annotation. Particularly multiple myeloma samples frequently exceed the normal memory requirements due to countless rearrangements in the immunoglobulin loci. In order to reduce the memory footprint, Arriba can be instructed to subsample reads, when an event has a sufficient number of supporting reads. By default, further reads are ignored, once an event has reached 300 supporting reads (see parameter `-U`).

