Intragenic deletions
--------------------

Arriba can detect intragenic inversions and duplications, but not deletions. This is because deletions within a gene are difficult to distinguish from ordinary splicing in RNA-Seq data. Moreover, Arriba's statistical model to find significant events is not applicable to the identification of a significant lack of exon coverage. These questions are better answered by indel callers, whole-genome sequencing, or algorithms to identify differential exon expression. For these reasons, Arriba does not report any intragenic deletions.

Small events (<fragment size)
-----------------------------

A major fraction of false positives has breakpoints of a distance of the fragment size or less. Presumably, these are artifacts introduced during library preparation. In order to maintain good specificity, Arriba currently applies strict filters to remove events with very close breakpoints. Note that Arriba uses the transcriptomic distance, i.e., introns are ignored. The genomic coordinates may therefore be much further apart because of splicing. Events which are so small that they can be represented as indel alignments, cannot be found by Arriba at all, because Arriba only considers supplementary alignments and discordant mates.

Memory consumption
------------------

Arriba usually consumes less than 10 GB of RAM. Samples with an extraordinary number of chimeric reads can require more memory. Approximately 1 GB of RAM is consumed per million chimeric reads, plus 4 GB of static overhead to load the assembly and gene annotation. Particularly multiple myeloma samples frequently exceed the normal memory requirements due to countless rearrangements in the immunoglobulin loci. In order to reduce the memory footprint, Arriba can be instructed to subsample reads, when an event has a sufficient number of supporting reads. By default, further reads are ignored, once an event has reached 300 supporting reads (see parameter `-U`).

