Multimapping reads
------------------

When a clipped segment aligns to multiple loci, it is not stored in the chimeric alignments file by STAR. Arriba cannot detect events that are solely based on multi-mapping reads, hence.

Viral integration sites
-----------------------

Fundamentally, fusion detection and detection of viral integration sites are based on the same methods (searching for split reads and discordant mates). Aligning a RNA-Seq sample against concatenated genomes consisting of human and pathogenic DNA is a conceivable approach. However, since there is a great degree of homology between related strains of pathogens, many of the supporting reads would map to multiple contigs. As multi-mapping chimeric reads are not reported by STAR, a lot of information would be invisible to Arriba and the detection rate would suffer.

At best, Arriba might be usable to detect the integration sites of a single pathogen of interest or a collection of pathogens with dissimilar genomes. This use case has not been tested with Arriba, but any feedback from adventurous users is welcome!

Intragenic deletions
--------------------

Arriba can detect intragenic inversions and duplications, but not deletions. This is because deletions within a gene are difficult to distinguish from ordinary splicing in RNA-Seq data. Moreover, Arriba's statistical model to find significant events is not applicable to the identification of a significant lack of exon coverage. These questions are better answered by indel callers, whole-genome sequencing, or algorithms to identify differential exon expression. For these reasons, Arriba does not report any intragenic deletions.

Small events (<fragment size)
-----------------------------

A major fraction of false positives has breakpoints of a distance of the fragment size or less. Presumably, these are artifacts introduced during library preparation. In order to maintain good specificity, Arriba currently does not report any events whose breakpoints are within a distance of the fragment size plus three standard deviations. Such events are removed by the filter `hairpin`. Note that Arriba uses the transcriptomic distance, i.e., introns are ignored. The genomic coordinates may therefore be much further apart because of splicing.

Memory consumption
------------------

Arriba usually consumes less than 10 GB of RAM. Samples with an extraordinary number of chimeric reads can require more memory. Approximately 1 GB of RAM is consumed per million reads in the chimeric alignments file or the read-through alignments file, plus 4 GB of static overhead to load the assembly and gene annotation. Particularly multiple myeloma samples frequently exceed the normal memory requirements due to countless rearrangements in the immunoglobulin loci. In order to reduce the memory footprint, Arriba can be instructed to subsample reads, when an event has a sufficient number of supporting reads. By default, further reads are ignored, once an event has reached 300 supporting reads (see parameter `-U`).

