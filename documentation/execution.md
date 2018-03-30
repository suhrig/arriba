Demo script
-----------

Arriba comes with a script `run_arriba.sh` that illustrates the usage of all components of the workflow and how they are meant to interact with each other. The script is deliberatly kept simple and lacks error checking and input validation. Moreover, indexing and sorting of BAM files could be optimized if a parallelized tool such as [sambamba](http://lomereiter.github.io/sambamba/) or [biobambam](https://github.com/gt1/biobambam) was used instead of samtools. The script is meant only as a guide which demonstrates how to integrate Arriba into your own STAR-based workflow. If done properly, then fusion detection adds only a few minutes of runtime to a regular STAR alignment workflow. The following sections explain the steps of the demo script.

STAR
----

In order for `STAR` to generate chimeric alignments, the parameter `--chimSegmentMin` must be specified. In addition, the following parameters are recommended to improve sensitivity:

```bash
--chimSegmentMin 10 --chimJunctionOverhangMin 10 --chimScoreMin 1 --chimScoreDropMax 30 --chimScoreJunctionNonGTAG 0 --chimScoreSeparation 1 --alignSJstitchMismatchNmax 5 -1 5 5 --chimSegmentReadGapMax 3
```

A complete call of `STAR` may look like this:

```bash
STAR \
	--runThreadN 8 \
	--genomeDir /path/to/STAR_index --genomeLoad NoSharedMemory \
	--readFilesIn read1.fastq.gz read2.fastq.gz --readFilesCommand zcat \
	--outSAMtype BAM SortedByCoordinate \
	--outSAMunmapped Within \
	--outFilterMultimapNmax 1 --outFilterMismatchNmax 3 \
	--chimSegmentMin 10 --chimJunctionOverhangMin 10 --chimScoreMin 1 --chimScoreDropMax 30 --chimScoreJunctionNonGTAG 0 --chimScoreSeparation 1 --alignSJstitchMismatchNmax 5 -1 5 5 --chimSegmentReadGapMax 3
```

Apart from the parameters concerning chimeric alignment, all parameters may be modified according to the user's preferences.

Note: By default, STAR stores chimeric alignments in a separate file named `Chimeric.out.sam`. Alternatively, recent versions of STAR can store chimeric alignments within the file containing normal alignments (`--chimOutType WithinBAM`). Arriba cannot process chimeric alignments stored within the normal alignments. If this is the case, the chimeric alignments need to be extracted to a separate file using the parameter `-c` of the utility `extract_reads` (see section [Command-line options](command-line-options.md#extract_reads)).

extract_reads
-------------

This utility takes the normal alignments as input and extracts all fragments which cross the boundaries of genes. The result is a subset of the normal alignments, stored in a smaller BAM file. `extract_reads` can be run stand-alone, for example:

```bash
extract_reads -x Aligned.out.bam -r read_through.bam -g annotation.gtf
```

The tool can process 100 million reads in approximately 1.5 minutes on a modern CPU. For optimal performance, it is not recommended to run it in stand-alone mode, though, because doing so extends the overall runtime of the workflow. It is more efficient to run the tool in tandem with `STAR`, as illustrated in the following example, because then the entire workflow is parallellized and single-threaded periods are avoided.

In the following example, `STAR` is instructed to output the normal alignments in BAM format (`--outSAMtype BAM`) twice: sorted by coordinate into the file `Aligned.sortedByCoord.out.bam` (`--outSAMtype SortedByCoordinate`) and collated by read names to `STDOUT` (`--outStd BAM_Unsorted`). `STDOUT` is directly piped to `extract_reads`. Like so, the tool runs in parallel to alignment and does not extend the runtime as in stand-alone mode.

```bash
STAR --outStd BAM_Unsorted --outSAMtype BAM Unsorted SortedByCoordinate [...] |
extract_reads -g annotation.gtf > read_through.bam
```

arriba
------

Finally, `arriba` is run on the output files of `STAR` (`Aligned.sortedByCoord.out.bam`, `Chimeric.out.sam`) and `extract_reads` (`read_through.bam`). The normal alignments need to be sorted by coordinate and indexed so that `arriba` can lookup specific positions.

For a detailed explanation of the parameters, please refer to section [Command-line options](command-line-options.md).

