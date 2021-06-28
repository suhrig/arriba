Alignments
----------

Arriba takes the main output file of STAR (`Aligned.out.bam`) as input (parameter `-x`). If STAR was run with the parameter `--chimOutType WithinBAM`, then this file contains all the information needed by Arriba to find fusions. When STAR was run with the parameter `--chimOutType SeparateSAMold`, the main output file lacks chimeric alignments. Instead, STAR writes them to a separate output file named `Chimeric.out.sam`. In this case, the file needs to be passed to Arriba via the parameter `-c` in addition to the main output file `Aligned.out.bam`.

Arriba extracts three types of reads from the alignment file(s):

1. Split-reads, i.e., reads composed of segments which map in a non-linear way. STAR stores such reads as supplementary alignments.

2. Discordant mates, i.e., paired-end reads which originate from the same fragment but which align in a non-linear way.

3. Alignments which cross the boundaries of annotated genes, because these alignments might arise from focal deletions. In RNA-Seq data deletions of up to several hundred kb are hard to distinguish from splicing. They are represented identically as gapped alignments, because the sizes of many introns are in fact of this order of magnitude. STAR applies a rather arbitrary measure to decide whether a gapped alignment arises from splicing or from a genomic deletion: The parameter `--alignIntronMax` determines what gap size is still assumed to be a splicing event and introns are used to represent these gaps. Only gaps larger than this limit are classified as potential evidence for genomic deletions and are stored as chimeric alignments. Most STAR-based fusion detection tools only consider chimeric alignments as evidence for gene fusions and are blind to focal deletions, hence. As a workaround, these tools recommend reducing the value of the parameter `--alignIntronMax`. But this impairs the quality of alignment, because it reduces the scope that STAR searches to find a spliced alignment. To avoid compromising the quality of alignment for the sake of fusion detection, the only solution would be to run STAR twice - once with settings optimized for regular alignment and once for fusion detection. This would double the runtime. In contrast, Arriba does not require to reduce the maximum intron size. It employs a more sensible criterion to distinguish splicing from deletions: Arriba considers all those reads as potential evidence for deletions that span the boundary of annotated genes.

The alignment files can be in SAM, BAM, and CRAM format. They need not be sorted for Arriba to accept them, but doing so comes with benefits: Often, this reduces the file size. And more importantly, the supporting reads of a fusion can be [inspected visually using a genome browser like IGV](visualization.md#inspection-of-events-using-igv), which typically requires BAM files to be sorted by coordinate.

Single-end and paired-end data and even mixtures are supported. Arriba automatically determines the data type on a read-by-read basis using the flag `BAM_FPAIRED`.

Assembly
--------

Arriba takes the assembly as input (parameter `-a`) to find mismatches between the chimeric reads and the reference genome, as well as to find alignment artifacts and homologous genes.

The script `download_references.sh` can be used to download the assembly. The available assemblies are listed when the script is run without parameters. The user is not restricted to these assemblies, however. Any assembly can be used as long as its coordinates are compatible with one of the supported assemblies (hg19/hs37d5/GRCh37 or hg38/GRCh38 or mm10/GRCm38 or mm39/GRCm39).

The assembly must be provided in FastA format and may be gzip-compressed. An index with the file extension `.fai` must exist only if CRAM files are processed.

Annotation
----------

The gene annotation (parameter `-g`) is used for multiple purposes:

- annotation of breakpoints with genes

- increased sensitivity for breakpoints at splice-sites

- calculation of transcriptomic distances

- determining the putative orientation of fused genes (i.e., 5' and 3' end)

GENCODE annotation is recommended over RefSeq annotation, because the former has a more comprehensive annotation of transcripts and splice-sites, which boosts the sensitivity. The file must be provided in GTF format and may be gzip-compressed. It does not need to be sorted.

The script `download_references.sh` can be used to download the annotation. The available annotation files are listed when the script is run without parameters. The user is not restricted to these annotation files, however. Any annotation can be used as long as its coordinates are compatible with one of the supported assemblies (hg19/hs37d5/GRCh37 or hg38/GRCh38 or mm10/GRCm38).

Blacklist
---------

It is strongly advised to run Arriba with a blacklist (parameter `-b`). Otherwise, the false positive rate increases by an order of magnitude. For this reason, using Arriba with assemblies or organisms which are not officially supported is not recommended. At the moment, the supported assemblies are: hg19/hs37d5/GRCh37, hg38/GRCh38, and mm10/GRCm38 (as well as any other assemblies that have compatible coordinates). The blacklists are contained in the [release tarballs](https://github.com/suhrig/arriba/releases) of Arriba.

The blacklist removes recurrent alignment artifacts and transcripts which are present in healthy tissue. This helps eliminate frequently observed transcripts, such as read-through fusions between neighboring genes, circular RNAs and other non-canonically spliced transcripts. It was trained on RNA-Seq samples from the [Human Protein Atlas](https://www.proteinatlas.org/), the [Illumina Human BodyMap2](https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-513/) , the [ENCODE project](https://www.encodeproject.org/) , the [Roadmap Epigenomics project](http://www.roadmapepigenomics.org/), and the [NCT MASTER cohort](https://doi.org/10.1002/ijc.30828), a heterogeneous cohort of cancer samples, from which highly recurrent artifacts were identified.

Blacklists for all supported assemblies are shipped with the download package of Arriba. They can be found in the package as `database/blacklist_*`.

The blacklist is a tab-separated file with two columns and may optionally be gzip-compressed. Lines starting with a hash (`#`) are treated as comments. Each line represents a pair of regions between which events are ignored. A region can be

- a 1-based coordinate in the format `CONTIG:POSITION`, optionally prefixed with the strand (example: `+9:56743754`). If `CONTIG` ends on an asterisk (`*`), the contig with the closest matching name is chosen.

- a range in the format `CONTIG:START-END`, optionally prefixed with a strand (example: `9:1000000-1100000`).

- the name of a gene given in the provided annotation.

In addition, special keywords are allowed for the second column:

- `any`: Discard all events if one of the breakpoints matches the given region.

- `split_read_donor`: Discard fusions only supported by split reads, if all of them have their anchor in the gene given in the first column. This filter is useful for highly mutable loci, which frequently trigger clipped alignments, such as the immunoglobulin loci or the T-cell receptor loci.

- `split_read_acceptor`: Discard events only supported by split reads, if all of them have their clipped segment in the given region.

- `split_read_any`: Discard events only supported by split reads, regardless of where the anchor is.

- `discordant_mates`: Discard fusions, if they are only supported by discordant mates (no split reads).

- `low_support`: Discard events, which have few supporting reads relative to expression (as determined by the filter `relative_support`), even if there is other evidence that the fusion might be a true positive, nonetheless. This keyword effectively prevents recovery of speculative events by filters such as `spliced` or `many_spliced`.

- `filter_spliced`: This keyword prevents the filter `spliced` from being applied to a given region. It is triggered under the same circumstances as the keyword `low_support`, but additionally requires that the breakpoints be at splice-sites for the event to be discarded. Some breakpoints produce recurrent artifacts, but the second breakpoint is always a different one, such that the pair of breakpoints is not recurrent and cannot be blacklisted. Often, such breakpoints are at splice-sites and the filter `spliced` tends to recover them. This keyword prevents the filter from doing so.

- `not_both_spliced`: This keyword discards events, unless both breakpoints are at splice-sites. This is a strict blacklist criterion, which makes sense to apply to genes which are prone to produce artifacts, because they are highly expressed, for example hemoglobins, collagens, or ribosomal genes.

- `read_through`: This keyword discards events, if they could arise from read-through transcription, i.e., the supporting reads are oriented like a deletion and are at most 400 kb apart.

Known fusions
-------------

Arriba can be instructed to be particularly sensitive towards events between certain gene pairs by supplying a list of gene pairs (parameter `-k`). A number of filters are not applied to these gene pairs. This is useful to improve the detection rate of expected or highly relevant events, such as recurrent fusions. Occassionally, this leads to false positive calls. But if high sensitivity is more important than specificity, this might be acceptable. Events which would be discarded by a filter and were recovered due to being listed in the known fusions list are usually assigned a low confidence.

Known fusions files for all supported assemblies are shipped with the download package of Arriba. They can be found in the package as `database/known_fusions_*`.

The file has two columns separated by a tab and may optionally be gzip-compressed. Lines starting with a hash (`#`) are treated as comments. Each line represents a pair of regions to which very sensitive filtering thresholds are applied. A region can be

- a 1-based coordinate in the format `CONTIG:POSITION`, optionally prefixed with the strand (example: `+9:56743754`). If `CONTIG` ends on an asterisk (`*`), the contig with the closest matching name is chosen.

- a range in the format `CONTIG:START-END`, optionally prefixed with a strand (example: `9:1000000-1100000`).

- the name of a gene given in the provided annotation.

The order of the given regions is important. The region given in the first column is assumed to denote the 5' end of the fusion and the region in the second column to be the 3' end. If Arriba cannot determine with confidence which gene constitutes the 5' and which the 3' end of a fusion prediction, then the order is ignored and the prediction is rescued in both cases.

Tags
----

Arriba can be supplied with a list of user-defined tags using the parameter `-t`. Whenever a fusion prediction matches the selection criteria for a tag, the column `tags` is populated with the respective tag. This feature is useful to annotate known oncogenic fusions, for example.

The known fusions file shipped with the download package of Arriba can be used for both known fusions and tags. It is constructed in a way that it can be passed as arguments to the parameters `-k` and `-t` alike. The former only uses the first two columns, the latter uses all three columns. If a user wants to separate filtering of known fusions and tagging of interesting fusions, different files may be used, however.

The file has three columns separated by a tab and may optionally be gzip-compressed. Lines starting with a hash (`#`) are treated as comments. Each line represents a pair of regions to be annotated. The first two columns specify the regions to be annotated; the third column the tag that is used for annotation. Some special characters in the tag are replaced with underscores (`_`) in Arriba's output file. A region can be

- a 1-based coordinate in the format `CONTIG:POSITION`, optionally prefixed with the strand (example: `+9:56743754`).

- a range in the format `CONTIG:START-END`, optionally prefixed with a strand (example: `9:1000000-1100000`).

- the name of a gene given in the provided annotation.

The order of the given regions is important. The region given in the first column is assumed to denote the 5' end of the fusion and the region in the second column to be the 3' end.

Protein domains
---------------

Protein domain annotation can be passed to Arriba via the parameter `-p`. The column `retained_protein_domains` of Arriba's output file is then populated accordingly.

Protein domain annotation files for all supported assemblies are shipped with the download package of Arriba. They can be found in the package as `database/protein_domains_*`.

The file must be in GFF3 format and may optionally be gzip-compressed. The ninth column must at least contain the following attributes:

```
Name=PROTEIN_DOMAIN_NAME;gene_id=GENE_ID;gene_name=GENE_NAME
```

The attribute `Name` is reported in the column `retained_protein_domains` of Arriba's output file. Some special characters in the name are replaced with underscores (`_`). The columns `gene_id` and `gene_name` are used to match the protein domains to the genes given in the [gene annotation](#annotation). If a match cannot be found, Arriba cannot determine the retained protein domains of the respective gene and a warning is issued. There may be many warnings if RefSeq annotation is used, because the protein domains file distributed with Arriba uses ENSEMBL gene names/IDs.

Structural variant calls from WGS
---------------------------------

If whole-genome sequencing (WGS) data is available, the sensitivity and specificity of Arriba can be improved by passing a list of structural variants detected from WGS to Arriba via the parameter `-d`. This has the following effects:

- Certain filters are overruled or run with extra sensitive settings, when an event is confirmed by WGS data.

- To reduce the false positive rate, Arriba does not report low-confidence events unless they can be matched with a structural variant found in the WGS data.

Both of these behaviors can be disabled by disabling the filters `genomic_support` and `no_genomic_support`, respectively. Providing Arriba with a list of structural variant calls then does not influence the calls, but it still has the benefit of filling the columns `closest_genomic_breakpoint1` and `closest_genomic_breakpoint2` with the breakpoints of the structural variant which is closest to a fusion. If the structural variant calls were obtained from whole-exome sequencing (WES) data rather than WGS data, the filter `no_genomic_support` should be disabled, since WES has poor coverage in most regions of the genome, such that many structural variants are missed.

Two file formats are accepted: a simple four-column format and the standard Variant Call Format (VCF). The format is detected automatically.

In case of the simple format, the file must contain four columns separated by tabs. The first two columns contain the breakpoints of the structural variants in the format `CONTIG:POSITION`. The last two columns contain the orientation of the breakpoints. The accepted values are:

- `downstream` or `+`: the fusion partner is fused downstream of the breakpoint, i.e., at a coordinate higher than the breakpoint

- `upstream` or `-`: the fusion partner is fused at a coordinate lower than the breakpoint

Example:

```
1:54420491	6:9248349	+	-
20:46703288	20:46734546	-	+
17:61499820	20:45133874	+	+
3:190967119	7:77868317	-	-
```

In case of the Variant Call Format, the file must comply with the [VCF specification for structural variants](https://samtools.github.io/hts-specs/VCFv4.2.pdf). In particular, Arriba requires that the `SVTYPE` field be present in the `INFO` column and specify one of the four values `BND`, `DEL`, `DUP`, `INV`. In addition, for all `SVTYPE`s other than `BND`, the `END` field must be present and specify the second breakpoint of the structural variant. Structural variants with single breakends are silently ignored.

Arriba checks if the orientation of the structural variant matches that of a fusion detected in the RNA-Seq data. If, for example, Arriba predicts the 5' end of a gene to be retained in a fusion, then a structural variant is expected to confirm this, or else the variant is not considered to be related.

Note: Arriba was designed for alignments from RNA-Seq data. It should not be run on WGS data directly. Many assumptions made by Arriba about the data (statistical models, blacklist, etc.) only apply to RNA-Seq data and are not valid for DNA-Seq data. For such data, a structural variant calling algorithm should be used and the results should be passed to Arriba.
