Conceptually, Arriba is nothing more than a collection of filters. The generation of fusion candidates is entirely handled by STAR, which collects all evidence about potential gene fusions in the chimeric alignments file (and the read-through alignments file). Most of the candidates in these files are alignment artifacts, in vitro-generated artifacts, or transcript variants that are also observed in healthy tissue. Arriba applies a set of filters which try to detect artifacts based on various features that are characteristic for artifacts. A subset of the filters recovers events that were discarded by previous filters, given that the event has convincing characteristics which suggest it was discarded erroneously. All filters are enabled by default. Filters can be turned off selectively using the option `-f`. The column `filters` in the output file lists the filters that discarded the event or part of the reads supporting an event.

The filters are listed in the order in which they are applied. Read-level filters, which assess candidates based on information contained in a single read (pair), come first and are followed by event-level filters, which integrate information from multiple reads.

Read-level filters
------------------

`duplicates`
: Arriba supports two modes of duplicate marking: internal and external. Internal duplicate marking is the default. Arriba then detects PCR duplicates automatically based on identical end coordinates of fragments. The library attribute (`LB`) of the read group tag (`RG`) in the BAM header is not respected by Arriba's internal duplicate marking. External duplicate marking can be enabled using the parameter `-u`. It requires that alignments be marked as duplicates via the `BAM_FDUP` flag. External duplicate marking is advisable when the `LB` attribute should be respected or when a library with unique molecular identifiers (UMIs) is used (and other scenarios where duplicates cannot be identified merely based on mapping coordinates).

`uninteresting_contigs`
: Apart from chromosomes, genome assemblies typically contain a few unassembled contigs (`GL...`) or decoy sequences (`hs37d5`). Events between the chromosomes and these contigs are not of interest, because they are merely an effect of the incomplete state of the assembly. This filter removes all events that concern contigs other than the chromosomes 1-22, X, and Y (as defined by parameter `-i`).

`viral_contigs`
: This filter removes all events that do not involve the host chromosomes 1-22, X, and Y, but only viral contigs (as defined by parameter `-v`).

`top_expressed_viral_contigs`
: If a tumor is truly infected with a virus, a substantial number of reads should map to the respective viral contig. This filter removes viral integration candidates, unless the reads map to the top N most highly expressed viral contigs (where N is defined by parameter `-t`). An exception is made for reads that map to viral contigs with a lot of intergenic insertions into the host genome. Viral integration is a mostly random process and due to the scarcity of gene-associated regions in the genome a higher number of intergenic integration sites is to be expected. Therefore, a high ratio of intergenic-to-genic integration sites should be an indication of true viral infiltration rather than in vitro contamination or alignment artifacts.

`low_coverage_viral_contigs`
: Some viral contigs attract alignment artifacts. Often, these misaligned reads map to a focal region on the viral contig. In contrast, the real presence of a virus in a tumor manifests as fairly homogeneous expression of most of the viral contig. This filter ignores viral contigs with high focal coverage and otherwise low coverage. The coverage is considered to be insufficient if it is less than 5% of the mean coverage of the whole viral contig. If the sufficient coverage fraction is less than 15% (see parameter `-C`), all fusion candidates pertaining to the viral contig are discarded.

`read_through`
: This filter removes a fragment, when the following conditions are met:

- one mate aligns to a gene and the other outside the gene

- the read outside the gene has a distance less than given by the parameter `-R` (by default 10 kb)

- the mates are oriented in such a way that they could arise from canonical splicing

The rationale for this filter is the same as for the filter `same_gene`: These alignments could indicate a small deletion near the end of a gene, but the more likely explanation is that they arise from unannotated introns near the end of the gene or from transcription beyond the end of the gene. Obviously, the filter precludes the detection of small deletions near the end of a gene, but such deletions often produce additional aberrant isoforms which can be recognized as alignments with gaps >10 kb.

`inconsistently_clipped`
: If the insert size is small enough, it can happen that both, the first and second mate, overlap the breakpoint. Both mates should then be clipped. Occasionally, STAR clips only one of the mates, while the other one is not clipped, because it matches the reference sequence. This indicates that the clipped read should in fact not be clipped, because the clipped segment matches the reference sequence, too. This filter discards mates, when both overlap the breakpoint but only one of them is clipped, whereas the other one matches the reference sequence for at least another 3 bp.

`homopolymer`
: Candidates with breakpoints adjacent to homopolymers are discarded, because homopolymers induce alignment artifacts.

`small_insert_size`
: When the insert size is very small, the two mates overlap completely. STAR puts these reads in the chimeric alignments file as potential evidence for a duplication. This filter removes all mates whose start and end coordinates are less than 5 bp apart and which are oriented as if they arose from a duplication.

`long_gap`
: A few genes in the genome have introns of more than 1 Mbps. Some users therefore prefer to set the STAR parameter `--alignIntronMax` to a value that is high enough to accomodate these introns. Occassionally, this leads to alignment artifacts where STAR inserts a long gap of >700 kbp to align a short segment of <15 nt that otherwise would have been clipped. The filter `long_gap` discards alignments with long gaps and short aligned segments.

`same_gene`
: This filter removes a fragment, when both of its mates align to the same gene in an orientation that could arise from canonical splicing. Potentially, these alignments could indicate small intragenic deletions, but more likely they arise from splicing and should be ignored, hence.

`hairpin`
: A large fraction of the candidates found by STAR are events with a distance between the breakpoints that is smaller than the fragment size. Presumably, these are artifacts introduced during extraction, library preparation, or sequencing that arise from molecules folding back on themselves. These hairpin structures might lead to spontaneous ligations within the molecule or serve as primers for polymerases and facilitate [template switching during reverse transcription](https://doi.org/10.1371/journal.pone.0012271). When a strand-specific library is used, a lot of small duplication events are produced; unstranded libraries produce predominantly small inversions. In order to filter these probable false positives, Arriba removes fragments with a transcriptomic distance (i.e., ignoring introns) of less than the mean fragment size plus three standard deviations. The fragment size is estimated automatically or - when single-end data is supplied - a size of 200 nt is assumed. The fragment size can be overwritten via the parameter `-F`.

`mismatches`
: This filter discards alignments with a high number of reference mismatches relative to the length of the aligned segment. A binomial model is employed to determine statistical significance. The sequencing error rate is assumed to be 1%. The significance cut-off can be adjusted via the parameter `-V` (default 1%).

`low_entropy`
: This filter assesses the entropy of reads by counting the occurrences of 3-mers. If any 3-mer makes up 60% or more of the aligned segment of a read, the read is discarded. The percentage can be adjusted via the parameter `-K`. This filter removes alignment artifacts caused by simple repeats as well as by poly-A stretches as often seen in poly-A-selected libraries.

Event-level filters
-------------------

`multimappers`
: Multi-mapping reads are reduced to a single alignment. For this purpose, only the alignment with the best alignment score is retained. When there are multiple alignments with the same alignment score, the filter searches all fusion candidates associated with the multi-mapping read and determines the fusion candidate that has the most supporting reads. All alignments not associated with the most supported fusion candidate are discarded.

`merge_adjacent`
: When multiple alternative alignments with equal quality are possible, STAR does not always choose the same one. This leads to multiple breakpoints in close proximity, all of which arise from one and the same event. This filter merges all breakpoints in a window of 5 bp to the breakpoint with the highest number of supporting reads.

`non_coding_neighbors`
: Genes which are not well studied suffer from incomplete annotation. Many exons are annotated as separate genes even though they might actually be part of one and the same gene. Predicted genes named `RP11-...` are common examples for this. When poorly understood genes lie next to each other on the same strand, this would frequently lead to false positive predictions of deletions, because the transcripts that span both genes give rise to reads, which resemble focal deletions. This filter discards deletions, which are predicted between two neighboring genes, if both genes are non-coding or one breakpoint is intergenic.

`intragenic_exonic`
: Since exons usually make up only a small fraction of a gene, it is more likely that a genomic rearrangement starts and ends in intronic regions. On the transcriptomic level, this manifests as breakpoints at splice-sites or in introns. Many candidates found by STAR have both breakpoints within exons of the same gene. This is particularly true for intragenic events, which are prone to in vitro artifacts. This filter removes intragenic events if both breakpoints are in exons and more than 80% of the region between the breakpoints is intronic, such that it should be very unlikely that both breakpoints are located inside exons (see parameter `-e`).

`min_support`
: This filter discards all events with fewer reads than specified by the parameter `-S` (default 2).

`relative_support`
: With increasing expression of a gene the number of chimeric reads mapping to the gene increases, too. Arriba assumes a polynomial relationship between the number of events and the number of supporting reads for a given event. This assumption is based on empirical evidence. An expected value (e-value) is calculated for every event, reflecting how many events with the given number of supporting reads are expected by chance. This filter selects only those events with a low e-value, i.e., that have a high number of supporting reads relative to the overall number of events in a gene. Multiple covariates are taken into account, such as whether a breakpoint is at a splice-site or whether the event is intragenic. The parameters of the polynomial relationship are fixed and were estimated from several hundred RNA-Seq samples of various cancer types. The cut-off for the calculated e-value can be adjusted using the parameter `-E`. Increasing the cut-off improves sensitivity at the loss of specificity.

`internal_tandem_duplication`
: Rearrangements with both breakpoints in the same exon are affected by multiple stringent filters, because the vast majority of them are artifacts. Internal tandem duplications (ITDs), such as BCOR or FLT3 ITDs, are therefore likely to be discarded. This filter rescues ITDs that were erroneously discarded by previous filters if the number of supporting reads is overwhelming.

`intronic`
: This filter discards an event, when none of the supporting reads overlaps an exon (i.e., all alignments are in introns or intergenic regions). True events most often have at least one breakpoint in an exon or at an exon boundary.

`known_fusions`
: When a list of highly recurrent fusions is supplied (see parameter `-k`), this filter recovers events which were discarded because of too few supporting reads, as long as there is no other indication that the event might be an artifact.

`in_vitro`
: In some tissues certain genes are expressed at very high levels, for example hemoglobin and fibrinogen in blood or collagens in connective tissue. Presumably, the abundance of fragments from such genes increases the chance of unrelated molecules sticking together, which [may cause the reverse transcriptase enzyme to switch templates](https://doi.org/10.1371/journal.pone.0012271). These processes generate a large amount of chimeric fragments in vitro. Such artifactual fusions can be recognized as an extraordinary number of events with breakpoints within exons (rather than at exon boundaries, which is more common for true predictions). This filter eliminates events with genes that are highly expressed (top 0.2%) and have an unbalanced number of split-reads vs. discordant mates or that have an excessive amount of intra-exonic breakpoints.

`spliced`
: This filter recovers events discarded due to a low number of supporting reads, given that both breakpoints of the event are at splice-sites and there is at least one additional event linking the same pair of genes.

`select_best`
: If there are multiple breakpoints detected between the same pair of genes, this filter discards all but the most credible one. Events with split reads in both genes are preferred over events with only discordant mates, because in the latter case, the precise breakpoint is unknown. Moreover, events with a higher number of supporting reads are favored.

`many_spliced`
: Occassionally, a genomic rearrangement produces multiple alternatively spliced transcripts, all of which are low expressed and therefore discarded by the filters `relative_support` or `min_support`. This filter recovers events that were discarded due to too few supporting reads, under the circumstances that there are many events between a pair of genes with at least one of the breakpoints at splice-sites. All events between the pair of genes with at least one spliced breakpoint are recovered.

`no_genomic_support`
: This filter removes events with low confidence, if they are not confirmed by structural variant calls obtained from whole-genome sequencing. A file with structural variant calls can be supplied via the parameter `-d`.

`blacklist`
: This filter removes events with breakpoint coordinates matching entries in the [blacklist](input-files.md#blacklist). If an event has no split-reads (but only discordant mates), the precise breakpoint coordinates are unknown. In this case the event is discarded, if the breakpoints are within a range of the insert size around the blacklisted coordinates.

`short_anchor`
: A chimeric read aligns to some part in one of the fused genes and to some part to the other gene. The anchor of a read is the longer aligned segment. It is more likely to be aligned correctly. True positives often have anchors in both fused genes, whereas alignment artifacts are frequently characterized by only a small segment aligning to one of the genes and all anchors to the other. If the cumulative length of aligned segments is short in one of the genes, the event is discarded. The minimum length can be set using the parameter `-A` (default 23 bp).

`end_to_end`
: Theoretically, it should be impossible to observe fusions which only retain the 3' ends of the fusion partners, because the promoter would be missing. In reality, end-to-end fused transcripts are observed, albeit rarely. Most of them are false positives. Therefore, this filter discards such events, unless they have a lot of supporting reads, i.e., split-reads in both fusion partners or split-reads and discordant mates.

`no_coverage`
: For intronic and intragenic breakpoints as well as read-through fusions, this filter checks, if there is some coverage in the vicinity of the breakpoint in the normal alignments. Only reads that are not chimeric are considered. When there are no non-chimeric reads near the breakpoint, this is indicative of an alignment artifact.

`homologs`
: This filter discards events between genes that have high sequence homology, which frequently leads to erroneous alignments. Homology is quantified by counting the number of shared 16-mers. If more than 30% of the k-mers are shared between the involved genes, the event is filtered. The threshold can be defined via the parameter `-L`.

`mismappers`
: Many alignment artifacts are caused by an excessive number of reference mismatches in close proximity due to sequencing errors or adjacent SNPs. STAR tends to clip reads when it encounters a cluster of mismatches. The clipped segment is then used for a chimeric alignment, which occassionally aligns elsewhere in the genome. The filter mismappers performs a sensitive realignment of both segments. If both segments can be aligned to the same gene (while allowing more mismatches than STAR does), the chimeric alignment is considered to be an artifact. When 80% or more of the supporting reads are classified as being aligned incorrectly, the event is discarded. The threshold can be adjusted using the parameter `-m`.

`genomic_support`
: This filter recovers events which were discarded by previous filters due to few supporting reads, but which can be explained by genomic rearrangements as evidenced by structural variant calls obtained from whole-genome sequencing data. Arriba considers structural variant calls to match with breakpoints seen in transcriptomic data, when the breakpoints are less then 100 kb apart (see parameter `-D`) and the orientation of the genomic and transcriptomic breakpoints are identical.

`isoforms`
: This filter searches for additional isoforms for those gene pairs that are predicted to be fused in proper orientation. Typically there is a major isoform which is expressed at a high level and a few low expressed isoforms.

