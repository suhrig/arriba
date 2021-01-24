fusions.tsv
-------------
The file `fusions.tsv` (as specified by the parameter `-o`) contains fusions which pass all of Arriba's filters. It should be highly enriched for true predictions. The predictions are listed from highest to lowest confidence. The following paragraphs describe the columns in detail:

`gene1` and `gene2`
: `gene1` contains the gene which makes up the 5' end of the transcript and `gene2` the gene which makes up the 3' end. The order is predicted on the basis of the strands that the supporting reads map to, how the reads are oriented, and splice patterns. Both columns may contain the same gene, if the event is intragenic. If a breakpoint is in an intergenic region, Arriba lists the closest genes upstream and downstream from the breakpoint, separated by a comma. The numbers in parentheses after the closest genes state the distance to the genes. If no genes are annotated for a contig (e.g., for viral genomes), the column contains a dot (`.`).

`strand1(gene/fusion)` and `strand2(gene/fusion)`
: Each of these columns contains two values seperated by a slash. The strand before the slash reflects the strand of the gene according to the gene annotation supplied to Arriba via the parameter `-g`. If the breakpoint is in an intergenic region, the value is `.`. The value after the slash reflects the strand that is transcribed. This does not necessarily match the strand of the gene, namely when the sense strand of a gene serves as the template for transcription. Occassionally, the strand that is transcribed cannot be predicted reliably. In this case, Arriba indicates the lack of information as a dot (`.`). Arriba uses splice-patterns of the alignments to assign a read to the appropriate originating gene. If a strand-specific library was used, Arriba also evaluates the strandedness in ambiguous situations, for example, when none of the supporting reads overlaps a splice-site.

`breakpoint1` and `breakpoint2`
: The columns contain the coordinates of the breakpoints in `gene1` and `gene2`, respectively. If an event is not supported by any split reads but only by discordant mates, the coordinates given here are those of the discordant mates which are closest to the true but unknown breakpoint.

`site1` and `site2`
: These columns add information about the location of the breakpoints. Possible values are: `5' UTR`, `3' UTR`, `UTR` (overlapping with a 5' UTR as well as a 3' UTR), `CDS` (coding sequence), `exon`, `intron`, and `intergenic`. The keyword `exon` is used for non-coding genes or for ambiguous situations where the breakpoint overlaps with both a coding exon and a UTR. If the breakpoint coincides with an exon boundary, the additional keyword `splice-site` is appended.

`type`
: Based on the orientation of the supporting reads and the coordinates of breakpoints, the type of event can be inferred. Possible values are: `translocation` (between different chromosomes), `duplication`, `inversion`, and `deletion`. If genes are fused head-to-head or tail-to-tail, this is indicated as `5'-5'` or `3'-3'` respectively. Genes fused in such an orientation cannot yield a chimeric protein, since one of the genes is transcribed from the wrong strand. This type of event is equivalent to the truncation of the genes. The following types of events are flagged with an extra keyword, because they are [frequent types of false positives](interpretation-of-results.md#frequent-types-of-false-positives) and/or it is not clear if they are somatic or germline variants: Deletions with a size in the range of introns (<400kb) are flagged as `read-through`, because there is a high chance that the fusion arises from read-through transcription rather than an underlying genomic deletion. Intragenic duplications with both breakpoints at splice-sites are flagged as `non-canonical-splicing`, because the supporting reads might originate from circular RNAs, which are very abundant even in normal tissue, but manifest as duplications in RNA-Seq data. Internal tandem duplications are flagged as `ITD`. It is not always clear whether the ITDs observable in RNA-Seq data are somatic or germline variants, because ITDs are abundant in the germline and germline variants cannot be filtered effectively due to lack of a normal control.

`split_reads1` and `split_reads2`
: The number of supporting split fragments with an anchor in `gene1` or `gene2`, respectively, is given in these columns. The gene to which the longer segment of the split read aligns is defined as the anchor.

`discordant_mates`
: This column contains the number of pairs (fragments) of discordant mates (a.k.a. spanning reads or bridge reads) supporting the fusion.

`coverage1` and `coverage2`
: These two columns show the coverage near breakpoint1 and breakpoint2, respectively. The coverage is calculated as the number of fragments near the breakpoint on the side of the breakpoint that is retained in the fusion transcript. Note that the coverage calculation counts all fragments (even duplicates), whereas the columns `split_reads1`, `split_reads2`, and `discordant_mates` only count non-discarded reads. Fragments discarded due to being duplicates or other types of artifacts can be found in the column `filters`.

`confidence`
: Each prediction is assigned one of the confidences `low`, `medium`, or `high`. Several characteristics are taken into account, including: the number of supporting reads, the balance of split reads and discordant mates, the distance between the breakpoints, the type of event, whether the breakpoints are intragenic or not, and whether there are other events which corroborate the prediction, e.g. multiple isoforms or balanced translocations. See section [Interpretation of results](interpretation-of-results.md) for further advice on judging the credibility of predictions.

`reading_frame`
: This column states whether the gene at the 3' end of the fusion is fused `in-frame` or `out-of-frame`. The value `stop-codon` indicates that there is a stop codon prior to the fusion junction, such that the 3' end is not translated, even if the reading frame is preserved across the junction. The prediction of the reading frame builds on the prediction of the peptide sequence. A dot (`.`) indicates that the peptide sequence cannot be predicted, for example, because the transcript sequence could not be determined or because the breakpoint of the 5' gene does not overlap a coding region.

`tags`
: When a user-defined list of tags is provided via the parameter `-t`, this column is populated with the provided tag whenever a fusion matches the coordinates specified for the respective tag. When multiple tags match, they are separated by a comma.

`retained_protein_domains`
: If Arriba is provided with protein domain annotation using the parameter `-p`, then this column is populated with protein domains retained in the fusion. Multiple protein domains are separated by a comma. Redundant protein domains are only listed once. After every domain the fraction that is retained is stated as a percentage value in parentheses. The protein domains of the 5' and 3' genes are separated by a pipe symbol (`|`).

`closest_genomic_breakpoint1` and `closest_genomic_breakpoint2`
: When a matched whole-genome sequencing sample is available, one can feed structural variant calls obtained therefrom into Arriba (see parameter `-d`). Arriba then considers this information during fusion calling, which improves the overall accuracy. These two columns contain the coordinates of the genomic breakpoints which are closest to the transcriptomic breakpoints given in the columns `breakpoint1` and `breakpoint2`. The values in parentheses are the distances between transcriptomic and genomic breakpoints.

`gene_id1` and `gene_id2`
: These two columns state the identifiers of the fused genes as given in the `gene_id` attribute in the GTF file.

`transcript_id1` and `transcript_id2`
: For both fused genes, Arriba determines the best matching isoform that is transcribed as part of the fusion. The isoform is selected by how well its annotated exons match the splice pattern of the supporting reads of a fusion.

`direction1` and `direction2`
: These columns indicate the orientation of the fusion. A value of `downstream` means that the partner is fused downstream of the breakpoint, i.e. at a coordinate higher than the breakpoint. A value of `upstream` means the partner is fused at a coordinate lower than the breakpoint. When the prediction of the strands or of the 5' gene fails, this information gives insight into which parts of the fused genes are retained in the fusion.

`filters`
: This column lists the filters which removed one or more of the supporting reads. The section [Internal algorithm](internal-algorithm.md) describes all filters in detail. The number of filtered reads is given in parentheses after the name of the filter. The total number of supporting reads can be obtained by summing up the reads given in the columns `split_reads1`, `split_reads2`, `discordant_mates`, and `filters`. If a filter discarded the event as a whole (all reads), the number of filtered reads is not stated.

`fusion_transcript`
: This column contains the fusion transcript sequence. The sequence is assembled from the supporting reads of the most highly expressed transcript. It represents the transcript isoform that is most likely expressed according to the splice patterns of the supporting reads. The column contains a dot (`.`), when the sequence could not be predicted. This is the case when the strands or the 5' end of the transcript could not be predicted reliably. The breakpoint is represented as a pipe symbol (`|`). When non-template bases are inserted between the fused genes, these bases are represented as lowercase letters between two pipes. Reference mismatches (SNPs or SNVs) are indicated as lowercase letters, insertions as bases between brackets (`[` and `]`), deleted bases as one or more dashes (`-`), introns as three underscores (`___`), and ambiguous positions, such as positions with diverse reference mismatches, are represented as `?`. Missing information due to insufficient coverage is denoted as an ellipsis (`...`). If the switch `-I` is used, then an attempt is made to fill missing information with the assembly sequence. A sequence stretch that was taken from the assembly sequence rather than the supporting reads is wrapped in parentheses (`(` and `)`). In addition, when `-I` is used, the sequence is trimmed to the boundaries of the fused transcripts. The coordinate of the fusion breakpoint relative to the start of the transcript can thus easily be inferred by counting the bases from the beginning of the fusion transcript to the breakpoint character (`|`). In case the full sequence could be constructed from the combined information of supporting reads and assembly sequence, the start of the fusion transcript is marked by a caret sign (`^`) and the end by a dollar sign (`$`). If the full sequence could not be constructed, these signs are missing.

`peptide_sequence`
: This column contains the fusion peptide sequence. The sequence is translated from the fusion transcript given in the column `fusion_transcript` and determines the reading frame of the fused genes according to the transcript isoforms given in the columns `transcript_id1` and `transcript_id2`. Translation starts at the start of the assembled fusion transcript or when the start codon is encountered in the 5' gene. Translation ends when either the end of the assembled fusion transcript is reached or when a stop codon is encountered. If the fusion transcript contains an ellipsis (`...`), the sequence beyond the ellipsis is trimmed before translation, because the reading frame cannot be determined reliably. The column contains a dot (`.`), when the transcript sequence could not be predicted or when the precise breakpoints are unknown due to lack of split reads or when the fusion transcript does not overlap any coding exons in the 5' gene or when no start codon could be found in the 5' gene or when there is a stop codon prior to the fusion junction (in which case the column `reading_frame` contains the value `stop-codon`). The breakpoint is represented as a pipe symbol (`|`). If a codon spans the breakpoint, the amino acid is placed on the side of the breakpoint where two of the three bases reside. Codons resulting from non-template bases are flanked by two pipes. Amino acids are written as lowercase characters in the following situations: non-silent SNVs/SNPs, insertions, frameshifts, codons spanning the breakpoint, non-coding regions (introns/intergenic regions/UTRs), and non-template bases. Codons which cannot be translated to amino acids, such as those having invalid characters, are represented as `?`.

`read_identifiers`
: This column contains the names of the supporting reads separated by commas.

fusions.discarded.tsv
-----------------------

The file `fusions.discarded.tsv` (as specified by the parameter `-O`) contains all events that Arriba classified as an artifact or that are also observed in healthy tissue. It has the same format as the file `fusions.tsv`. This file may be useful, if one suspects that an event should be present, but was erroneously discarded by Arriba.

