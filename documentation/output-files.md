fusions.tsv
-------------
The file `fusions.tsv` (as specified by the parameter `-o`) contains fusions which pass all of Arriba's filters. It should be highly enriched for true predictions. The predictions are listed from highest to lowest confidence. The following paragraphs describe the columns in detail:

`gene1` and `gene2`
: `gene1` contains the gene which makes up the 5' end of the transcript and `gene2` the gene which makes up the 3' end. The order is predicted on the basis of the strands that the supporting reads map to, how the reads are oriented, and splice patterns. Both columns may contain the same gene, if the event is intragenic. If a breakpoint is in an intergenic region, Arriba lists the closest genes upstream and downstream from the breakpoint, separated by a comma. The numbers in parantheses after the closest genes state the distance to the genes.

`strand1(gene/fusion)` and `strand2(gene/fusion)`
: Each of these columns contains two values seperated by a slash. The strand before the slash reflects the strand of the gene according to the gene annotation supplied to Arriba via the parameter `-g`. If the breakpoint is in an intergenic region, the value is `.`. The value after the slash reflects the strand that is transcribed. This does not necessarily match the strand of the gene, namely when the opposite strand of a gene is transcribed. Occassionally, the strand that is transcribed cannot be predicted reliably. In this case, Arriba indicates the lack of information as a `.`. Arriba uses splice-patterns of the alignments to assign a read to the appropriate originating gene. If a strand-specific library was used, Arriba also evaluates the strandedness in ambiguous situations, for example, when none of the supporting reads overlaps a splice-site.

`breakpoint1` and `breakpoint2`
: The columns contain the coordinates of the breakpoints in gene1 and gene2, respectively. If an event is not supported by any split reads but only by discordant mates, the coordinates given here are those of the discordant mates which are closest to the true but unknown breakpoint.

`site1` and `site2`
: These columns add information about the location of the breakpoints. Possible values are: `splice-site` (at an exon boundary and oriented in way that the transcript has likely been spliced), `exon` (inside an exon, but not at an exon boundary), `intron`, `5' UTR`, `3' UTR`, `UTR` (overlapping with a 5' UTR as well as a 3' UTR), and `intergenic`.

`type`
: Based on the orientation of the supporting reads and the coordinates of breakpoints, the type of event can be inferred. Possible values are: `translocation` (between different chromosomes), `duplication`, `inversion`, and `deletion`. If genes are fused head-to-head or tail-to-tail, this is indicated as `5'-5'` or `3'-3'` respectively. Genes fused in such an orientation cannot yield a chimeric protein, since one of the genes is transcribed from the wrong strand. This type of event is equivalent to the truncation of the genes. Deletions with a size in the range of introns (<400kb) are flagged as `read-through`, because there is a high chance that the fusion arises from [read-through transcription](interpretation-of-results.md#frequent-types-of-false-positives) rather than an underlying genomic deletion. Intragenic duplications with both breakpoints at splice-sites are flagged as `non-canonical-splicing`, because the supporting reads might originate from [circular RNA](interpretation-of-results.md#frequent-types-of-false-positives), which are very abundant even in normal tissue, but manifest as duplications in RNA-Seq data.

`direction1` and `direction2`
: These columns indicate the orientation of the fusion. A value of `downstream` means that the partner is fused downstream of the breakpoint, i.e. at a coordinate higher than the breakpoint. A value of `upstream` means the partner is fused at a coordinate lower than the breakpoint. When the prediction of strands or 5'- and 3'-genes fails, this information gives insight into which parts of the fused genes are retained in the fusion.

`split_reads1` and `split_reads2`
: The number of suppporting split reads (fragments) with an anchor in gene1 or gene2, respectively, is given in these columns. The gene to which the longer segment of the split read aligns is defined as the anchor.

`discordant_mates`
: This column contains the number of pairs (fragments) of discordant mates (a.k.a. spanning reads or bridge reads) supporting the fusion.

`confidence`
: Each prediction is assigned one of the confidences `low`, `medium`, or `high`. Several characteristics are taken into account, including: the number of supporting reads, the balance of split reads and discordant mates, the distance between the breakpoints, the type of event, whether the breakpoints are intragenic or not, and whether there are other events which corroborate the prediction, e.g. multiple isoforms or balanced translocations. See section [Interpretation of results](interpretation-of-results.md) for further advice on judging the credibility of predictions.


`closest_genomic_breakpoint1` and `closest_genomic_breakpoint2`
: When a matched whole-genome sequencing sample is available, one can feed structural variant calls obtained therefrom into Arriba (see parameter `-d`). Arriba then considers this information during fusion calling, which improves the overall accuracy. These two columns contain the coordinates of the genomic breakpoints which are closest to the transcriptomic breakpoints given in the columns `breakpoint1` and `breakpoint2`.

`filters`
: This column lists the filters which removed one or more of the supporting reads. The section [Internal algorithm](internal-algorithm.md) describes all filters in detail. The number of filtered reads is given in parantheses after the name of the filter. The total number of supporting reads can be obtained by summing up the reads given in the columns `split_reads1`, `split_reads2`, `discordant_mates`, and `filters`. If a filter discarded the event as a whole (all reads), the number of filtered reads is missing.

`fusion_transcript`
: If the parameter `-T` is is set, Arriba puts the transcript sequence in this column. The sequence is assembled from the supporting reads of the most highly expressed transcript. The column contains a `.`, when the sequence could not be predicted. This is the case, when the strands or the 5' end of the transcript could not be predicted reliably. The breakpoint is represented via a pipe symbol (`|`). When non-template bases are inserted between the fused genes, these bases are represented as lowercase letters between two pipes. Reference mismatches (SNPs or SNVs) are indicated as lowercase letters, insertions as bases between brackets (`[` and `]`), deleted bases as one or more dashes (`-`), introns as three underscores (`___`), and missing information due to insufficient coverage as an ellipsis (`...`). Ambiguous positions, such as positions with diverse reference mismatches, are represented as `n`.

`read_identifiers`
: When the parameter `-I` is set, this column contains the names of the supporting reads separated by commas.

discarded_fusions.tsv
-----------------------

The file `discarded_fusions.tsv` (as specified by the parameter `-O`) contains all events that Arriba classified as an artifact or that are also observed in healthy tissue. It has the same format as the file `fusions.tsv`. This file may be useful, if one suspects that an event should be present, but was erroneously discarded by Arriba.

