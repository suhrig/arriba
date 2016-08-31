#include <iostream>
#include <unordered_map>
#include <vector>
#include <string>
#include <sstream>
#include <algorithm>
#include "common.hpp"
#include "annotation.hpp"
#include "options_ariba.hpp"
#include "read_chimeric_alignments.hpp"
#include "filter_multi_mappers.hpp"
#include "filter_uninteresting_contigs.hpp"
#include "filter_inconsistently_clipped.hpp"
#include "filter_homopolymer.hpp"
#include "filter_duplicates.hpp"
#include "filter_same_gene.hpp"
#include "fusions.hpp"
#include "filter_both_intronic.hpp"
#include "filter_proximal_read_through.hpp"
#include "filter_low_entropy.hpp"
#include "filter_promiscuous_genes.hpp"
#include "filter_min_support.hpp"
#include "recover_known_fusions.hpp"
#include "recover_both_spliced.hpp"
#include "filter_blacklisted_ranges.hpp"
#include "filter_end_to_end.hpp"
#include "filter_pcr_fusions.hpp"
#include "merge_adjacent_fusions.hpp"
#include "select_best.hpp"
#include "filter_short_anchor.hpp"
#include "filter_nonexpressed.hpp"
#include "filter_mismappers.hpp"
#include "output_fusions.hpp"

using namespace std;

unordered_map<string,filter_t> FILTERS({
	{"inconsistently_clipped", NULL},
	{"homopolymer", NULL},
	{"duplicates", NULL},
	{"same_gene", NULL},
	{"intronic", NULL},
	{"read_through", NULL},
	{"mismappers", NULL},
	{"promiscuous_genes", NULL},
	{"min_support", NULL},
	{"known_fusions", NULL},
	{"spliced", NULL},
	{"blacklist", NULL},
	{"end_to_end", NULL},
	{"pcr_fusions", NULL},
	{"merge_adjacent", NULL},
	{"select_best", NULL},
	{"short_anchor", NULL},
	{"non_expressed", NULL},
	{"uninteresting_contigs", NULL},
	{"low_entropy", NULL}
});


int main(int argc, char **argv) {

	// initialize filter names
	for (auto i = FILTERS.begin(); i != FILTERS.end(); ++i)
		i->second = &i->first; // filters are represented by pointers to the name of the filter (this saves memory compared to storing strings)

	// parse command-line options
	options_t options = parse_arguments(argc, argv);

	// load chimeric alignments
	chimeric_alignments_t chimeric_alignments;
	contigs_t contigs;
	cout << "Reading chimeric alignments from '" << options.chimeric_bam_file << "'" << flush;
	cout << " (total=" << read_chimeric_alignments(options.chimeric_bam_file, chimeric_alignments, contigs) << ")" << endl;

	// load read-through alignments
	if (options.read_through_bam_file != "") {
		cout << "Reading read-through alignments from '" << options.read_through_bam_file << "'" << flush;
		cout << " (total=" << read_chimeric_alignments(options.read_through_bam_file, chimeric_alignments, contigs, true) << ")" << endl;
	}

	if (chimeric_alignments.size() == 0) {
		cerr << "ERROR: empty input files" << endl;
		exit(1);
	}

	// map contig IDs to names
	vector<string> contigs_by_id(contigs.size());
	for (contigs_t::iterator i = contigs.begin(); i != contigs.end(); ++i)
		contigs_by_id[i->second] = i->first;

	cout << "Filtering multi-mappers and single mates" << flush;
	cout << " (remaining=" << filter_multi_mappers(chimeric_alignments) << ")" << endl;

	vector<bool> interesting_contigs(contigs.size());
	if (options.filters.at("uninteresting_contigs") && !options.interesting_contigs.empty()) {
		istringstream iss(options.interesting_contigs);
		while (iss) {
			string contig;
			iss >> contig;
			contig = removeChr(contig);
			if (contigs.find(contig) != contigs.end())
				interesting_contigs[contigs[contig]] = true;
		}
		cout << "Filtering mates which do not map to interesting contigs (" << options.interesting_contigs << ")" << flush;
		cout << " (remaining=" << filter_uninteresting_contigs(chimeric_alignments, interesting_contigs) << ")" << endl;
	} else { // all contigs are interesting
		for (vector<bool>::iterator i = interesting_contigs.begin(); i != interesting_contigs.end(); ++i)
			*i = true;
	}

	if (options.filters.at("inconsistently_clipped")) {
		cout << "Filtering inconsistently clipped mates" << flush;
		cout << " (remaining=" << filter_inconsistently_clipped_mates(chimeric_alignments) << ")" << endl;
	}

	if (options.filters.at("homopolymer")) {
		cout << "Filtering breakpoints adjacent to homopolymers >=" << options.homopolymer_length << "nt" << flush;
		cout << " (remaining=" << filter_homopolymer(chimeric_alignments, options.homopolymer_length) << ")" << endl;
	}

	if (options.filters.at("duplicates")) {
		cout << "Filtering duplicates" << flush;
		cout << " (remaining=" << filter_duplicates(chimeric_alignments) << ")" << endl;
	}

	cout << "Annotating alignments using genes from '" << options.gene_annotation_file << "'" << endl << flush;
	// load GTF file
	annotation_t exon_annotation;
	annotation_t gene_annotation;
	read_annotation_gtf(options.gene_annotation_file, gene_annotation, exon_annotation, contigs);

	// translate exons to genes
	unordered_map<string,gene_t> genes;
	for (unsigned int i = 0; i < gene_annotation.size(); ++i)
		genes[gene_annotation[i].name] = i;
	for (annotation_t::iterator i = exon_annotation.begin(); i != exon_annotation.end(); ++i) {
		auto gene = genes.find(i->name);
		if (gene == genes.end()) {
			cerr << "ERROR: exon belongs to unknown gene: " << i->name << endl;
			exit(1);
		} else {
			i->id = gene->second;
		}
	}

	// sort genes and exons by coordinate (make index)
	annotation_index_t exon_annotation_index;
	make_annotation_index(exon_annotation, exon_annotation_index, contigs);
	annotation_index_t gene_annotation_index;
	make_annotation_index(gene_annotation, gene_annotation_index, contigs);

	// calculate sum of the lengths of all exons for each gene
	// we will need this to normalize the number of events over the gene length
	for (annotation_index_t::iterator contig = exon_annotation_index.begin(); contig != exon_annotation_index.end(); ++contig) {
		position_t region_start;
		for (contig_annotation_index_t::iterator region = contig->begin(); region != contig->end(); ++region) {
			for (gene_set_t::iterator overlapping_gene = region->second.begin(); overlapping_gene != region->second.end(); overlapping_gene = region->second.upper_bound(*overlapping_gene))
				gene_annotation[*overlapping_gene].exonic_length += region->first - region_start;
			region_start = region->first;
		}
	}

	// first, try to annotate with exons
	for (chimeric_alignments_t::iterator i = chimeric_alignments.begin(); i != chimeric_alignments.end(); ++i) {
		for (mates_t::iterator j = i->second.begin(); j != i->second.end(); ++j) {
			get_annotation_by_coordinate(j->contig, j->start, j->end, j->genes, exon_annotation_index);
			j->exonic = !j->genes.empty();
		}
	}

	// if the alignment does not map to an exon, try to map it to a gene
	for (chimeric_alignments_t::iterator i = chimeric_alignments.begin(); i != chimeric_alignments.end(); ++i) {
		for (mates_t::iterator j = i->second.begin(); j != i->second.end(); ++j) {
			if (j->genes.empty())
				get_annotation_by_coordinate(j->contig, j->start, j->end, j->genes, gene_annotation_index);
		}
		// try to resolve ambiguous mappings using mapping information from mate
		if (i->second.size() == 3) {
			gene_set_t combined;
			combine_annotations(i->second[SPLIT_READ].genes, i->second[MATE1].genes, combined);
			if (i->second[SPLIT_READ].genes.empty() || combined.size() < i->second[SPLIT_READ].genes.size())
				i->second[SPLIT_READ].genes = combined;
		}
	}

	// if the alignment maps neither to an exon nor to a gene, make a dummy gene which subsumes all alignments with a distance of 10kb
	annotation_t unmapped_alignments;
	for (chimeric_alignments_t::iterator i = chimeric_alignments.begin(); i != chimeric_alignments.end(); ++i) {
		for (mates_t::iterator j = i->second.begin(); j != i->second.end(); ++j) {
			if (j->genes.empty()) { // put all unmapped alignments in a annotation_t structure for sorting by coordinate
				annotation_record_t annotation_record = { 0, "", "", j->contig, j->start, j->end, FORWARD, 0, false };
				unmapped_alignments.push_back(annotation_record);
			}
		}
	}
	if (unmapped_alignments.size() > 0) {
		sort(unmapped_alignments.begin(), unmapped_alignments.end());
		annotation_record_t annotation_record = { 0, "", "", unmapped_alignments[0].contig, unmapped_alignments[0].start, unmapped_alignments[0].end, FORWARD, 10000 /*TODO more exact estimation of exonic_length*/, true };
		for (int i = 1; i <= unmapped_alignments.size(); ++i) {
			// subsume all unmapped alignments in a range of 10kb into a dummy gene with the generic name "contig:start-end"
			if (i == unmapped_alignments.size() || annotation_record.end+10000 < unmapped_alignments[i].start || unmapped_alignments[i].contig != annotation_record.contig) {
				annotation_record.name = contigs_by_id[annotation_record.contig] + ":" + to_string(annotation_record.start) + "-" + to_string(annotation_record.end);
				annotation_record.id = gene_annotation.size();
				gene_annotation.push_back(annotation_record);
				annotation_record.contig = unmapped_alignments[i].contig;
				annotation_record.start = unmapped_alignments[i].start;
			}
			annotation_record.end = unmapped_alignments[i].end;	
		}
	}

	// map yet unmapped alignments to the newly created dummy genes
	gene_annotation_index.clear();
	make_annotation_index(gene_annotation, gene_annotation_index, contigs); // index needs to be regenerated after adding dummy genes
	for (chimeric_alignments_t::iterator i = chimeric_alignments.begin(); i != chimeric_alignments.end(); ++i) {
		for (mates_t::iterator j = i->second.begin(); j != i->second.end(); ++j) {
			if (j->genes.empty())
				get_annotation_by_coordinate(j->contig, j->start, j->end, j->genes, gene_annotation_index);
		}
	}
	
	if (options.filters.at("same_gene")) {
		cout << "Filtering fragments with both mates in the same gene" << flush;
		cout << " (remaining=" << filter_same_gene(chimeric_alignments) << ")" << endl;
	}

	if (options.filters.at("read_through")) {
		cout << "Filtering intergenic read-through fragments with a distance <=" << options.min_read_through_distance << "bp" << flush;
		cout << " (remaining=" << filter_proximal_read_through(chimeric_alignments, gene_annotation, options.min_read_through_distance) << ")" << endl;
	}

	if (options.filters.at("low_entropy")) {
		cout << "Filtering reads with low entropy (k-mer content >=" << (options.max_kmer_content*100) << "%)" << flush;
		cout << " (remaining=" << filter_low_entropy(chimeric_alignments, 3, options.max_kmer_content) << ")" << endl;
	}

	cout << "Finding fusions and counting supporting reads" << flush;
	fusions_t fusions; // list of detected (preliminary) fusions
	cout << " (total=" << find_fusions(chimeric_alignments, fusions, exon_annotation_index) << ")" << endl;

	unsigned long int mapped_reads = 0;
	if (!options.low_tumor_content) {
		cout << "Counting mapped reads on interesting contigs" << flush;
		mapped_reads = count_mapped_reads(options.rna_bam_file, interesting_contigs);
		cout << " (total=" << mapped_reads << ")" << endl;
	}

	cout << "Estimating expected number of fusions by random chance (e-value)" << endl << flush;
	estimate_expected_fusions(fusions, gene_annotation, mapped_reads);

	if (options.filters.at("intronic")) {
		cout << "Filtering fusions with both breakpoints in intronic/intergenic regions" << flush;
		cout << " (remaining=" << filter_both_intronic(fusions) << ")" << endl;
	}

	if (options.filters.at("promiscuous_genes")) {
		cout << "Filtering fusions with an e-value >=" << options.evalue_cutoff << flush;
		cout << " (remaining=" << filter_promiscuous_genes(fusions, options.evalue_cutoff) << ")" << endl;
	}

	// this step must come after the 'promiscuous_genes' filter,
	// because fusions with low support are used to find promiscuous genes
	if (options.filters.at("min_support")) {
		cout << "Filtering fusions with <" << options.min_support << " supporting reads" << flush;
		cout << " (remaining=" << filter_min_support(fusions, options.min_support) << ")" << endl;
	}

	// this step must come right after the 'promiscuous_genes' and 'min_support' filters
	if (!options.known_fusions_file.empty() && options.filters.at("known_fusions")) {
		cout << "Searching for known fusions in '" << options.known_fusions_file << "'" << flush;
		cout << " (remaining=" << recover_known_fusions(fusions, options.known_fusions_file, genes, options.low_tumor_content) << ")" << endl;
	}

	// this step must come right after the 'promiscuous_genes' and 'min_support' filters
	if (options.filters.at("spliced")) {
		cout << "Searching for fusions with spliced split reads" << flush;
		cout << " (remaining=" << recover_both_spliced(fusions, gene_annotation, options.low_tumor_content) << ")" << endl;
	}

	if (options.filters.at("end_to_end")) {
		cout << "Filtering end-to-end fusions with low support" << flush;
		cout << " (remaining=" << filter_end_to_end_fusions(fusions, gene_annotation) << ")" << endl;
	}

	if (options.filters.at("merge_adjacent")) {
		cout << "Merging adjacent fusion breakpoints" << flush;
		cout << " (remaining=" << merge_adjacent_fusions(fusions, 2) << ")" << endl;
	}

	// this step must come after the 'merge_adjacent' filter,
	// or else adjacent breakpoints will be counted several times
	if (options.filters.at("pcr_fusions")) {
		cout << "Filtering PCR fusions" << flush;
		cout << " (remaining=" << filter_pcr_fusions(fusions, gene_annotation, 20, 4, 4, 3) << ")" << endl;
	}

	// this step must come after the 'merge_adjacent' filter,
	// because merging might yield a different best breakpoint
	if (options.filters.at("select_best")) {
		cout << "Selecting best breakpoints from genes with multiple breakpoints" << flush;
		cout << " (remaining=" << select_most_supported_breakpoints(fusions) << ")" << endl;
	}

	if (options.filters.at("short_anchor")) {
		cout << "Filtering fusions with anchors <=" << options.min_anchor_length << "nt" << flush;
		cout << " (remaining=" << filter_short_anchor(fusions, options.min_anchor_length) << ")" << endl;
	}

	// this step must come after the 'select_best' filter, because this filter removes breakpoints which are
	// only supported by discordant mates, if there is a fusion with breakpoints also supported by split reads
	if (options.filters.at("blacklist") && !options.blacklist_file.empty()) {
		cout << "Filtering blacklisted fusions in '" << options.blacklist_file << "'" << flush;
		cout << " (remaining=" << filter_blacklisted_ranges(fusions, options.blacklist_file, contigs, genes) << ")" << endl;
	}

	if (!options.assembly_file.empty()) {
		cout << "Fetching sequences of genes from '" << options.assembly_file << "'" << endl << flush;
		fetch_gene_sequences_from_fasta(options.assembly_file, fusions, gene_annotation, contigs_by_id);

		// this step must come near the end, because it is computationally expensive
		if (options.filters.at("mismappers")) {
			cout << "Re-aligning chimeric reads to filter fusions with >=" << (options.max_mismapper_fraction*100) << "% mis-mappers" << flush;
			cout << " (remaining=" << filter_mismappers(fusions, gene_annotation, options.max_mismapper_fraction) << ")" << endl;
		}
	}

	// this step must come near the end, because random BAM file accesses are slow
	if (options.filters.at("non_expressed")) {
		cout << "Filtering fusions with no expression in '" << options.rna_bam_file << "'" << flush;
		cout << " (remaining=" << filter_nonexpressed(fusions, options.rna_bam_file, chimeric_alignments, gene_annotation, exon_annotation) << ")" << endl;
	}

	cout << "Writing fusions to file '" << options.output_file << "'" << endl;
	write_fusions_to_file(fusions, options.output_file, gene_annotation, gene_annotation_index, exon_annotation_index, contigs_by_id, options.print_supporting_reads, options.print_fusion_sequence, false);

	if (options.discarded_output_file != "") {
		cout << "Writing discarded fusions to file '" << options.discarded_output_file << "'" << endl;
		write_fusions_to_file(fusions, options.discarded_output_file, gene_annotation, gene_annotation_index, exon_annotation_index, contigs_by_id, options.print_supporting_reads_for_discarded_fusions, options.print_fusion_sequence_for_discarded_fusions, true);
	}

	return 0;
}
