#include <algorithm>
#include <cmath>
#include <iostream>
#include <map>
#include <vector>
#include <string>
#include <sstream>
#include <set>
#include <unordered_map>
#include "sam.h"
#include "common.hpp"
#include "annotation.hpp"
#include "read_compressed_file.hpp"
#include "htslib/faidx.h"

using namespace std;

bool parse_gtf_features(string gtf_features_string, gtf_features_t& gtf_features) {
	replace(gtf_features_string.begin(), gtf_features_string.end(), ',', ' ');
	replace(gtf_features_string.begin(), gtf_features_string.end(), '=', ' ');
	istringstream iss(gtf_features_string);
	while (iss) {
		string feature, value;
		iss >> feature >> value;
		if (feature != "" && value == "")
			return false;
		if (feature == "gene_name") {
			gtf_features.gene_name = value;
		} else if (feature == "gene_id") {
			gtf_features.gene_id = value;
		} else if (feature == "transcript_id") {
			gtf_features.transcript_id = value;
		} else if (feature == "gene_status") {
			gtf_features.gene_status = value;
		} else if (feature == "status_KNOWN") {
			gtf_features.keyword_known = value;
		} else if (feature == "feature_exon") {
			gtf_features.feature_exon = value;
		} else if (feature == "feature_UTR") {
			gtf_features.feature_utr = value;
		} else if (feature == "feature_gene") {
			gtf_features.feature_gene = value;
		} else if (feature == "") {
			// the last feature has been processed
		} else {
			return false;
		}
	}

	return true;
}


string removeChr(string contig) {
	if (contig.substr(0, 3) == "chr")
		contig = contig.substr(3);
	if (contig == "M")
		contig = "MT";
	return contig;
}

string addChr(string contig) {
	if (contig.substr(0, 3) != "chr")
		contig = "chr" + contig;
	if (contig == "chrMT")
		contig = "chrM";
	return contig;
}

bool get_gtf_attribute(const string& attributes, const string& attribute_name, string& attribute_value) {

	// find start of attribute
	size_t start = attributes.find(attribute_name + " \"");
	if (start != string::npos)
		start = attributes.find('"', start);
	if (start == string::npos) {
		cerr << "WARNING: failed to extract " << attribute_name << " from line in GTF file: " << attributes << endl;
		return false;
	}
	start++;

	// find end of attribute
	size_t end = attributes.find('"', start);
	if (end == string::npos) {
		cerr << "WARNING: failed to extract " << attribute_name << " from line in GTF file: " << attributes << endl;
		return false;
	}
	attribute_value = attributes.substr(start, end - start);

	return true;
}

void read_annotation_gtf(const string& filename, const contigs_t& contigs, const string& gtf_features_string, gene_annotation_t& gene_annotation, exon_annotation_t& exon_annotation, unordered_map<string,gene_t>& gene_names) {

	gtf_features_t gtf_features;
	parse_gtf_features(gtf_features_string, gtf_features);

	unordered_map<string,transcript_t> transcripts; // translates transcript IDs to numeric IDs
	unordered_map<string,gene_t> gene_by_gene_id; // maps gene IDs to genes (used to map exons to genes)
	unordered_map<exon_t,string> gene_id_by_exon; // maps exons to gene IDs (used to map exons to genes)

	stringstream gtf_file;
	autodecompress_file(filename, gtf_file);
	string line;
	while (getline(gtf_file, line)) {
		if (!line.empty() && line[0] != '#') { // skip comment lines

			istringstream iss(line);
			annotation_record_t annotation_record;
			string contig, strand, feature, attributes, trash, gene_name, gene_id;

			// parse line
			iss >> contig >> trash >> feature >> annotation_record.start >> annotation_record.end >> trash >> strand >> trash;
			if (contig.empty() || feature.empty() || strand.empty()) {
				cerr << "WARNING: failed to parse line in GTF file: " << line << endl;
				continue;
			}
			getline(iss, attributes);

			// extract gene name and ID from attributes
			if (!get_gtf_attribute(attributes, gtf_features.gene_name, gene_name) ||
			    !get_gtf_attribute(attributes, gtf_features.gene_id, gene_id))
				continue;

			contig = removeChr(contig);
			if (contigs.find(contig) == contigs.end()) {
				cerr << "WARNING: unknown contig in GTF file: " << contig << endl;

			} else {

				// make annotation record
				annotation_record.contig = contigs.at(contig);
				annotation_record.start--; // GTF files are one-based
				annotation_record.end--; // GTF files are one-based
				annotation_record.strand = (strand[0] == '+') ? FORWARD : REVERSE;

				if (feature == gtf_features.feature_gene) {

					// make gene annotation record
					gene_annotation_record_t gene_annotation_record;
					gene_annotation_record.copy(annotation_record);

					gene_annotation_record.name = gene_name;
					gene_annotation_record.exonic_length = 0; // is calculated later in ariba.cpp
					gene_annotation_record.is_dummy = false;
					string gene_status;
					if (!get_gtf_attribute(attributes, gtf_features.gene_status, gene_status))
						continue;
					gene_annotation_record.is_known = gene_status == gtf_features.keyword_known;
					gene_annotation.push_back(gene_annotation_record);

					// remember ID of gene, so we can map exons to genes later
					gene_by_gene_id[gene_id] = &(gene_annotation[gene_annotation.size()-1]);

				} else if (feature == gtf_features.feature_exon || feature == gtf_features.feature_utr) {

					// make exon annotation record
					exon_annotation_record_t exon_annotation_record;
					exon_annotation_record.copy(annotation_record);
					exon_annotation_record.is_transcript_start = false; // is set further down
					exon_annotation_record.is_transcript_end = false; // is set further down

					// extract transcript ID from attributes
					string transcript_id;
					if (!get_gtf_attribute(attributes, gtf_features.transcript_id, transcript_id))
						continue;
					exon_annotation_record.transcript = transcripts[transcript_id];
					if (exon_annotation_record.transcript == 0) // this is the first time we encounter this transcript ID => give it a numeric ID
						exon_annotation_record.transcript = transcripts[transcript_id] = transcripts.size();

					// give the record an ID (= a consecutive number)
					exon_annotation.push_back(exon_annotation_record);

					// remember ID of gene, so we can map exons to genes later
					gene_id_by_exon[&(exon_annotation[exon_annotation.size()-1])] = gene_id;
				}
			}
		}
	}

	// make a map of gene_name -> gene
	//TODO this can cause collisions, because gene names are not unique
	for (gene_annotation_t::iterator gene = gene_annotation.begin(); gene != gene_annotation.end(); ++gene)
		gene_names[gene->name] = &(*gene);

	// assign exons to genes
	for (exon_annotation_t::iterator exon = exon_annotation.begin(); exon != exon_annotation.end();) {
		auto gene_id = gene_by_gene_id.find(gene_id_by_exon[&(*exon)]);
		if (gene_id == gene_by_gene_id.end()) {
			cerr << "WARNING: exon belongs to unknown gene with ID: " << gene_id_by_exon[&(*exon)] << endl;
			exon = exon_annotation.erase(exon);
		} else {
			exon->gene = gene_id->second;
			++exon;
		}
	}

	// mark first/last exons of a transcript as transcript_start/transcript_end
	unordered_map< transcript_t,vector<exon_annotation_record_t*> > first_exons_by_transcript;
	unordered_map< transcript_t,vector<exon_annotation_record_t*> > last_exons_by_transcript;
	for (exon_annotation_t::iterator exon = exon_annotation.begin(); exon != exon_annotation.end(); ++exon) {

		// for each transcript, find the first exon
		auto old_first_exons = &(first_exons_by_transcript[exon->transcript]);
		if (old_first_exons->empty() ||
		    exon->strand == FORWARD && (**old_first_exons->begin()).start > exon->start ||
		    exon->strand == REVERSE && (**old_first_exons->begin()).end   < exon->end) {
			old_first_exons->clear();
			old_first_exons->push_back(&(*exon));
		} else if (exon->strand == FORWARD && (**old_first_exons->begin()).start == exon->start ||
		           exon->strand == REVERSE && (**old_first_exons->begin()).end   == exon->end) {
			old_first_exons->push_back(&(*exon));
		}

		// for each transcript, find the last exon
		auto old_last_exons = &(last_exons_by_transcript[exon->transcript]);
		if (old_last_exons->empty() ||
		    exon->strand == FORWARD && (**old_last_exons->begin()).end   < exon->end ||
		    exon->strand == REVERSE && (**old_last_exons->begin()).start > exon->start) {
			old_last_exons->clear();
			old_last_exons->push_back(&(*exon));
		} else if (exon->strand == FORWARD && (**old_last_exons->begin()).end   == exon->end ||
		           exon->strand == REVERSE && (**old_last_exons->begin()).start == exon->start) {
			old_last_exons->push_back(&(*exon));
		}

	}

	// mark all first exons
	for (auto first_exons = first_exons_by_transcript.begin(); first_exons != first_exons_by_transcript.end(); ++first_exons)
		for (auto first_exon = first_exons->second.begin(); first_exon != first_exons->second.end(); ++first_exon)
			(**first_exon).is_transcript_start = true;

	// mark all last exons
	for (auto last_exons = last_exons_by_transcript.begin(); last_exons != last_exons_by_transcript.end(); ++last_exons)
		for (auto last_exon = last_exons->second.begin(); last_exon != last_exons->second.end(); ++last_exon)
			(**last_exon).is_transcript_end = true;

	// don't mark exons of genes with just one exon as transcript_start/transcript_end (i.e., undo marking)
	// these are often predicted genes, which are often involved in splicing
	for (exon_annotation_t::iterator exon = exon_annotation.begin(); exon != exon_annotation.end(); ++exon)
		if (exon->is_transcript_start && exon->is_transcript_end)
			exon->is_transcript_start = exon->is_transcript_end = false;
}

// given a coordinate, return all exons, if the coordinate is a splice-site
void get_exons_from_splice_site(const gene_t gene, const direction_t direction, const contig_t contig, const position_t breakpoint, const exon_annotation_index_t& exon_annotation_index, exon_multiset_t& exons) {

	// nothing to do, if there are no exons on the given contig
	if (exon_annotation_index[contig].empty()) {
		exons.clear();
		return;
	}

	// find exons at the breakpoint
	exon_contig_annotation_index_t::const_iterator exons_at_breakpoint = exon_annotation_index[contig].lower_bound(breakpoint);

	// find exons before the breakpoint
	exon_contig_annotation_index_t::const_iterator exons_before_breakpoint = exons_at_breakpoint;
	if (exons_before_breakpoint != exon_annotation_index[contig].begin())
		--exons_before_breakpoint;

	if (breakpoint - exons_before_breakpoint->first <= MAX_SPLICE_SITE_DISTANCE || exons_at_breakpoint == exon_annotation_index[contig].begin()) {

		// get all exons at the breakpoint that are not a transcript start and belong to <gene>
		exon_multiset_t exons_of_gene_at_breakpoint;
		if (exons_at_breakpoint != exon_annotation_index[contig].end())
			for (exon_multiset_t::const_iterator exon = exons_at_breakpoint->second.begin(); exon != exons_at_breakpoint->second.end(); ++exon)
				if ((**exon).gene == gene && // exon must belong to given gene
				    !(((**exon).is_transcript_start && (**exon).strand == FORWARD || (**exon).is_transcript_end && (**exon).strand == REVERSE) && (**exon).start > exons_before_breakpoint->first && direction == UPSTREAM)) // when orientation is upstream, exon must not be terminal
					exons_of_gene_at_breakpoint.insert(*exon);

		// get all exons before the breakpoint that are not a transcript end and belong to <gene>
		exon_multiset_t exons_of_gene_before_breakpoint;
		for (exon_multiset_t::const_iterator exon = exons_before_breakpoint->second.begin(); exon != exons_before_breakpoint->second.end(); ++exon)
			if ((**exon).gene == gene && // exon must belong to given gene
			    !(((**exon).is_transcript_end && (**exon).strand == FORWARD || (**exon).is_transcript_start && (**exon).strand == REVERSE) && (**exon).end <= exons_before_breakpoint->first && direction == DOWNSTREAM)) // when orientation is downstream, exon must not be terminal
				exons_of_gene_before_breakpoint.insert(*exon);

		// we are at a splice-site, if the number of exons before vs. at the breakpoint differ
		if (direction == UPSTREAM && exons_of_gene_at_breakpoint.size() > exons_of_gene_before_breakpoint.size()) {
			exons = exons_of_gene_at_breakpoint;
			return;
		}
		if (direction == DOWNSTREAM && exons_of_gene_before_breakpoint.size() > exons_of_gene_at_breakpoint.size()) {
			exons = exons_of_gene_before_breakpoint;
			return;
		}
	}

	if (exons_at_breakpoint->first - breakpoint <= MAX_SPLICE_SITE_DISTANCE && exons_at_breakpoint != exon_annotation_index[contig].end()) {

		// get all exons at the breakpoint that are not a transcript end and belong to <gene>
		exon_multiset_t exons_of_gene_at_breakpoint;
		for (exon_multiset_t::const_iterator exon = exons_at_breakpoint->second.begin(); exon != exons_at_breakpoint->second.end(); ++exon)
			if ((**exon).gene == gene && // exon must belong to given gene
			    !(((**exon).is_transcript_end && (**exon).strand == FORWARD || (**exon).is_transcript_start && (**exon).strand == REVERSE) && (**exon).end <= exons_at_breakpoint->first && direction == DOWNSTREAM)) // when orientation is downstream, exon must not be terminal
				exons_of_gene_at_breakpoint.insert(*exon);

		// get all exons after the breakpoint that are not a transcript start and belong to <gene>
		exon_contig_annotation_index_t::const_iterator exons_after_breakpoint = exons_at_breakpoint;
		++exons_after_breakpoint;
		exon_multiset_t exons_of_gene_after_breakpoint;
		if (exons_after_breakpoint != exon_annotation_index[contig].end())
			for (exon_multiset_t::const_iterator exon = exons_after_breakpoint->second.begin(); exon != exons_after_breakpoint->second.end(); ++exon)
				if ((**exon).gene == gene && // exon must belong to given gene
				    !(((**exon).is_transcript_start && (**exon).strand == FORWARD || (**exon).is_transcript_end && (**exon).strand == REVERSE) && (**exon).start > exons_at_breakpoint->first && direction == UPSTREAM)) // when orientation is upstream, exon must not be terminal
					exons_of_gene_after_breakpoint.insert(*exon);

		// we are at a splice-site, if the number of exons after vs. at the breakpoint differ
		if (direction == UPSTREAM && exons_of_gene_after_breakpoint.size() > exons_of_gene_at_breakpoint.size()) {
			exons = exons_of_gene_after_breakpoint;
			return;
		}
		if (direction == DOWNSTREAM && exons_of_gene_at_breakpoint.size() > exons_of_gene_after_breakpoint.size()) {
			exons = exons_of_gene_at_breakpoint;
			return;
		}
	}

	exons.clear(); // if we get here, the breakpoint is not at a splice-site
}

// check if a breakpoint is near an annotated splice site
bool is_breakpoint_spliced(const gene_t gene, const direction_t direction, const contig_t contig, const position_t breakpoint, const exon_annotation_index_t& exon_annotation_index) {
	exon_multiset_t exons;
	get_exons_from_splice_site(gene, direction, contig, breakpoint, exon_annotation_index, exons);
	return !exons.empty();
}

void annotate_alignment(alignment_t& alignment, gene_set_t& gene_set, const exon_annotation_index_t& exon_annotation_index) {

	// first, try to annotate based on the boundaries (start+end) of the alignment
	exon_set_t exon_set;
	get_annotation_by_coordinate(alignment.contig, alignment.start, alignment.end, exon_set, exon_annotation_index);

	// translate exons to genes
	for (auto exon = exon_set.begin(); exon != exon_set.end(); ++exon)
		gene_set.insert((**exon).gene);

	// try to resolve ambiguity and strand by looking for splice sites that are specific for one gene
	if (alignment.cigar.size() > 1 && // if there are no CIGAR operations, there is nothing we can do
	    (gene_set.size() > 1 || alignment.predicted_strand_ambiguous)) { // if the strand and mapped genes are clear, there is nothing we need to do

		// cycle through CIGAR string and look for introns
		gene_set_t gene_set_supported_by_splicing;
		position_t reference_position = alignment.start;
		for (unsigned int i = 0; i < alignment.cigar.size() && gene_set_supported_by_splicing.empty(); ++i) {

			switch (alignment.cigar.operation(i)) {
				case BAM_CSOFT_CLIP:
				case BAM_CHARD_CLIP:
				case BAM_CREF_SKIP:

					// whenever we hit an intron (or clipped segment) in the CIGAR string, check if it aligns with splice sites for each gene
					gene_set_supported_by_splicing = gene_set;
					for (gene_set_t::iterator gene = gene_set_supported_by_splicing.begin(); gene != gene_set_supported_by_splicing.end();) {
						if (((alignment.cigar.operation(i) == BAM_CSOFT_CLIP || alignment.cigar.operation(i) == BAM_CHARD_CLIP) &&
						     (i == 0 && !is_breakpoint_spliced(*gene, UPSTREAM, alignment.contig, reference_position, exon_annotation_index) || // preclipped segment aligns with exon start
						      i != 0 && !is_breakpoint_spliced(*gene, DOWNSTREAM, alignment.contig, reference_position, exon_annotation_index)) || // postclipped segment aligns with exon end
						     alignment.cigar.operation(i) == BAM_CREF_SKIP &&
						     !is_breakpoint_spliced(*gene, DOWNSTREAM, alignment.contig, reference_position, exon_annotation_index) && // intron aligns with exon start
						     !is_breakpoint_spliced(*gene, UPSTREAM, alignment.contig, reference_position + alignment.cigar.op_length(i), exon_annotation_index))) // intron aligns with exon end
							gene = gene_set_supported_by_splicing.erase(gene);
						else
							++gene;
					}
			}

			switch (alignment.cigar.operation(i)) {
				case BAM_CREF_SKIP:
				case BAM_CMATCH:
				case BAM_CDEL:
					reference_position += alignment.cigar.op_length(i);
			}
		}

		if (!gene_set_supported_by_splicing.empty()) {

			// if none of the genes match the splice pattern of the alignment, return all genes,
			// otherwise, only return the matching genes
			if (gene_set_supported_by_splicing.size() < gene_set.size())
				gene_set = gene_set_supported_by_splicing;

			// try to get the strand of the alignment from the strand of the gene,
			// if this is unstranded data
			if (alignment.predicted_strand_ambiguous) { // this is unstranded data
				strand_t predicted_strand = (**gene_set_supported_by_splicing.begin()).strand;
				alignment.predicted_strand_ambiguous = false;
				for (gene_set_t::iterator gene = gene_set_supported_by_splicing.begin(); gene != gene_set_supported_by_splicing.end() && !alignment.predicted_strand_ambiguous; ++gene)
					if ((**gene).strand != predicted_strand)
						alignment.predicted_strand_ambiguous = true;

				if (!alignment.predicted_strand_ambiguous)
					alignment.predicted_strand = predicted_strand;
			}
		}
	}

}

void annotate_alignments(mates_t& mates, const exon_annotation_index_t& exon_annotation_index) {

	// annotate each mate individually
	for (mates_t::iterator mate = mates.begin(); mate != mates.end(); ++mate) {
		annotate_alignment(*mate, mate->genes, exon_annotation_index);
		mate->exonic = !mate->genes.empty();
	}

	// try to resolve ambiguous strand of one mate by infering from other mate
	if (mates[MATE1].predicted_strand_ambiguous && !mates[MATE2].predicted_strand_ambiguous) { // infer strand of MATE1 from MATE2
		mates[MATE1].predicted_strand = complement_strand_if(mates[MATE2].predicted_strand, mates[MATE1].strand == mates[MATE2].strand);
		mates[MATE1].predicted_strand_ambiguous = false;
	} else if (!mates[MATE1].predicted_strand_ambiguous && mates[MATE2].predicted_strand_ambiguous) { // infer strand of MATE2 from MATE1
		mates[MATE2].predicted_strand = complement_strand_if(mates[MATE1].predicted_strand, mates[MATE1].strand == mates[MATE2].strand);
		mates[MATE2].predicted_strand_ambiguous = false;
	} else if (!mates[MATE1].predicted_strand_ambiguous && !mates[MATE2].predicted_strand_ambiguous) { // both mates claim they know from which strand they come
		if ((mates[MATE1].predicted_strand != mates[MATE2].predicted_strand) != /*xor*/ (mates[MATE1].strand == mates[MATE2].strand)) { // contradiction => set to ambiguous
			mates[MATE1].predicted_strand_ambiguous = true;
			mates[MATE2].predicted_strand_ambiguous = true;
		}
	}

	if (mates.size() == 3) { // split read

		// try to resolve ambiguous mappings using mapping information from mate
		gene_set_t combined;
		combine_annotations(mates[SPLIT_READ].genes, mates[MATE1].genes, combined);
		if (mates[MATE1].genes.empty() || combined.size() < mates[MATE1].genes.size())
			mates[MATE1].genes = combined;
		if (mates[SPLIT_READ].genes.empty() || combined.size() < mates[SPLIT_READ].genes.size())
			mates[SPLIT_READ].genes = combined;

		// try to resolve ambiguous strand of split read or supplementary by infering from each other
		if (mates[SPLIT_READ].predicted_strand_ambiguous && !mates[SUPPLEMENTARY].predicted_strand_ambiguous) { // infer MATE1 and SPLIT_READ from SUPPLEMENTARY
			mates[MATE1].predicted_strand = complement_strand_if(mates[SUPPLEMENTARY].predicted_strand, mates[SUPPLEMENTARY].strand != mates[SPLIT_READ].strand);
			mates[MATE1].predicted_strand_ambiguous = false;
			mates[SPLIT_READ].predicted_strand = mates[MATE1].predicted_strand;
			mates[SPLIT_READ].predicted_strand_ambiguous = false;
		} else if (!mates[SPLIT_READ].predicted_strand_ambiguous && mates[SUPPLEMENTARY].predicted_strand_ambiguous) { // infer SUPPLEMENTARY from SPLIT_READ
			mates[SUPPLEMENTARY].predicted_strand = complement_strand_if(mates[SPLIT_READ].predicted_strand, mates[SUPPLEMENTARY].strand != mates[SPLIT_READ].strand);
			mates[SUPPLEMENTARY].predicted_strand_ambiguous = false;
		} else if (!mates[SPLIT_READ].predicted_strand_ambiguous && !mates[SUPPLEMENTARY].predicted_strand_ambiguous) { // both mates claim they know from which strand they come
			if ((mates[SPLIT_READ].predicted_strand != mates[SUPPLEMENTARY].predicted_strand) != /*xor*/ (mates[SPLIT_READ].strand != mates[SUPPLEMENTARY].strand)) { // contradiction => set to ambiguous
				mates[MATE1].predicted_strand_ambiguous = true;
				mates[SPLIT_READ].predicted_strand_ambiguous = true;
				mates[SUPPLEMENTARY].predicted_strand_ambiguous = true;
			}
		}

	}
}

// when a read overlaps with multiple genes, this function returns the boundaries of the biggest one
void get_boundaries_of_biggest_gene(gene_set_t& genes, position_t& start, position_t& end) {
	start = -1;
	end = -1;
	for (gene_set_t::iterator i = genes.begin(); i != genes.end(); ++i) {
		if (start == -1 || start > (**i).start)
			start = (**i).start;
		if (end == -1 || end < (**i).end)
			end = (**i).end;
	}
}

void fetch_gene_sequences_from_fasta(const string& assembly_file_path, fusions_t& fusions, const vector<string>& contigs_by_id) {

	// open FastA file and index
	faidx_t* genome_fa_index = fai_load(assembly_file_path.c_str());

	// fetch sequences of non-discarded fusions from FastA file
	for (fusions_t::iterator i = fusions.begin(); i != fusions.end(); ++i) {

		if (i->second.filter != NULL)
			continue; // skip discarded fusions

		vector<gene_t> genes = { i->second.gene1, i->second.gene2 };
		for (vector<gene_t>::iterator gene = genes.begin(); gene != genes.end(); ++gene) {

			if ((**gene).sequence.empty()) { // if we have not already fetched the sequence
				char* sequence;
				int length;
				string region_string = contigs_by_id[(**gene).contig] + ":" + to_string((**gene).start+1) + "-" + to_string((**gene).end+1);
				sequence = fai_fetch(genome_fa_index, region_string.c_str(), &length);
				if (length == 0) { // fetching failed
					// maybe the FastA file has "chr" in contig names => try again
					sequence = fai_fetch(genome_fa_index, addChr(region_string).c_str(), &length);
					if (length == 0) { // fetching failed again => give up
						cerr << "WARNING: Failed to fetch sequence of feature '" << (**gene).name << "' from '" << assembly_file_path << "'.";
						return;
					}
				}
				(**gene).sequence = string(sequence, sequence + length);
				// convert to uppercase
				std::transform((**gene).sequence.begin(), (**gene).sequence.end(), (**gene).sequence.begin(), (int (*)(int))std::toupper);
				free(sequence);
			}
		}
	}

	// close FastA file and index
	fai_destroy(genome_fa_index);
}

void dna_to_reverse_complement(string& dna, string& reverse_complement) {
	if (!reverse_complement.empty())
		reverse_complement.clear();
	reverse_complement.reserve(dna.length());
	for (string::reverse_iterator i = dna.rbegin(); i != dna.rend(); ++i)
		if (*i == 'a') {
			reverse_complement += 't';
		} else if (*i == 't') {
			reverse_complement += 'a';
		} else if (*i == 'c') {
			reverse_complement += 'g';
		} else if (*i == 'g') {
			reverse_complement += 'c';
		} else if (*i == 'A') {
			reverse_complement += 'T';
		} else if (*i == 'T') {
			reverse_complement += 'A';
		} else if (*i == 'C') {
			reverse_complement += 'G';
		} else if (*i == 'G') {
			reverse_complement += 'C';
		} else if (*i == '[') {
			reverse_complement += ']';
		} else if (*i == ']') {
			reverse_complement += '[';
		} else
			reverse_complement += *i;
}

string dna_to_reverse_complement(string& dna) {
	string reverse_complement;
	dna_to_reverse_complement(dna, reverse_complement);
	return reverse_complement;
}


void get_unique_transcripts_from_exons(const exon_multiset_t& exons, set<transcript_t>& transcripts) {
	transcripts.clear();
	for (exon_multiset_t::const_iterator exon = exons.begin(); exon != exons.end(); ++exon)
		transcripts.insert((**exon).transcript);
}

gene_multiset_t get_genes_from_exons(const exon_multiset_t& exons) {
	gene_multiset_t result;
	for (exon_multiset_t::const_iterator exon = exons.begin(); exon != exons.end(); ++exon)
		result.insert((**exon).gene);
	return result;
}

// get the distance between two positions after splicing (i.e., ignoring introns)
int get_spliced_distance(const contig_t contig, position_t position1, position_t position2, direction_t direction1, direction_t direction2, const gene_t gene, const exon_annotation_index_t& exon_annotation_index) {

	int negate = 1;
	if (position1 > position2) {
		swap(position1, position2);
		swap(direction1, direction2);
		negate = -1;
	}

	// find exon/intron of position1
	exon_contig_annotation_index_t::const_iterator p1 = exon_annotation_index[contig].lower_bound(position1);
	// find exon/intron of position2
	exon_contig_annotation_index_t::const_iterator p2 = exon_annotation_index[contig].lower_bound(position2);

	if (p1 == p2) // position1 and position2 are in the same intron/exon => distance = difference
		return (position2 - position1) * negate;

	// the exon/intron where position1/position2 is located must only be counted partially
	// (i.e., from position1/position2 to next exon boundary)
	position_t distance = 0;
	exon_multiset_t exons_at_position1, exons_at_position2;
	get_exons_from_splice_site(gene, direction1, contig, position1, exon_annotation_index, exons_at_position1);
	if (exons_at_position1.size() == 0) { // position is not at a splice site
		if (p1 != exon_annotation_index[contig].end())
			exons_at_position1 = p1->second;
		// add distance from position1 to next higher exon boundary
		distance += p1->first - position1;
	}
	get_exons_from_splice_site(gene, direction2, contig, position2, exon_annotation_index, exons_at_position2);
	if (exons_at_position2.size() > 0) {
		--p2;
	} else {
		if (p2 != exon_annotation_index[contig].end())
			exons_at_position2 = p2->second;
		// add distance from position2 to next lower exon boundary
		distance += position2;
		--p2;
		distance -= p2->first;
	}

	// check if the positions are in exons which belong to the same transcript
	// if so, calculate the distance considering only the exons of this transcript
	// if not, calculate the distance considering all the exons between position1 and 2
	set<transcript_t> transcripts1, transcripts2;
	get_unique_transcripts_from_exons(exons_at_position1, transcripts1);
	get_unique_transcripts_from_exons(exons_at_position2, transcripts2);
	set<transcript_t> common_transcripts;
	set_intersection(transcripts1.begin(), transcripts1.end(), transcripts2.begin(), transcripts2.end(), inserter(common_transcripts, common_transcripts.begin()));
	unordered_map<transcript_t,int> distance_by_transcript;
	for (auto common_transcript = common_transcripts.begin(); common_transcript != common_transcripts.end(); ++common_transcript)
		distance_by_transcript[*common_transcript] = distance;

	// sum up exon sizes between position1 and position2
	while (p1 != p2) {
		position_t boundary = p1->first;
		++p1;
		unsigned int exons_of_gene = 0;
		for (exon_multiset_t::const_iterator exon = p1->second.begin(); exon != p1->second.end() && exons_of_gene == 0; ++exon)
			if ((**exon).gene == gene)
				exons_of_gene++;
		if (exons_of_gene > 0) {
			// calculate distance considering all exons
			distance += p1->first - boundary;
			// calculate distance considering only exons of transcripts common to position1 and 2
			set<transcript_t> transcripts_at_current_position;
			get_unique_transcripts_from_exons(p1->second, transcripts_at_current_position);
			for (auto transcript = transcripts_at_current_position.begin(); transcript != transcripts_at_current_position.end(); ++transcript) {
				auto d = distance_by_transcript.find(*transcript);
				if (d != distance_by_transcript.end())
					d->second += p1->first - boundary;
			}
		}
	}

	// select transcript with shortest distance
	for (auto d = distance_by_transcript.begin(); d != distance_by_transcript.end(); ++d)
		if (d->second < distance)
			distance = d->second;

	return distance * negate;
}

