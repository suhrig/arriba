#include <algorithm>
#include <cmath>
#include <iostream>
#include <iterator>
#include <map>
#include <vector>
#include <string>
#include <sstream>
#include <set>
#include <tuple>
#include <unordered_map>
#include "sam.h"
#include "common.hpp"
#include "annotation.hpp"
#include "read_compressed_file.hpp"

using namespace std;

void split_string(const string& unsplit, const char separator, vector<string>& split) {
	stringstream ss(unsplit);
	string value;
	while (getline(ss, value, separator))
		if (!value.empty())
			split.push_back(value);
}

bool parse_gtf_features(string gtf_features_string, gtf_features_t& gtf_features) {
	replace(gtf_features_string.begin(), gtf_features_string.end(), ',', ' ');
	istringstream iss(gtf_features_string);
	while (iss) {
		string feature_value_pair;
		iss >> feature_value_pair;
		tsv_stream_t tsv(feature_value_pair, '=');
		string feature, value;
		tsv >> feature >> value;
		if (feature != "" && value == "")
			return false;
		if (feature == "gene_name") {
			split_string(value, '|', gtf_features.gene_name);
		} else if (feature == "gene_id") {
			split_string(value, '|', gtf_features.gene_id);
		} else if (feature == "transcript_id") {
			split_string(value, '|', gtf_features.transcript_id);
		} else if (feature == "feature_exon") {
			split_string(value, '|', gtf_features.feature_exon);
		} else if (feature == "feature_CDS") {
			split_string(value, '|', gtf_features.feature_cds);
		} else if (feature == "") {
			// the last feature has been processed
		} else {
			return false;
		}
	}

	return !gtf_features.gene_name.empty() &&
	       !gtf_features.gene_id.empty() &&
	       !gtf_features.transcript_id.empty() &&
	       !gtf_features.feature_exon.empty() &&
	       !gtf_features.feature_cds.empty();
}

void remove_gene(const gene_t gene_to_remove, gene_annotation_t& gene_annotation, exon_annotation_t& exon_annotation) {
	// remove all exons belonging to gene
	for (exon_annotation_t::iterator exon = exon_annotation.begin(); exon != exon_annotation.end();) {
		if (exon->gene == gene_to_remove) {
			exon = exon_annotation.erase(exon);
		} else {
			++exon;
		}
	}
	// remove gene
	for (gene_annotation_t::iterator gene = gene_annotation.begin(); gene != gene_annotation.end(); ++gene) {
		if (&(*gene) == gene_to_remove) {
			gene = gene_annotation.erase(gene);
			break;
		}
	}
}

void remove_transcript(const transcript_t transcript_to_remove, gene_annotation_t& gene_annotation, exon_annotation_t& exon_annotation) {
	gene_t gene = NULL;

	// remove all exons belonging to transcript
	for (exon_annotation_t::iterator exon = exon_annotation.begin(); exon != exon_annotation.end();) {
		if (exon->transcript == transcript_to_remove) {
			gene = exon->gene;
			exon = exon_annotation.erase(exon);
		} else {
			++exon;
		}
	}

	// shrink gene to boundaries of remaining transcripts
	position_t new_start = -1;
	position_t new_end = -1;
	for (exon_annotation_t::iterator exon = exon_annotation.begin(); exon != exon_annotation.end(); ++exon) {
		if (exon->gene == gene) {
			if (new_start == -1 || new_start > exon->start)
				new_start = exon->start;
			if (new_end == -1 || new_end < exon->end)
				new_end = exon->end;
		}
	}
	if (new_start == -1) {
		remove_gene(gene, gene_annotation, exon_annotation); // there are no exons left in the gene
	} else {
		gene->start = new_start;
		gene->end = new_end;
	}
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

bool get_gtf_attribute(const string& attributes, const vector<string>& attribute_names, string& attribute_value) {

	// find start of attribute
	size_t start = string::npos;
	for (auto attribute_name = attribute_names.begin(); attribute_name != attribute_names.end() && start == string::npos; ++attribute_name)
		start = attributes.find(*attribute_name + " \"");
	if (start != string::npos)
		start = attributes.find('"', start);
	if (start == string::npos) {
		cerr << "WARNING: failed to extract ";
		for (auto attribute_name = attribute_names.begin(); attribute_name != attribute_names.end(); ++attribute_name) {
			if (attribute_name != attribute_names.begin())
				cerr << "|";
			cerr << *attribute_name;
		}
		cerr << " from line in GTF file: " << attributes << endl;
		return false;
	}
	start++;

	// find end of attribute
	size_t end = attributes.find('"', start);
	if (end == string::npos) {
		cerr << "WARNING: failed to extract ";
		for (auto attribute_name = attribute_names.begin(); attribute_name != attribute_names.end(); ++attribute_name) {
			if (attribute_name != attribute_names.begin())
				cerr << "|";
			cerr << *attribute_name;
		}
		cerr << " from line in GTF file: " << attributes << endl;
		return false;
	}
	attribute_value = attributes.substr(start, end - start);

	return true;
}

bool sort_exons_by_coordinate(const exon_annotation_record_t* exon1, const exon_annotation_record_t* exon2) {
	return *exon1 < *exon2;
}

struct coding_region_t {
	position_t start, end;
	string transcript_id;
};

void read_annotation_gtf(const string& filename, const string& gtf_features_string, contigs_t& contigs, gene_annotation_t& gene_annotation, transcript_annotation_t& transcript_annotation, exon_annotation_t& exon_annotation, unordered_map<string,gene_t>& gene_names) {

	gtf_features_t gtf_features;
	parse_gtf_features(gtf_features_string, gtf_features);

	unordered_map<string,transcript_t> transcripts; // translates transcript IDs to numeric IDs
	unordered_map<tuple<string,contig_t,strand_t>,gene_t> gene_by_id; // maps gene IDs to genes (used to map exons to genes)
	unordered_map<string,vector<exon_annotation_record_t*> > exons_by_transcript_id; // maps transcript IDs to exons (used to map coding regions to exons)
	vector<coding_region_t> coding_regions; // keeps track of coding regions (used to map coding regions to exons)

	gene_set_t bogus_genes; // genes with bogus annotation are ignored

	autodecompress_file_t gtf_file(filename);
	string line;
	unsigned int new_id = 0; // ID generator for genes and transcripts
	while (gtf_file.getline(line)) {
		if (!line.empty() && line[0] != '#') { // skip comment lines

			tsv_stream_t tsv(line);
			annotation_record_t annotation_record;
			string contig, strand, feature, attributes, trash, gene_name, gene_id;

			// parse line
			tsv >> contig >> trash >> feature >> annotation_record.start >> annotation_record.end >> trash >> strand >> trash >> attributes;
			if (tsv.fail() || contig.empty() || feature.empty() || strand.empty()) {
				cerr << "WARNING: failed to parse line in GTF file: " << line << endl;
				continue;
			}

			// extract gene name and ID from attributes
			if (!get_gtf_attribute(attributes, gtf_features.gene_name, gene_name) ||
			    !get_gtf_attribute(attributes, gtf_features.gene_id, gene_id))
				continue;

			// if Gencode, remove version number from gene ID
			string::size_type trim_position = string::npos;
			if (gene_id.substr(0, 3) == "ENS" && (trim_position = gene_id.find_last_of('.', string::npos)) != string::npos)
				gene_id = gene_id.substr(0, trim_position);

			// convert string representation of contig to numeric ID
			contig = removeChr(contig);
			pair<contigs_t::iterator,bool> find_contig_by_name = contigs.insert(pair<string,contig_t>(contig, contigs.size())); // this adds a new contig only if it does not yet exist

			// make annotation record
			annotation_record.contig = find_contig_by_name.first->second;
			annotation_record.start--; // GTF files are one-based
			annotation_record.end--; // GTF files are one-based
			annotation_record.strand = (strand[0] == '+') ? FORWARD : REVERSE;

			if (find(gtf_features.feature_exon.begin(), gtf_features.feature_exon.end(), feature) != gtf_features.feature_exon.end()) {

				// make exon annotation record
				exon_annotation_record_t exon_annotation_record;
				exon_annotation_record.copy(annotation_record);
				exon_annotation_record.coding_region_start = -1;
				exon_annotation_record.coding_region_end = -1;

				// extract transcript ID from attributes
				string transcript_id;
				if (!get_gtf_attribute(attributes, gtf_features.transcript_id, transcript_id))
					continue;
				// if Gencode, remove version number from transcript ID
				string short_transcript_id = transcript_id;
				if (short_transcript_id.substr(0, 3) == "ENS" && (trim_position = short_transcript_id.find_last_of('.', string::npos)) != string::npos)
					short_transcript_id = short_transcript_id.substr(0, trim_position);
				exon_annotation_record.transcript = transcripts[short_transcript_id];
				if (exon_annotation_record.transcript == NULL) { // this is the first time we encounter this transcript ID => make a new transcript_annotation_record_t
					transcript_annotation_record_t transcript_annotation_record;
					transcript_annotation_record.id = new_id++;
					transcript_annotation_record.name = transcript_id;
					transcript_annotation_record.first_exon = NULL; // is set once we have loaded all exons
					transcript_annotation_record.last_exon = NULL; // is set once we have loaded all exons
					transcript_annotation.push_back(transcript_annotation_record);
					exon_annotation_record.transcript = transcripts[short_transcript_id] = &(*transcript_annotation.rbegin());
				}

				// make a gene annotation record, if this is the first exon of a gene
				gene_t gene = gene_by_id[make_tuple(gene_id, annotation_record.contig, annotation_record.strand)];
				if (gene == NULL) {
					gene_annotation_record_t gene_annotation_record;
					gene_annotation_record.copy(annotation_record);
					gene_annotation_record.name = gene_name;
					gene_annotation_record.id = new_id++;
					gene_annotation_record.exonic_length = 0; // is calculated later in arriba.cpp
					gene_annotation_record.is_dummy = false;
					gene_annotation_record.is_protein_coding = false;
					gene_annotation.push_back(gene_annotation_record);
					gene = &(*gene_annotation.rbegin());
					gene_by_id[make_tuple(gene_id, annotation_record.contig, annotation_record.strand)] = gene;
				} else { // gene has already been seen previously
					// expand the boundaries of the gene, so that all exons fit inside
					if (gene->start > exon_annotation_record.start)
						gene->start = exon_annotation_record.start;
					if (gene->end < exon_annotation_record.end)
						gene->end = exon_annotation_record.end;
					// check if annotation is sensible
					if (gene->contig != annotation_record.contig || gene->end - gene->start > 3000000) {
						cout << "WARNING: gene ID '" << gene_id << "' appears to be non-unique and will be ignored" << endl;
						bogus_genes.insert(gene);
					}
				}
				exon_annotation_record.gene = gene;

				exon_annotation.push_back(exon_annotation_record);

				// keep track of all exons of a transcript, so we can map coding regions to exons later
				exons_by_transcript_id[transcript_id].push_back(&(*exon_annotation.rbegin()));

			} else if (find(gtf_features.feature_cds.begin(), gtf_features.feature_cds.end(), feature) != gtf_features.feature_cds.end()) {

				// remember which regions of an exon are coding
				coding_region_t coding_region;
				coding_region.start = annotation_record.start;
				coding_region.end = annotation_record.end;
				if (!get_gtf_attribute(attributes, gtf_features.transcript_id, coding_region.transcript_id))
					continue;
				coding_regions.push_back(coding_region);
			}
		}
	}

	if (gene_annotation.empty()) {
		cerr << "ERROR: failed to parse GTF file, please consider using -G" << endl;
		exit(1);
	}

	// map coding regions to exons
	for (auto coding_region = coding_regions.begin(); coding_region != coding_regions.end(); ++coding_region) {

		// find exons of the same transcript as the coding region
		auto transcript = exons_by_transcript_id.find(coding_region->transcript_id);
		if (transcript == exons_by_transcript_id.end()) {
			cerr << "WARNING: CDS record has unknown transcript ID: " << coding_region->transcript_id << endl;
			continue;
		}

		// check if coding region overlaps with any of the exons of its transcript
		for (auto exon = transcript->second.begin(); exon != transcript->second.end(); ++exon)
			if ((**exon).start <= coding_region->start && (**exon).end >= coding_region->start || // exon overlaps start of coding region
			    (**exon).start <= coding_region->end && (**exon).end >= coding_region->end || // exon overlaps end of coding region
			    (**exon).start >= coding_region->start && (**exon).end <= coding_region->end) { // coding region encompasses exon
				(**exon).coding_region_start = max(coding_region->start, (**exon).start);
				(**exon).coding_region_end = min(coding_region->end, (**exon).end);
				(**exon).gene->is_protein_coding = true;
			}
	}

	// make double-linked list of exons according to their position in the transcript
	for (auto transcript = exons_by_transcript_id.begin(); transcript != exons_by_transcript_id.end(); ++transcript) {
		sort(transcript->second.begin(), transcript->second.end(), sort_exons_by_coordinate);
		for (auto exon = transcript->second.begin(); exon != transcript->second.end(); ++exon) {
			(**exon).previous_exon = (exon != transcript->second.begin()) ? *(exon-1) : NULL;
			(**exon).next_exon = (exon+1 != transcript->second.end()) ? *(exon+1) : NULL;
		}
	}

	// remove bogus genes
	for (gene_set_t::iterator bogus_gene = bogus_genes.begin(); bogus_gene != bogus_genes.end(); ++bogus_gene)
		remove_gene(*bogus_gene, gene_annotation, exon_annotation);

	// fix some errors in the Gencode annotation
	// remove fusion transcript FIP1L1:PDGFRA
	if (transcripts.find("ENST00000507166") != transcripts.end())
		remove_transcript(transcripts.at("ENST00000507166"), gene_annotation, exon_annotation);
	// remove fusion transcript GOPC:ROS1
	if (transcripts.find("ENST00000467125") != transcripts.end())
		remove_transcript(transcripts.at("ENST00000467125"), gene_annotation, exon_annotation);
	// remove fusion transcripts MTAP:CDKN2B-AS1
	if (transcripts.find("ENST00000404796") != transcripts.end())
		remove_transcript(transcripts.at("ENST00000404796"), gene_annotation, exon_annotation);
	if (transcripts.find("ENST00000577563") != transcripts.end())
		remove_transcript(transcripts.at("ENST00000577563"), gene_annotation, exon_annotation);
	if (transcripts.find("ENST00000580900") != transcripts.end())
		remove_transcript(transcripts.at("ENST00000580900"), gene_annotation, exon_annotation);

	// make a map of gene_name -> gene
	//TODO this can cause collisions, because gene names are not unique
	for (gene_annotation_t::iterator gene = gene_annotation.begin(); gene != gene_annotation.end(); ++gene)
		gene_names[gene->name] = &(*gene);

	// compute starts and ends of all transcripts
	for (exon_annotation_t::iterator exon = exon_annotation.begin(); exon != exon_annotation.end(); ++exon) {
		if (exon->transcript->first_exon == NULL || exon->start < exon->transcript->first_exon->start)
			exon->transcript->first_exon = &(*exon);
		if (exon->transcript->last_exon == NULL || exon->end > exon->transcript->last_exon->end)
			exon->transcript->last_exon = &(*exon);
	}

}

bool filter_exons_near_splice_site(const gene_t gene, const direction_t direction, const position_t breakpoint, const exon_set_t& exons_near_splice_site) {
	// return only exons which
	// - belong to given gene
	// - have a boundary within MAX_SPLICE_SITE_DISTANCE from the breakpoint
	// - are not the first/last exon in a transcript
	// - unless:
	//   - the transcript has only one exon
	//   - the gene misses a start/stop codon (=> indicates incomplete annotation)
	for (exon_set_t::const_iterator exon = exons_near_splice_site.begin(); exon != exons_near_splice_site.end(); ++exon)
		if ((**exon).gene == gene)
			if (direction == UPSTREAM &&
			    abs((**exon).start - breakpoint) <= MAX_SPLICE_SITE_DISTANCE &&
			    ((**exon).previous_exon != NULL || // exon is not a terminal one
			     (**exon).previous_exon == NULL && (**exon).next_exon == NULL || // unless transcript has only one exon
			     (**exon).start == (**exon).coding_region_start) || // or unless the first base of the exon is coding (=> gene is not annotated properly and misses preceeding exons (see TCR genes))
			    direction == DOWNSTREAM &&
			    abs((**exon).end - breakpoint) <= MAX_SPLICE_SITE_DISTANCE &&
			    ((**exon).next_exon != NULL || // exon is not a terminal one
			     (**exon).previous_exon == NULL && (**exon).next_exon == NULL || // unless transcript has only one exon
			     (**exon).end == (**exon).coding_region_end)) // or unless the last base of the exon is coding (=> gene is not annotated properly and misses following exons (see TCR genes))
				return true;
	return false;
}

// check if a breakpoint is near an annotated splice site
bool is_breakpoint_spliced(const gene_t gene, const direction_t direction, const position_t breakpoint, const exon_annotation_index_t& exon_annotation_index) {

	// nothing to do, if there are no exons on the given contig
	if ((unsigned int) gene->contig >= exon_annotation_index.size() || exon_annotation_index[gene->contig].empty())
		return false;

	// find exons in the vicinity of the breakpoint
	exon_contig_annotation_index_t::const_iterator exons_at_breakpoint = exon_annotation_index[gene->contig].lower_bound(breakpoint);
	if (exons_at_breakpoint != exon_annotation_index[gene->contig].end()) {
		if (filter_exons_near_splice_site(gene, direction, breakpoint, exons_at_breakpoint->second))
			return true;
		exon_contig_annotation_index_t::const_iterator exons_after_breakpoint = exons_at_breakpoint;
		++exons_after_breakpoint;
		if (exons_after_breakpoint != exon_annotation_index[gene->contig].end() &&
		    filter_exons_near_splice_site(gene, direction, breakpoint, exons_after_breakpoint->second))
			return true;
	}
	if (exons_at_breakpoint != exon_annotation_index[gene->contig].begin()) {
		exon_contig_annotation_index_t::const_iterator exons_before_breakpoint = exons_at_breakpoint;
		--exons_before_breakpoint;
		if (filter_exons_near_splice_site(gene, direction, breakpoint, exons_before_breakpoint->second))
			return true;
	}

	return false;
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
						     (i == 0 && !is_breakpoint_spliced(*gene, UPSTREAM, reference_position, exon_annotation_index) || // preclipped segment aligns with exon start
						      i != 0 && !is_breakpoint_spliced(*gene, DOWNSTREAM, reference_position, exon_annotation_index)) || // postclipped segment aligns with exon end
						     alignment.cigar.operation(i) == BAM_CREF_SKIP &&
						     !is_breakpoint_spliced(*gene, DOWNSTREAM, reference_position, exon_annotation_index) && // intron aligns with exon start
						     !is_breakpoint_spliced(*gene, UPSTREAM, reference_position + alignment.cigar.op_length(i), exon_annotation_index))) { // intron aligns with exon end
							gene = gene_set_supported_by_splicing.erase(gene);
						} else {
							++gene;
						}
					}
			}

			switch (alignment.cigar.operation(i)) {
				case BAM_CREF_SKIP:
				case BAM_CMATCH:
				case BAM_CDIFF:
				case BAM_CEQUAL:
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
	for (gene_set_t::iterator gene = genes.begin(); gene != genes.end(); ++gene) {
		if (start == -1 || start > (**gene).start)
			start = (**gene).start;
		if (end == -1 || end < (**gene).end)
			end = (**gene).end;
	}
}

// get the distance between two positions after splicing (i.e., ignoring introns)
int get_spliced_distance(const contig_t contig, position_t position1, position_t position2, direction_t direction1, direction_t direction2, const gene_t gene, const exon_annotation_index_t& exon_annotation_index) {

	// make sure position1 contains the smaller coordinate
	if (position1 > position2) {
		swap(position1, position2);
		swap(direction1, direction2);
	}

	// take the plain distance, if no exons are annotated for the given contig
	if ((unsigned int) contig >= exon_annotation_index.size() || exon_annotation_index[contig].empty())
		return position2 - position1;

	// find exon/intron of position1 & 2
	exon_contig_annotation_index_t::const_iterator p1 = exon_annotation_index[contig].lower_bound(position1);
	exon_contig_annotation_index_t::const_iterator p2 = exon_annotation_index[contig].lower_bound(position2);

	// take the plain distance, if no exons are annotated for the given region or if the positions are in the same exon
	if (p1 == exon_annotation_index[contig].end() || p1 == p2)
		return position2 - position1;

	int distance = 0;

	// check if position1 is at a splice-site
	// if not, add the distance from the position to the next exon boundary
	bool position1_is_spliced = direction1 == DOWNSTREAM && is_breakpoint_spliced(gene, direction1, position1, exon_annotation_index);
	if (!position1_is_spliced)
		distance += p1->first - position1;

	// move p1 towards p2 and sum up the lengths of all exons in-between for each transcript
	unordered_map<transcript_t,int> distance_by_transcript;
	unordered_map<transcript_t,position_t> boundary_by_transcript;
	position_t boundary = p1->first;
	for (; p1 != p2; p1++) {
		for (auto exon = p1->second.begin(); exon != p1->second.end(); exon++) {
			auto transcript = distance_by_transcript.insert(make_pair((**exon).transcript, distance)); // if this is the first exon of a transcript, initialize with distance incl. all exons
			transcript.first->second += p1->first - boundary;
			boundary_by_transcript[(**exon).transcript] = p1->first;
		}
		distance += p1->first - boundary; // measure distance between positions considering all exons of all transcripts
		boundary = p1->first;
	}

	// add the remaining distance from position2 to the next exon boundary
	bool position2_is_spliced = direction2 == UPSTREAM && is_breakpoint_spliced(gene, direction2, position2, exon_annotation_index);
	for (auto transcript = distance_by_transcript.begin(); transcript != distance_by_transcript.end(); ++transcript)
		if (!position2_is_spliced)
			transcript->second += position2 - boundary_by_transcript[transcript->first];
	if (!position2_is_spliced)
		distance += position2 - boundary;

	// find transcript with shortest spliced distance between positions
	for (auto transcript = distance_by_transcript.begin(); transcript != distance_by_transcript.end(); ++transcript)
		if (transcript->second < distance)
			distance = transcript->second;

	return distance;
}

