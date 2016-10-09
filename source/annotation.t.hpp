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
#include "annotation.hpp"
#include "read_compressed_file.hpp"
#include "htslib/faidx.h"

using namespace std;

// split overlapping genes into disjunct regions
// for each region, collect the genes overlapping it and store them in a gene set
// for example, if the annotation contains two regions:
// - gene1 chr1:10,000-20,000
// - gene2 chr1:12,000-13,000
// then the resulting index will contain the following items:
// - chr1:10,000-11,999 gene1
// - chr1:12,000-13,000 gene1+gene2
// - chr1:13,001-20,000 gene1
template <class T> void make_annotation_index(annotation_t<T>& annotation, annotation_index_t<T*>& annotation_index, const contigs_t& contigs) {
	annotation_index.resize(contigs.size()); // create a contig_annotation_index_t for each contig
	for (typename annotation_t<T>::iterator feature = annotation.begin(); feature != annotation.end(); ++feature) {

		typename contig_annotation_index_t<T*>::const_iterator overlapping_features = annotation_index[feature->contig].lower_bound(feature->end);
		if (overlapping_features == annotation_index[feature->contig].end())
			annotation_index[feature->contig][feature->end]; // this creates an empty gene set, if it does not exist yet
		else
			annotation_index[feature->contig][feature->end] = overlapping_features->second;

		overlapping_features = annotation_index[feature->contig].lower_bound(feature->start-1);
		if (overlapping_features == annotation_index[feature->contig].end())
			annotation_index[feature->contig][feature->start-1]; // this creates an empty gene set, if it does not exist yet
		else
			annotation_index[feature->contig][feature->start-1] = overlapping_features->second;

		// add the gene to all gene sets between start and end of the gene
		for (typename contig_annotation_index_t<T*>::iterator i = annotation_index[feature->contig].lower_bound(feature->end); i->first >= feature->start; --i)
			i->second.insert(&(*feature));
	}
}

//TODO use const reference for annotation_multiset?
template <class T> void annotation_multiset_to_set(annotation_multiset_t<T> annotation_multiset, annotation_set_t<T>& annotation_set) {
	for (typename annotation_multiset_t<T>::iterator i = annotation_multiset.begin(); i != annotation_multiset.end(); i = annotation_multiset.upper_bound(*i))
		annotation_set.insert(*i);
}

//TODO use const reference for annotation_multiset?
template <class T> annotation_set_t<T> annotation_multiset_to_set(annotation_multiset_t<T> annotation_multiset) {
	annotation_set_t<T> annotation_set;
	annotation_multiset_to_set(annotation_multiset, annotation_set);
	return annotation_set;
}

template <class T> void combine_annotations(const annotation_set_t<T>& genes1, const annotation_set_t<T>& genes2, annotation_set_t<T>& combined, bool make_union) {
	// when the two ends of a read map to different genes, the mapping is ambiguous
	// in this case, we try to resolve the ambiguity by taking the gene that both - start and end - overlap with
	set_intersection(genes1.begin(), genes1.end(), genes2.begin(), genes2.end(), inserter(combined, combined.begin()));
	if (combined.empty() && make_union)
		set_union(genes1.begin(), genes1.end(), genes2.begin(), genes2.end(), inserter(combined, combined.begin()));
}

template <class T> void get_annotation_by_coordinate(const contig_t contig, const position_t start, const position_t end, annotation_set_t<T>& annotation_set, const annotation_index_t<T>& annotation_index) {
//TODO support strand-specific libraries
	if (contig < annotation_index.size()) {
		typename contig_annotation_index_t<T>::const_iterator result_start = annotation_index[contig].lower_bound(start);
		annotation_set_t<T> empty_set;
		if (start == end) {
			if (result_start != annotation_index[contig].end())
				annotation_multiset_to_set(result_start->second, annotation_set);
			else
				annotation_set = empty_set;
		} else {
			typename contig_annotation_index_t<T>::const_iterator result_end = annotation_index[contig].lower_bound(end);
			combine_annotations(
				(result_start != annotation_index[contig].end()) ? annotation_multiset_to_set(result_start->second) : empty_set,
				(result_end != annotation_index[contig].end()) ? annotation_multiset_to_set(result_end->second) : empty_set,
				annotation_set
			);
		}
	}
}

