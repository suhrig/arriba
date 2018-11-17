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
template <class T> void make_annotation_index(annotation_t<T>& annotation, annotation_index_t<T*>& annotation_index) {
	annotation_index.resize(annotation.size()); // create a contig_annotation_index_t for each contig
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
		for (typename contig_annotation_index_t<T*>::iterator annotation_set = annotation_index[feature->contig].lower_bound(feature->end); annotation_set->first >= feature->start; --annotation_set)
			annotation_set->second.insert(&(*feature));
	}
}

template <class T> void combine_annotations(const annotation_set_t<T>& genes1, const annotation_set_t<T>& genes2, annotation_set_t<T>& combined, bool make_union) {
	// when the two ends of a read map to different genes, the mapping is ambiguous
	// in this case, we try to resolve the ambiguity by taking the gene that both - start and end - overlap with
	set_intersection(genes1.begin(), genes1.end(), genes2.begin(), genes2.end(), back_inserter(combined));
	if (combined.empty() && make_union)
		set_union(genes1.begin(), genes1.end(), genes2.begin(), genes2.end(), back_inserter(combined));
}

template <class T> void get_annotation_by_coordinate(const contig_t contig, position_t start, position_t end, annotation_set_t<T>& annotation_set, const annotation_index_t<T>& annotation_index) {
	if ((unsigned int) contig >= annotation_index.size()) {
		annotation_set.clear(); // return empty set
		return;
	}

	if (start == end) {

		// get all features at position
		typename contig_annotation_index_t<T>::const_iterator position = annotation_index[contig].lower_bound(start);
		if (position != annotation_index[contig].end())
			annotation_set = position->second;
		else
			annotation_set.clear(); // return empty set

	} else {
		if (start > end)
			swap(start, end);

		// get all features at start (+ 2bp)
		annotation_set_t<T> result_start;
		typename contig_annotation_index_t<T>::const_iterator position_start = annotation_index[contig].lower_bound(start);
		if (position_start != annotation_index[contig].end()) {
			result_start = position_start->second;
			if (position_start->first - start <= 2) {
				++position_start;
				if (position_start != annotation_index[contig].end())
					result_start.insert(position_start->second.begin(), position_start->second.end());
			}
		}

		// get all features at end (- 2 bp)
		annotation_set_t<T> result_end;
		typename contig_annotation_index_t<T>::const_iterator position_end = annotation_index[contig].lower_bound(end);
		if (position_end != annotation_index[contig].end())
			result_end = position_end->second;
		if (position_end != annotation_index[contig].begin() && annotation_index[contig].size() > 0) {
			--position_end;
			if (end - position_end->first <= 2)
				result_end.insert(position_end->second.begin(), position_end->second.end());
		}

		// take intersection of genes at start and end
		combine_annotations(result_start, result_end, annotation_set);
	}
}

