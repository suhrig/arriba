#ifndef _READ_CHIMERIC_ALIGNMENTS_H
#define _READ_CHIMERIC_ALIGNMENTS_H 1

#include <string>
#include "common.hpp"

using namespace std;

unsigned int read_chimeric_alignments(const string& bam_file_path, chimeric_alignments_t& chimeric_alignments, contigs_t& contigs, const bool read_through = false);

void assign_strands_from_strandedness(chimeric_alignments_t& chimeric_alignments, const strandedness_t strandedness);

#endif /* _READ_CHIMERIC_ALIGNMENTS_H */

