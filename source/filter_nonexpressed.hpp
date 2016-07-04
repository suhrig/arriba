#ifndef _FILTER_NONEXPRESSED_H
#define _FILTER_NONEXPRESSED_H 1

#include <string>
#include "common.hpp"

using namespace std;

unsigned int filter_nonexpressed(fusions_t& fusions, const string& bam_file_path, chimeric_alignments_t& chimeric_alignments);

#endif /* _FILTER_NONEXPRESSED_H */
