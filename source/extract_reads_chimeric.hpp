#ifndef _EXTRACT_READS_CHIMERIC_H
#define _EXTRACT_READS_CHIMERIC_H 1

#include "sam.h"

using namespace std;

void extract_chimeric(samFile* chimeric_file, bam_hdr_t* bam_header, const bam1_t* const read1, const bam1_t* const read2);

#endif /* _EXTRACT_READS_CHIMERIC_H */
