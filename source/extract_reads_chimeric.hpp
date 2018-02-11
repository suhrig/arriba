#ifndef _EXTRACT_READS_CHIMERIC_H
#define _EXTRACT_READS_CHIMERIC_H 1

#include "sam.h"

using namespace std;

void extract_chimeric(BGZF* chimeric_file, const bam1_t* const read1, const bam1_t* const read2);

#endif /* _EXTRACT_READS_CHIMERIC_H */
