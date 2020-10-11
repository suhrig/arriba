FROM ubuntu:bionic
MAINTAINER Sebastian Uhrig @ DKFZ

# install dependencies
RUN export DEBIAN_FRONTEND=noninteractive && \
apt-get update -y && \
apt-get install -y --no-install-recommends build-essential samtools r-base wget ca-certificates libcurl4-openssl-dev libxml2-dev && \
Rscript -e 'install.packages("circlize", repos="http://cran.r-project.org"); install.packages("BiocManager"); BiocManager::install(c("GenomicRanges", "GenomicAlignments"))'

# install version of STAR that supports --chimMultimapNmax and --chimOutType WithinBAM
RUN wget -q -O - 'https://github.com/alexdobin/STAR/archive/2.7.6a.tar.gz' | tar -xzf - && \
cp -p STAR-2.7.6a/bin/Linux_x86_64/STAR /usr/local/bin/

# install arriba
RUN wget -q -O - "https://github.com/suhrig/arriba/releases/download/v1.2.0/arriba_v1.2.0.tar.gz" | tar -xzf -

# make wrapper script for download_references.sh
RUN echo '#!/bin/bash\n\
cd /references\n\
/arriba*/download_references.sh $1 && \\\n\
ASSEMBLY=$(sed -e "s/viral+.*//" -e "s/+.*//" <<<"$1") && \\\n\
cp /arriba*/database/*$ASSEMBLY* /references' > /usr/local/bin/download_references.sh && \
chmod a+x /usr/local/bin/download_references.sh

# make wrapper script for run_arriba.sh
RUN echo '#!/bin/bash\n\
cd /output && \\\n\
/arriba*/run_arriba.sh /references/STAR_index_* /references/*.gtf /references/*.fa /references/blacklist_*.tsv.gz /references/known_fusions_*.tsv.gz /references/protein_domains_*.gff3 ${THREADS-8} /read1.fastq.gz $(ls /read2.fastq.gz 2> /dev/null)' > /usr/local/bin/arriba.sh && \
chmod a+x /usr/local/bin/arriba.sh

# make wrapper script for draw_fusions.R
RUN echo '#!/bin/bash\n\
Rscript /arriba*/draw_fusions.R --annotation=$(ls /references/*.gtf) --fusions=/fusions.tsv --output=/output/fusions.pdf --proteinDomains=$(ls /references/protein_domains_*.gff3) --alignments=/Aligned.sortedByCoord.out.bam --cytobands=$(ls /references/cytobands_*.tsv)' > /usr/local/bin/draw_fusions.sh && \
chmod a+x /usr/local/bin/draw_fusions.sh

