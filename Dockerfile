FROM ubuntu:focal
MAINTAINER Sebastian Uhrig @ DKFZ

# install dependencies
RUN export DEBIAN_FRONTEND=noninteractive && \
apt-get update -y && \
apt-get install -y --no-install-recommends wget ca-certificates samtools r-base r-cran-circlize r-bioc-genomicalignments r-bioc-genomicranges libxml2 && \
apt-get clean && rm -rf /var/lib/apt/lists/*

# install version of STAR that supports --chimMultimapNmax and --chimOutType WithinBAM
RUN wget -qO - 'https://github.com/alexdobin/STAR/archive/2.7.10a.tar.gz' | \
tar --strip-components=3 -C /usr/local/bin -xzf - 'STAR-2.7.10a/bin/Linux_x86_64/STAR'

# install arriba
RUN wget -qO - 'https://github.com/suhrig/arriba/releases/download/v2.3.0/arriba_v2.3.0.tar.gz' | tar -xzf - --exclude='arriba*/.git'

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

