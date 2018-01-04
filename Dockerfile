FROM ubuntu:xenial
MAINTAINER Sebastian Uhrig @ DKFZ

# install utilities
RUN apt-get update -y && apt-get install -y wget zip

# install samtools
RUN apt-get install -y samtools

# install STAR
RUN wget -q -O STAR-master.zip https://github.com/alexdobin/STAR/archive/master.zip && \
unzip -q STAR-master.zip && \
cp STAR-master/bin/Linux_x86_64_static/STAR /usr/local/bin/

# install arriba
RUN URL=$(wget -q -O - https://api.github.com/repos/suhrig/arriba/releases/latest | sed -n -e 's/.*"browser_download_url":\s*"\([^"]*\)".*/\1/p') && \
wget -q -O - "$URL" | tar -xzf -

ENV THREADS=8

# make script to download references
RUN echo '#!/bin/bash\n\
wget -q -O - ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz | \n\
gunzip -c | sed -e "s/^chrM\t/MT\t/" -e "s/^chr//" > /references/gencode.v19.annotation.gtf \n\
wget -q -O - ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz | \n\
gunzip -c > /references/hs37d5.fa \n\
mkdir /references/STAR_index_hs37d5_gencode19 \n\
STAR --runMode genomeGenerate --genomeDir /references/STAR_index_hs37d5_gencode19 --genomeFastaFiles /references/hs37d5.fa --sjdbGTFfile /references/gencode.v19.annotation.gtf --runThreadN $THREADS --sjdbOverhang 150' > /usr/local/bin/download_references.sh && \
chmod a+x /usr/local/bin/download_references.sh

# make script to run arriba
RUN echo '#!/bin/bash\n\
cd /output\n\
/arriba*/run_arriba.sh /STAR_index /annotation.gtf /assembly.fa /read1.fastq.gz /read2.fastq.gz $THREADS' > /usr/local/bin/arriba.sh && \
chmod a+x /usr/local/bin/arriba.sh

