Bootstrap: docker
From: ubuntu:xenial

%labels
Maintainer Sebastian Uhrig @ DKFZ

%post
# install utilities
apt-get update -y && apt-get install -y wget

# install samtools
apt-get install -y samtools

# install STAR
wget -q -O - https://github.com/alexdobin/STAR/archive/master.tar.gz | \
tar -x -z --strip-components=3 -C /usr/local/bin -f - STAR-master/bin/Linux_x86_64_static/STAR

# install arriba
URL=$(wget -q -O - https://api.github.com/repos/suhrig/arriba/releases/latest | sed -n -e 's/.*"browser_download_url":\s*"\([^"]*\)".*/\1/p') && \
wget -q -O - "$URL" | tar -xzf -

# install dependencies of draw_fusions.R
apt-get install -y r-base && \
Rscript -e 'install.packages("circlize", repos="http://cran.r-project.org"); source("https://bioconductor.org/biocLite.R"); biocLite(c("GenomicRanges", "GenomicAlignments"))'

echo '[ -z "$THREADS" ] && export THREADS=8' >> $SINGULARITY_ENVIRONMENT

# make script to download references
echo '#!/bin/bash
wget -q -O - ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz |
gunzip -c | sed -e "s/^chrM\t/MT\t/" -e "s/^chr//" > /references/gencode.v19.annotation.gtf
wget -q -O - ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz |
gunzip -c > /references/hs37d5.fa
mkdir /references/STAR_index_hs37d5_gencode19
STAR --runMode genomeGenerate --genomeDir /references/STAR_index_hs37d5_gencode19 --genomeFastaFiles /references/hs37d5.fa --sjdbGTFfile /references/gencode.v19.annotation.gtf --runThreadN $THREADS --sjdbOverhang 150' > /usr/local/bin/download_references.sh && \
chmod a+x /usr/local/bin/download_references.sh

# make script to run arriba
echo '#!/bin/bash
cd /output
/arriba*/run_arriba.sh /references/STAR_index_hs37d5_gencode19 /references/gencode.v19.annotation.gtf /references/hs37d5.fa /arriba*/database/blacklist_hg19_hs37d5_GRCh37_*.tsv.gz /read1.fastq.gz /read2.fastq.gz $THREADS' > /usr/local/bin/arriba.sh && \
chmod a+x /usr/local/bin/arriba.sh

# make script to run draw_fusions.R
echo '#!/bin/bash
PROTEIN_DOMAINS=/arriba*/database/protein_domains_hg19_hs37d5_GRCh37_*.gff3
CYTOBANDS=/arriba*/database/cytobands_hg19_hs37d5_GRCh37_*.tsv
Rscript /arriba*/draw_fusions.R --annotation=/references/gencode.v19.annotation.gtf --fusions=/fusions.tsv --output=/output/fusions.pdf --proteinDomains=$PROTEIN_DOMAINS --alignments=/Aligned.sortedByCoord.out.bam --cytobands=$CYTOBANDS --optimizeDomainColors=TRUE' > /usr/local/bin/draw_fusions.sh && \
chmod a+x /usr/local/bin/draw_fusions.sh

