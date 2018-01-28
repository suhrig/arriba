Manual installation
-------------------

Arriba requires [STAR](https://github.com/alexdobin/STAR) (version >=2.5.3a recommended) and [samtools](http://www.htslib.org/). Download and install the two tools according to the developers' instructions and make them available in your `$PATH`.

Build a STAR index from your genome assembly and annotation. The following commands use the hs37d5 assembly and GencodeV19 annotation. Currently, the only supported assemblies are hg19/hs37d5/GRCh37 and hg38/GRCh38. Support for mm10 is in development. If you use another assembly, then the coordinates in the blacklist will not match and the predictions will contain many false positives. There are no restrictions on the annotation, but Gencode is recommended over RefSeq due to more comprehensive annotation of splice-sites, which improves sensitivity.

```bash
# download and extract annotation
wget ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz
gunzip -c gencode.v19.annotation.gtf.gz | sed -e 's/^chrM\t/MT\t/' -e 's/^chr//' > gencode.v19.annotation.gtf
# download and extract assembly
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz
gunzip hs37d5.fa.gz
# build STAR index
mkdir STAR_index_hs37d5_gencode19
STAR --runMode genomeGenerate --genomeDir STAR_index_hs37d5_gencode19 \
     --genomeFastaFiles hs37d5.fa --sjdbGTFfile gencode.v19.annotation.gtf \
     --runThreadN 8 --sjdbOverhang 150
```

Compile the latest stable version of Arriba or use the precompiled binaries in the download file. **Note: You should not use `git clone` to download Arriba, because the git repository does not include the blacklist! Instead, download the latest tarball from the [releases page](https://github.com/suhrig/arriba/releases/) as shown here:**

```bash
wget https://github.com/suhrig/arriba/releases/download/v0.12.0/arriba_v0.12.0.tar.gz
tar -xzf arriba_v0.12.0.tar.gz
cd arriba_v0.12.0 && make && cd .. # or use precompiled binaries
```

The download file contains a script `run_arriba.sh`, which demonstrates the usage of Arriba (see also section [Execution](execution.md)). We recommend that you use this as a guide to integrate Arriba into your existing STAR-based RNA-Seq pipeline. When Arriba is integrated properly, fusion detection only adds a few minutes to the regular alignment workflow, since Arriba utilizes the alignments produced by STAR during a normal RNA-Seq workflow and does not require alignment solely for the sake of fusion detection.

Run the demo script with 8 threads:
```bash
arriba_v0.12.0/run_arriba.sh STAR_index_hs37d5_gencode19/ gencode.v19.annotation.gtf hs37d5.fa arriba_v0.12.0/database/blacklist_hg19_hs37d5_GRCh37_2018-01-13.tsv.gz read1.fastq.gz read2.fastq.gz 8
```

Installation using Docker
-------------------------

Install [Docker](https://www.docker.com/) according to the developers' instructions.

Build the Docker image:

```bash
docker build --tag arriba:latest https://raw.githubusercontent.com/suhrig/arriba/master/Dockerfile
```

If you have not already downloaded the annotation/assembly and built a STAR index, run the `download_references.sh` script inside the container. Note that this step requires ~30 GB of RAM and 8 cores (can be adjusted with `--env THREADS=...`). The script downloads the assembly hs37d5 and GencodeV19 annotation. Please refer to the manual installation instructions or modify the `Dockerfile`, if you wish to use a different assembly/annotation. The files will be extracted to the directory `/path/to/references` in the following example:

```bash
docker run --rm -t -v /path/to/references:/references arriba:latest download_references.sh
```

Use the following Docker command to run Arriba from the container. Replace `/path/to/` with the path to the respective input file. Leave the paths after the colons unmodified - these are the paths inside the Docker container.

```bash
docker run --rm -t \
       -v /path/to/output:/output \
       -v /path/to/references/STAR_index_hs37d5_gencode19:/STAR_index:ro \
       -v /path/to/references/gencode.v19.annotation.gtf:/annotation.gtf:ro \
       -v /path/to/references/hs37d5.fa:/assembly.fa:ro \
       -v /path/to/read1.fastq.gz:/read1.fastq.gz:ro \
       -v /path/to/read2.fastq.gz:/read2.fastq.gz:ro \
       arriba:latest \
       arriba.sh
```

Output files
------------

Arriba creates two output files:

- `fusions.tsv`: This is the main output file of Arriba and contains a list of fusion candidates sorted from highest to lowest confidence. Low-confidence fusions should be considered speculative unless there is additional evidence, e.g., when the fusion is a recurrent event in the given entity. For cohort-analyses, you might want to ignore low-confidence fusions.

- `discarded_fusions.tsv`: This file lists all events which were classified as artifacts or which are observed in healthy tissue (e.g., read-through transcripts and non-canonically spliced transcripts).

For a detailed description of the columns, refer to section [Output files](output-files.md).
