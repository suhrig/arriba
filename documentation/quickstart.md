Manual installation
-------------------

Arriba requires [STAR](https://github.com/alexdobin/STAR) (version >=2.5.3a recommended). If you want to visualize the results, you will also need [samtools](http://www.htslib.org/) for sorting and indexing of BAM files. Download and install the two tools according to the developers' instructions and make them available in your `$PATH`.

Compile the latest stable version of Arriba or use the precompiled binaries in the download file. **Note: You should not use `git clone` to download Arriba, because the git repository does not include the blacklist! Instead, download the latest tarball from the [releases page](https://github.com/suhrig/arriba/releases/) as shown here:**

```bash
wget https://github.com/suhrig/arriba/releases/download/v0.11.0/arriba_v0.11.0.tar.gz
tar -xzf arriba_v0.11.0.tar.gz
cd arriba_v0.11.0 && make # or use precompiled binaries
```

Arriba requires an assembly in FastA format, gene annotation in GTF format, and a STAR index built from the two. You can use your preferred assembly and annotation, as long as their coordinates are compatible with hg19/hs37d5/GRCh37 or hg38/GRCh38. Support for mm10 is in development. If you use another assembly, then the coordinates in the blacklist will not match and the predictions will contain many false positives. Gencode annotation is recommended over RefSeq due to more comprehensive annotation of splice-sites, which improves sensitivity. If you do not already have the files and a STAR index, you can use the script `download_references.sh`. It downloads the files to the current working directory and builds a STAR index. Run the script without arguments to see a list of available files. Note that this step requires ~30 GB of RAM and 8 cores (or whatever number of cores you pass as the second argument).

```bash
./download_references.sh hs37d5+GENCODE19
```

The download file contains a script `run_arriba.sh`, which demonstrates the usage of Arriba (see also section [Execution](execution.md)). We recommend that you use this as a guide to integrate Arriba into your existing STAR-based RNA-Seq pipeline. When Arriba is integrated properly, fusion detection only adds a few minutes to the regular alignment workflow, since Arriba utilizes the alignments produced by STAR during a normal RNA-Seq workflow and does not require alignment solely for the sake of fusion detection.

Run the demo script with 8 threads:

```bash
./run_arriba.sh STAR_index_hs37d5_GENCODE19/ GENCODE19.gtf hs37d5.fa database/blacklist_hg19_hs37d5_GRCh37_2018-04-04.tsv.gz read1.fastq.gz read2.fastq.gz 8
```

Installation using Docker
-------------------------

Install [Docker](https://www.docker.com/) according to the developers' instructions.

Build the Docker image:

```bash
docker build --tag arriba:latest https://raw.githubusercontent.com/suhrig/arriba/master/Dockerfile
```

Run the script `download_references.sh` inside the container. It downloads the assembly and gene annotation to the directory `/path/to/references` and builds a STAR index. Run the script without arguments to see a list of available files. Note that this step requires ~30 GB of RAM and 8 cores (or whatever number of cores you pass as the second argument).

```bash
docker run --rm -v /path/to/references:/references arriba:latest download_references.sh hs37d5+GENCODE19
```

Use the following Docker command to run Arriba from the container. Replace `/path/to/` with the path to the respective input file. Leave the paths after the colons unmodified - these are the paths inside the Docker container.

```bash
docker run --rm \
       -v /path/to/output:/output \
       -v /path/to/references:/references:ro \
       -v /path/to/read1.fastq.gz:/read1.fastq.gz:ro \
       -v /path/to/read2.fastq.gz:/read2.fastq.gz:ro \
       arriba:latest \
       arriba.sh
```

Installation using Singularity
------------------------------

Install [Singularity](http://singularity.lbl.gov/) according to the developers' instructions.

Build the Singularity image:

```bash
wget https://raw.githubusercontent.com/suhrig/arriba/master/Singularityfile
singularity build arriba.img Singularityfile
```

Run the script `download_references.sh` inside the container. It downloads the assembly and gene annotation to the directory `/path/to/references` and builds a STAR index. Run the script without arguments to see a list of available files. Note that this step requires ~30 GB of RAM and 8 cores (or whatever number of cores you pass as the second argument).

```bash
mkdir /path/to/references
singularity exec -B /path/to/references:/references arriba.img download_references.sh hs37d5+GENCODE19
```

Use the following Singularity command to run Arriba from the container. Replace `/path/to/` with the path to the respective input file. Leave the paths after the colons unmodified - these are the paths inside the Singularity container.

```bash
singularity exec \
       -B /path/to/output:/output \
       -B /path/to/references:/references:ro \
       -B /path/to/read1.fastq.gz:/read1.fastq.gz:ro \
       -B /path/to/read2.fastq.gz:/read2.fastq.gz:ro \
       arriba.img \
       arriba.sh
```

Output files
------------

Arriba creates two output files:

- `fusions.tsv`: This is the main output file of Arriba and contains a list of fusion candidates sorted from highest to lowest confidence. Low-confidence fusions should be considered speculative unless there is additional evidence, e.g., when the fusion is a recurrent event in the given entity. For cohort-analyses, you might want to ignore low-confidence fusions.

- `fusions.discarded.tsv`: This file lists all events which were classified as artifacts or which are observed in healthy tissue (e.g., read-through transcripts and non-canonically spliced transcripts).

For a detailed description of the columns, refer to section [Output files](output-files.md).
