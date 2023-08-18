Manual installation
-------------------

Arriba has only a single prerequisite: [STAR](https://github.com/alexdobin/STAR) (version >=2.7.10a recommended). Download and install the tool according to the developers' instructions and make it available in your `$PATH`. If you want to make use of Arriba's [visualization tools](visualization.md), a few additional software components [need to be installed](visualization.md#publication-quality-figures).

Compile the latest stable version of Arriba or use the precompiled binaries in the download file. **Note: You should not use `git clone` to download Arriba, because the git repository does not include the blacklist! Instead, download the latest tarball from the [releases page](https://github.com/suhrig/arriba/releases/) as shown here:**

```bash
wget https://github.com/suhrig/arriba/releases/download/v2.4.0/arriba_v2.4.0.tar.gz
tar -xzf arriba_v2.4.0.tar.gz
cd arriba_v2.4.0 && make # or use precompiled binaries
```

Arriba requires an assembly in FastA format, gene annotation in GTF format, and a STAR index built from the two. You can use your preferred assembly and annotation, as long as their coordinates are compatible with hg19/hs37d5/GRCh37 or hg38/GRCh38 or mm10/GRCm38 or mm39/GRCm39. If you use another assembly, then the coordinates in the blacklist will not match and the predictions will contain many false positives. GENCODE annotation is recommended over RefSeq due to more comprehensive annotation of immunoglobulin/T-cell receptor loci and splice sites, which improves sensitivity. If you do not already have the files and a STAR index, you can use the script `download_references.sh`. It downloads the files to the current working directory and builds a STAR index. Run the script without arguments to see a list of available files. Choose a file with the keyword `viral` if Arriba is supposed to detect viral integration sites. Note that this step requires ~45 GB of RAM and 8 cores (can be adjusted by setting the environment variable `THREADS`).

```bash
./download_references.sh hs37d5viral+GENCODE19
```

The download file contains a script `run_arriba.sh`, which demonstrates the usage of Arriba (see also section [Workflow](workflow.md#demo-script)). We recommend that you use this as a guide to integrate Arriba into your existing STAR-based RNA-Seq pipeline. When Arriba is integrated properly, fusion detection only adds a few minutes to the regular alignment workflow, since Arriba utilizes the alignments produced by STAR during a standard RNA-Seq workflow and does not require alignment solely for the sake of fusion detection.

Run the demo script with 8 threads. In case of single-end data, the second FastQ file is omitted.

```bash
./run_arriba.sh STAR_index_hs37d5viral_GENCODE19/ GENCODE19.gtf hs37d5viral.fa database/blacklist_hg19_hs37d5_GRCh37_v2.4.0.tsv.gz database/known_fusions_hg19_hs37d5_GRCh37_v2.4.0.tsv.gz database/protein_domains_hg19_hs37d5_GRCh37_v2.4.0.gff3 8 test/read1.fastq.gz test/read2.fastq.gz
```

Installation using Docker
-------------------------

Install [Docker](https://www.docker.com/) according to the developers' instructions.

Run the script `download_references.sh` inside the Docker container. It downloads the assembly and gene annotation to the directory `/path/to/references` and builds a STAR index. Run the script without arguments to see a list of available files. Choose a file with the keyword `viral` if Arriba is supposed to detect viral integration sites. Note that this step requires ~45 GB of RAM and 8 cores (can be adjusted by passing the parameter `--env=THREADS=...`).

```bash
docker run --rm -v /path/to/references:/references uhrigs/arriba:2.4.0 download_references.sh hs37d5viral+GENCODE19
```

Use the following Docker command to run Arriba from the container. Replace `/path/to/` with the path to the respective input file. Leave the paths after the colons unmodified - these are the paths inside the Docker container. In case of single-end data, the second FastQ file is omitted. Running Arriba requires ~45 GB of RAM and 8 cores (can be adjusted by passing the parameter `--env=THREADS=...`).

```bash
docker run --rm \
       -v /path/to/output:/output \
       -v /path/to/references:/references:ro \
       -v /path/to/read1.fastq.gz:/read1.fastq.gz:ro \
       -v /path/to/read2.fastq.gz:/read2.fastq.gz:ro \
       uhrigs/arriba:2.4.0 \
       arriba.sh
```

Installation using Singularity
------------------------------

Install [Singularity](https://www.sylabs.io/) according to the developers' instructions.

The Docker container is compatible with Singularity. If desired, it can be converted to a Singularity image, but executing the Docker container directly works, too. Run the script `download_references.sh` inside the Docker container/Singularity image. It downloads the assembly and gene annotation to the directory `/path/to/references` and builds a STAR index. Run the script without arguments to see a list of available files. Choose a file with the keyword `viral` if Arriba is supposed to detect viral integration sites. Note that this step requires ~45 GB of RAM and 8 cores (can be adjusted by setting the environment variable `SINGULARITYENV_THREADS`).

```bash
mkdir /path/to/references
singularity exec -B /path/to/references:/references docker://uhrigs/arriba:2.4.0 download_references.sh hs37d5viral+GENCODE19
```

Use the following Singularity command to run Arriba from the container. Replace `/path/to/` with the path to the respective input file. Leave the paths after the colons unmodified - these are the paths inside the Singularity container. In case of single-end data, the second FastQ file is omitted. Running Arriba requires ~45 GB of RAM and 8 cores (can be adjusted by setting the environment variable `SINGULARITYENV_THREADS`).

```bash
singularity exec \
       -B /path/to/output:/output \
       -B /path/to/references:/references:ro \
       -B /path/to/read1.fastq.gz:/read1.fastq.gz:ro \
       -B /path/to/read2.fastq.gz:/read2.fastq.gz:ro \
       docker://uhrigs/arriba:2.4.0 \
       arriba.sh
```

Installation using Bioconda
---------------------------

Install [Miniconda](https://conda.io/) according to the developers' instructions.

Install the `arriba` package:

```bash
conda install -c conda-forge -c bioconda arriba=2.4.0
```

Run the script `download_references.sh`, which is installed inside the conda environment. It downloads the assembly and gene annotation to the current working directory and builds a STAR index. Run the script without arguments to see a list of available files. Choose a file with the keyword `viral` if Arriba is supposed to detect viral integration sites. Note that this step requires ~45 GB of RAM and 8 cores (can be adjusted by setting the environment variable `THREADS`). Replace `$CONDA_PREFIX` with the path to your conda environment.

```bash
$CONDA_PREFIX/var/lib/arriba/download_references.sh hs37d5viral+GENCODE19
```

To process FastQ files, run the script `run_arriba.sh`, which is installed inside the conda environment. In case of single-end data, the second FastQ file is omitted. Replace `$CONDA_PREFIX` with the path to your conda environment.

```bash
ARRIBA_FILES=$CONDA_PREFIX/var/lib/arriba
run_arriba.sh STAR_index_hs37d5viral_GENCODE19/ GENCODE19.gtf hs37d5viral.fa $ARRIBA_FILES/blacklist_hg19_hs37d5_GRCh37_v2.4.0.tsv.gz $ARRIBA_FILES/known_fusions_hg19_hs37d5_GRCh37_v2.4.0.tsv.gz $ARRIBA_FILES/protein_domains_hg19_hs37d5_GRCh37_v2.4.0.gff3 8 $ARRIBA_FILES/read1.fastq.gz $ARRIBA_FILES/read2.fastq.gz
```

Output files
------------

Arriba creates two output files:

- `fusions.tsv`: This is the main output file of Arriba and contains a list of fusion candidates sorted from highest to lowest confidence. Low-confidence fusions should be considered speculative unless there is additional evidence, e.g., when the fusion is a recurrent event in the given entity. For cohort-analyses, you might want to ignore low-confidence fusions.

- `fusions.discarded.tsv`: This file lists all events which were classified as artifacts or which are observed in healthy tissue (e.g., read-through transcripts and non-canonically spliced transcripts).

For a detailed description of the columns, refer to section [Output files](output-files.md).
