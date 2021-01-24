About
=====

Arriba is a command-line tool for the detection of gene fusions from RNA-Seq data. It was developed for the use in a clinical research setting. Therefore, short runtimes and high sensitivity were important design criteria. It is based on the ultrafast [STAR aligner](https://github.com/alexdobin/STAR) and the post-alignment runtime is typically just ~2 minutes. In contrast to many other fusion detection tools which build on STAR, Arriba does not require to reduce the `alignIntronMax` parameter of STAR to detect fusions arising from focal deletions.

Apart from gene fusions, Arriba can detect other structural rearrangements with potential clinical relevance, such as viral integration sites, internal tandem duplications, whole exon duplications, truncations of genes (i.e., breakpoints in introns and intergenic regions).

Arriba is the winner of the [DREAM SMC-RNA Challenge](https://www.synapse.org/SMC_RNA), an international competition organized by ICGC, TCGA, IBM, and Sage Bionetworks to determine the current gold standard for the detection of gene fusions from RNA-Seq data. The final results of the challenge are posted on the [Round 5 Leaderboard](https://www.synapse.org/#!Synapse:syn2813589/wiki/588511).

Get help
========

Use the [GitHub issue tracker](https://github.com/suhrig/arriba/issues) to get help or to report bugs.

Citation
========

Sebastian Uhrig, Julia Ellermann, Tatjana Walther, Pauline Burkhardt, Martina Fröhlich, Barbara Hutter, Umut H. Toprak, Olaf Neumann, Albrecht Stenzinger, Claudia Scholl, Stefan Fröhling and Benedikt Brors: *Accurate and efficient detection of gene fusions from RNA sequencing data.* Genome Research. Published in Advance January 13, 2021. doi: [10.1101/gr.257246.119](https://doi.org/10.1101/gr.257246.119)

License
-------

The code, software and database files of Arriba are distributed under the MIT/Expat License, with the exception of the script `draw_fusions.R`, which is distributed under the GNU GPL v3 due to dependencies on GPL-licensed R packages. The terms and conditions of both licenses can be found in the [LICENSE file](https://raw.githubusercontent.com/suhrig/arriba/master/LICENSE).

User manual
===========

Please refer to the [user manual](http://arriba.readthedocs.io/en/latest/) for installation instructions and information about usage. **Note: You should not use `git clone` to download Arriba, because the git repository does not include the blacklist and other database files!**

1. [Quickstart](https://arriba.readthedocs.io/en/latest/quickstart/)

   - [Manual installation](https://arriba.readthedocs.io/en/latest/quickstart/#manual-installation)
   - [Installation using Docker](https://arriba.readthedocs.io/en/latest/quickstart/#installation-using-docker)
   - [Installation using Singularity](https://arriba.readthedocs.io/en/latest/quickstart/#installation-using-singularity)
   - [Installation using Bioconda](https://arriba.readthedocs.io/en/latest/quickstart/#installation-using-bioconda)
   - [Output files](https://arriba.readthedocs.io/en/latest/quickstart/#output-files)

2. [Workflow](https://arriba.readthedocs.io/en/latest/workflow/)

   - [Demo script](https://arriba.readthedocs.io/en/latest/workflow/#demo-script)

3. [Input files](https://arriba.readthedocs.io/en/latest/input-files/)

   - [Alignments](https://arriba.readthedocs.io/en/latest/input-files/#alignments)
   - [Assembly](https://arriba.readthedocs.io/en/latest/input-files/#assembly)
   - [Annotation](https://arriba.readthedocs.io/en/latest/input-files/#annotation)
   - [Blacklist](https://arriba.readthedocs.io/en/latest/input-files/#blacklist)
   - [Known fusions](https://arriba.readthedocs.io/en/latest/input-files/#known-fusions)
   - [Tags](https://arriba.readthedocs.io/en/latest/input-files/#tags)
   - [Protein domains](https://arriba.readthedocs.io/en/latest/input-files/#protein-domains)
   - [Structural variant calls from WGS](https://arriba.readthedocs.io/en/latest/input-files/#structural-variant-calls-from-wgs)

4. [Output files](https://arriba.readthedocs.io/en/latest/output-files/)
   
   - [fusions.tsv](https://arriba.readthedocs.io/en/latest/output-files/#fusionstsv)
   - [fusions.discarded.tsv](https://arriba.readthedocs.io/en/latest/output-files/#fusionsdiscardedtsv)

5. [Visualization](https://arriba.readthedocs.io/en/latest/visualization/)

   - [Publication-quality figures](https://arriba.readthedocs.io/en/latest/visualization/#publication-quality-figures)
   - [Inspection of events using IGV](https://arriba.readthedocs.io/en/latest/visualization/#inspection-of-events-using-igv)

6. [Command line options](https://arriba.readthedocs.io/en/latest/command-line-options/)

   - [Arriba](https://arriba.readthedocs.io/en/latest/command-line-options/#arriba)
   - [draw_fusions.R](https://arriba.readthedocs.io/en/latest/command-line-options/#draw_fusionsr)

7. [Interpretation of results](https://arriba.readthedocs.io/en/latest/interpretation-of-results/)

   - [Confidence scoring](https://arriba.readthedocs.io/en/latest/interpretation-of-results/#confidence-scoring)
   - [Activating fusions](https://arriba.readthedocs.io/en/latest/interpretation-of-results/#activating-fusions)
   - [Inactivating fusions](https://arriba.readthedocs.io/en/latest/interpretation-of-results/#inactivating-fusions)
   - [Supporting read count](https://arriba.readthedocs.io/en/latest/interpretation-of-results/#supporting-read-count)
   - [Frequent types of false positives](https://arriba.readthedocs.io/en/latest/interpretation-of-results/#frequent-types-of-false-positives)
   - [Multiple transcript variants](https://arriba.readthedocs.io/en/latest/interpretation-of-results/#multiple-transcript-variants)
   - [Cohort analysis](https://arriba.readthedocs.io/en/latest/interpretation-of-results/#cohort-analysis)

8. [Current limitations](https://arriba.readthedocs.io/en/latest/current-limitations/)
   
   - [Intragenic deletions](https://arriba.readthedocs.io/en/latest/current-limitations/#intragenic-deletions)
   - [Memory consumption](https://arriba.readthedocs.io/en/latest/current-limitations/#memory-consumption)

9. [Internal algorithm](https://arriba.readthedocs.io/en/latest/internal-algorithm/)

   - [Read-level filters](https://arriba.readthedocs.io/en/latest/internal-algorithm/#read-level-filters)
   - [Event-level filters](https://arriba.readthedocs.io/en/latest/internal-algorithm/#event-level-filters)

