About
-----

Arriba is a command-line tool for the detection of gene fusions from RNA-Seq data. It was developed for the use in a clinical research setting. Therefore, short runtimes and high sensitivity were important design criteria. It is based on the ultrafast [STAR aligner](https://github.com/alexdobin/STAR), and the post-alignment runtime is typically just ~2 minutes. Arriba's workflow produces fully reusable alignments, which can serve as input to other common analyses, such as quantification of gene expression. In contrast to many other fusion detection tools which build on STAR, Arriba does not require to reduce the STAR parameter `--alignIntronMax` to detect fusions arising from focal deletions. Reducing this parameter impairs mapping of reads to genes with long introns and may affect expression quantification, hence.

Apart from gene fusions, Arriba can detect other structural rearrangements with potential clinical relevance, including viral integration sites, internal tandem duplications, whole exon duplications, intragenic inversions, enhancer hijacking events involving immunoglobulin/T-cell receptor loci, translocations affecting genes with many paralogs such as DUX4, and truncations of genes (i.e., breakpoints in introns or intergenic regions).

Arriba is the winner of the [DREAM SMC-RNA Challenge](https://www.synapse.org/SMC_RNA), an international competition organized by ICGC, TCGA, IBM, and Sage Bionetworks to determine the current gold standard for the detection of gene fusions from RNA-Seq data. The final results of the challenge are posted on the [Round 5 Leaderboard](https://www.synapse.org/#!Synapse:syn2813589/wiki/588511) and discussed in the accompanying [publication](https://doi.org/10.1016/j.cels.2021.05.021).

Get help
--------

Use the [GitHub issue tracker](https://github.com/suhrig/arriba/issues) to get help or to report bugs.

Citation
--------

Sebastian Uhrig, Julia Ellermann, Tatjana Walther, Pauline Burkhardt, Martina Fröhlich, Barbara Hutter, Umut H. Toprak, Olaf Neumann, Albrecht Stenzinger, Claudia Scholl, Stefan Fröhling and Benedikt Brors: *Accurate and efficient detection of gene fusions from RNA sequencing data.* Genome Research. March 2021 31: 448-460; Published in Advance January 13, 2021. doi: [10.1101/gr.257246.119](https://doi.org/10.1101/gr.257246.119)

License
-------

The code, software and database files of Arriba are distributed under the MIT/Expat License, with the exception of the script `draw_fusions.R`, which is distributed under the GNU GPL v3 due to dependencies on GPL-licensed R packages. The terms and conditions of both licenses can be found in the [LICENSE file](https://raw.githubusercontent.com/suhrig/arriba/master/LICENSE).

User manual
-----------

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

8. [Utility scripts](https://arriba.readthedocs.io/en/latest/utility-scripts/)

   - [Extract fusion-supporting alignments](https://arriba.readthedocs.io/en/latest/utility-scripts/#extract-fusion-supporting-alignments)
   - [Convert fusions.tsv to VCF](https://arriba.readthedocs.io/en/latest/utility-scripts/#convert-fusionstsv-to-vcf)
   - [Run Arriba on prealigned BAM file](https://arriba.readthedocs.io/en/latest/utility-scripts/#run-arriba-on-prealigned-bam-file)
   - [Quantify virus expression](https://arriba.readthedocs.io/en/latest/utility-scripts/#quantify-virus-expression)

9. [Current limitations](https://arriba.readthedocs.io/en/latest/current-limitations/)
   
   - [Intragenic deletions](https://arriba.readthedocs.io/en/latest/current-limitations/#intragenic-deletions)
   - [RefSeq annotation](https://arriba.readthedocs.io/en/latest/current-limitations/#refseq-annotation)
   - [Memory consumption](https://arriba.readthedocs.io/en/latest/current-limitations/#memory-consumption)
   - [Adapter trimming](https://arriba.readthedocs.io/en/latest/current-limitations/#adapter-trimming)
   - [Viral detection](https://arriba.readthedocs.io/en/latest/current-limitations/#viral-detection)
   - [Targeted sequencing](https://arriba.readthedocs.io/en/latest/current-limitations/#targeted-sequencing)
   - [Supporting read count vs. coverage](https://arriba.readthedocs.io/en/latest/current-limitations/#supporting-read-count-vs-coverage)

10. [Internal algorithm](https://arriba.readthedocs.io/en/latest/internal-algorithm/)

   - [Read-level filters](https://arriba.readthedocs.io/en/latest/internal-algorithm/#read-level-filters)
   - [Event-level filters](https://arriba.readthedocs.io/en/latest/internal-algorithm/#event-level-filters)

