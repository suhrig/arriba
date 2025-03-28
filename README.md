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

Please refer to the [wiki pages](http://github.com/suhrig/arriba/wiki/01-Home) for installation instructions and information about usage. **Note: You should not use `git clone` to download Arriba, because the git repository does not include the blacklist and other database files!**

