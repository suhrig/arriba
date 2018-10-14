About
=====

Arriba is a command-line tool for the detection of gene fusions from RNA-Seq data. It was developed for the use in a clinical research setting.
Therefore, short runtimes and high sensitivity were important design criteria. It is based on the ultrafast [STAR aligner](https://github.com/alexdobin/STAR) and the post-alignment runtime is typically just ~2 minutes. In contrast to many other fusion detection tools which build on STAR, Arriba does not require to reduce the `alignIntronMax` parameter of STAR to detect small deletions.

Apart from gene fusions, Arriba can detect other structural rearrangements with potential clinical relevance, such as exon duplications or truncations of genes (i.e., breakpoints in introns and intergenic regions).

Arriba has been submitted to the [DREAM SMC RNA Challenge](https://www.synapse.org/SMC_RNA), an international competition organized by ICGC, TCGA, IBM, and Sage Bionetworks to determine the current gold standard for the detection of gene fusions from RNA-Seq data. As of [round 4](https://www.synapse.org/#!Synapse:syn2813589/wiki/423306), Arriba is the best-performing algorithm.

Installation
============

Please refer to the [quickstart guide](http://arriba.readthedocs.io/en/latest/quickstart/). **Note: You should not use `git clone` to download Arriba, because the git repository does not include the blacklist!**

