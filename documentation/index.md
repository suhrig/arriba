About
-----

Arriba is a command-line tool for the detection of gene fusions from RNA-Seq data. It was developed for the use in a clinical setting.
Therefore, short runtimes and high sensitivity were important design criteria. It is based on the ultrafast [STAR aligner](https://github.com/alexdobin/STAR) and the post-alignment runtime is typically just ~2 minutes. In contrast to many other fusion detection tools which build on STAR, Arriba does not require to reduce the parameter `--alignIntronMax` of STAR to detect small deletions.

Apart from gene fusions, Arriba can detect other structural rearrangements with potential clinical relevance, such as exon duplications or truncations of genes (i.e., breakpoints in introns and intergenic regions).

Arriba has been submitted to the [DREAM SMC RNA Challenge](https://www.synapse.org/#!Synapse:syn2813589), an international competition organized by ICGC, TCGA, IBM, and Sage Bionetworks to determine the current gold standard for the detection of gene fusions from RNA-Seq data. As of [round 4](https://www.synapse.org/#!Synapse:syn2813589/wiki/423306), Arriba is the best-performing algorithm.

License
-------

Apart from the script `draw_fusions.R` all software/code of Arriba is disributed under the MIT/Expat License. The script `draw_fusions.R` is distributed under the GNU GPL v3 due to dependencies on GPL-licensed R packages. The terms and conditions of both licenses can be found in the [LICENSE file](https://raw.githubusercontent.com/suhrig/arriba/master/LICENSE).

Citing
------

A dedicated publication about Arriba has not been released yet. Until then, please refer to Arriba in your methods section as follows (or similar):
> We used Arriba (https://github.com/suhrig/arriba/) to detect gene fusions from RNA-Seq data.
