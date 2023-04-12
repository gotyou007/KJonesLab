# Detect SS18-SSX in patient-derived PDX 
About
-----

Arriba is a command-line tool for the detection of gene fusions from RNA-Seq data. It was developed for the use in a clinical research setting. Therefore, short runtimes and high sensitivity were important design criteria. It is based on the ultrafast [STAR aligner](https://github.com/alexdobin/STAR), and the post-alignment runtime is typically just ~2 minutes. Arriba's workflow produces fully reusable alignments, which can serve as input to other common analyses, such as quantification of gene expression. In contrast to many other fusion detection tools which build on STAR, Arriba does not require to reduce the STAR parameter `--alignIntronMax` to detect fusions arising from focal deletions. Reducing this parameter impairs mapping of reads to genes with long introns and may affect expression quantification, hence.

Apart from gene fusions, Arriba can detect other structural rearrangements with potential clinical relevance, including viral integration sites, internal tandem duplications, whole exon duplications, intragenic inversions, enhancer hijacking events involving immunoglobulin/T-cell receptor loci, translocations affecting genes with many paralogs such as DUX4, and truncations of genes (i.e., breakpoints in introns or intergenic regions).

Arriba is the winner of the [DREAM SMC-RNA Challenge](https://www.synapse.org/SMC_RNA), an international competition organized by ICGC, TCGA, IBM, and Sage Bionetworks to determine the current gold standard for the detection of gene fusions from RNA-Seq data. The final results of the challenge are posted on the [Round 5 Leaderboard](https://www.synapse.org/#!Synapse:syn2813589/wiki/588511) and discussed in the accompanying [publication](https://doi.org/10.1016/j.cels.2021.05.021).

1. Download reference: GRCm38+RefSeq
2. run_arriba.sh: a. detect chimeric reads; b. Arriba extracts the supplementary alignments from the given input file(s)
3. draw_fusions.R: Visualization
