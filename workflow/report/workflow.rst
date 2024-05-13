Overview

1. Determine regions of homology using BLAST to align the pf15 intergenic
   regions to celegans transcripts (see "homology" section").

2. The workflow diverges into two sets of outputs: "neuron" and
   "neuron-germline". The homologous regions are annotated with the wormbase
   genes ids that are associated with the aligned celegans transcript, and then
   filtered into two lists. One contains only homologous regions annotated with
   genes that are in the "expressed in neuron" list, the other in the
   "expressed in neuron and germline" list (lists provided).

3. The fastq samples provided are then aligned to the pf15 genome and regions
   of the pf15 genome are identified as "putative small rnas" if they are
   covered by at least {{ snakemake.config["params"]["min_coverage"] }} read for all samples
   (see "putative small RNA" section).

4. The lists of homologous regions are then intersected with the putative small rnas.

5. Finally, the lists are further filtered to contain only regions associated
   with genes in the down regulated in C. elegans after PF15 lawn training for
   24 h (list provided).
