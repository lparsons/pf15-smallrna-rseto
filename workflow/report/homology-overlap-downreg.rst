Homologous regions in genes expressed in "{{ snakemake.wildcards.expression }}"
intersected with putative small RNAs (regions of the pf15 genome covered by at
least {{ snakemake.wildcards.threshold }} reads in at least {{ snakemake.params.min_samples }} samples).

Further filterd to only genes downregulated in c. elegans after pf15 lawn
training for 24h (determined by list provided).
