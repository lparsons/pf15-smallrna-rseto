rule homology_overlap:
    input:
        left="results/homology/blast-pf15-intergenic-regions-celegans-transcripts-{min_length}-bp-{min_pct_id}-pctid-genes-expressed-{expression}.bed",
        right="results/aligned-pf15-coverage/putative-smallrna-regions.bed",
    output:
        "results/homology-overlap/homology-{min_length}-bp-{min_pct_id}-pctid-genes-expressed-{expression}-overlap-putative-smallrna.bed",
    log:
        "logs/homology-overlap/homology-{min_length}-bp-{min_pct_id}-pctid-genes-expressed-{expression}-overlap-putative-smallrna.log",
    wrapper:
        "v3.10.2/bio/bedtools/intersect"


rule homology_overlap_downreg:
    input:
        homology_regions="results/homology-overlap/homology-{min_length}-bp-{min_pct_id}-pctid-genes-expressed-{expression}-overlap-putative-smallrna.sorted-col4.bed",
        downreg_genes="resources/downreg-genes-pf15-p0.sorted-col1.tsv",
    output:
        "results/homology-overlap/homology-{min_length}-bp-{min_pct_id}-pctid-genes-expressed-{expression}-overlap-putative-smallrna-downreg.bed",
    log:
        "logs/homology-overlap/homology-{min_length}-bp-{min_pct_id}-pctid-genes-expressed-{expression}-overlap-putative-smallrna-downreg.log",
    conda:
        "../envs/coreutils.yaml"
    shell:
        """
        LC_ALL=C join -t $'\\t' --header -1 4 -2 1 {input.homology_regions:q} {input.downreg_genes:q} > {output:q} 2>> {log:q}
        """
