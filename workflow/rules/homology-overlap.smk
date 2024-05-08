rule homology_overlap:
    input:
        left="results/homology/blast-pf15-intergenic-regions-celegans-transcripts-{min_length}-bp-{min_pct_id}-pctid-genes-expressed-{expression}.bed",
        right="results/aligned-pf15-coverage/putative-smallrna-regions.bed",
    output:
        "results/homology-overlap/homology-{min_length}-bp-{min_pct_id}-pctid-genes-expressed-{expression}-overlap-putative-smallrna.bed"
    log:
        "logs/homology-overlap/homology-{min_length}-bp-{min_pct_id}-pctid-genes-expressed-{expression}-overlap-putative-smallrna.log"
    wrapper:
        "v3.10.2/bio/bedtools/intersect"
