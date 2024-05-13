rule homology_overlap:
    input:
        left="results/homology/blast-pf15-intergenic-regions-celegans-transcripts-{min_length}-bp-{min_pct_id}-pctid-genes-expressed-{expression}.bed",
        right="results/aligned-pf15-coverage/putative-smallrna-regions.bed",
    output:
        report(
            "results/homology-overlap/homology-{min_length}-bp-{min_pct_id}-pctid-genes-expressed-{expression}-overlap-putative-smallrna.bed",
            caption="../report/homology-overlap.rst",
            category="{expression}",
        ),
    log:
        "logs/homology-overlap/homology-{min_length}-bp-{min_pct_id}-pctid-genes-expressed-{expression}-overlap-putative-smallrna.log",
    wrapper:
        "v3.10.2/bio/bedtools/intersect"


rule homology_overlap_genelist:
    input:
        "results/homology-overlap/homology-{min_length}-bp-{min_pct_id}-pctid-genes-expressed-{expression}-overlap-putative-smallrna.bed",
    output:
        report(
            "results/homology-overlap/homology-{min_length}-bp-{min_pct_id}-pctid-genes-expressed-{expression}-overlap-putative-smallrna.unique-genes.txt",
            caption="../report/homology-overlap-genelist.rst",
            category="{expression}",
        ),
    log:
        "logs/homology-overlap/homology-{min_length}-bp-{min_pct_id}-pctid-genes-expressed-{expression}-overlap-putative-smallrna.unique-genes.log",
    conda:
        "../envs/coreutils.yaml"
    shell:
        """
        awk '{{print $4}}' {input:q} | LC_ALL=C sort | uniq -c > {output:q} 2> {log:q}
        """


rule homology_overlap_downreg:
    input:
        homology_regions="results/homology-overlap/homology-{min_length}-bp-{min_pct_id}-pctid-genes-expressed-{expression}-overlap-putative-smallrna.sorted-col4.bed",
        downreg_genes="resources/downreg-genes-pf15-p0.sorted-col1.tsv",
    output:
        report(
            "results/homology-overlap/homology-{min_length}-bp-{min_pct_id}-pctid-genes-expressed-{expression}-overlap-putative-smallrna-downreg.bed",
            caption="../report/homology-overlap-downreg.rst",
            category="{expression}",
        ),
    log:
        "logs/homology-overlap/homology-{min_length}-bp-{min_pct_id}-pctid-genes-expressed-{expression}-overlap-putative-smallrna-downreg.log",
    conda:
        "../envs/coreutils.yaml"
    shell:
        """
        LC_ALL=C join -t $'\\t' --header -1 4 -2 1 {input.homology_regions:q} {input.downreg_genes:q} > {output:q} 2>> {log:q}
        """


rule homology_overlap_downreg_genelist:
    input:
        "results/homology-overlap/homology-{min_length}-bp-{min_pct_id}-pctid-genes-expressed-{expression}-overlap-putative-smallrna-downreg.bed",
    output:
        report(
            "results/homology-overlap/homology-{min_length}-bp-{min_pct_id}-pctid-genes-expressed-{expression}-overlap-putative-smallrna-downreg.unique-genes.txt",
            caption="../report/homology-overlap-downreg-genelist.rst",
            category="{expression}",
        ),
    log:
        "logs/homology-overlap/homology-{min_length}-bp-{min_pct_id}-pctid-genes-expressed-{expression}-overlap-putative-smallrna-downreg.unique-genes.log",
    conda:
        "../envs/coreutils.yaml"
    shell:
        """
        awk '{{print $1}}' {input:q} | LC_ALL=C sort | uniq -c > {output:q} 2> {log:q}
        """
