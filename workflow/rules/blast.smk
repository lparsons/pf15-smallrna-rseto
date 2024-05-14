rule blast_pf15_intergenic_celegans_transcripts:
    input:
        celegrans_transcripts="results/celegans-transcripts.fasta",
        pf15_intergenic_regions=f"{config['reference']['pf15']['fasta']}.intergenic-regions.fasta",
    output:
        "results/homology/blast-pf15-intergenic-regions-celegans-transcripts.tsv",
    log:
        "logs/homology/blast-pf15-intergenic-regions-celegans-transcripts.log",
    conda:
        "../envs/blast.yaml"
    shell:
        "blastn "
        "-subject {input.celegrans_transcripts:q} "
        "-query {input.pf15_intergenic_regions:q} "
        "-task blastn-short "
        "-outfmt 6 "
        "-ungapped "
        "> {output:q} "
        "2> {log:q}"


rule filter_homology:
    input:
        "results/homology/blast-pf15-intergenic-regions-celegans-transcripts.tsv",
    output:
        "results/homology/blast-pf15-intergenic-regions-celegans-transcripts-{minlen}-bp-{minpctid}-pctid.tsv",
    log:
        "logs/homology/blast-pf15-intergenic-regions-celegans-transcripts-{minlen}-bp-{minpctid}-pctid.log",
    conda:
        "../envs/coreutils.yaml"
    shell:
        """
        echo -e "qseqid\\tsseqid\\tpident\\tlength\\tmismatch\\tgapopen\\tqstart\\tqend\\tsstart\\tsend\\tevalue\\tbitscore" > {output:q} 2> {log:q}
        awk 'BEGIN{{FS="\\t"; OFS=FS}} $3 >= {wildcards.minpctid} && $4 >= {wildcards.minlen} {{print $0}}' {input:q} >> {output:q} 2>> {log:q}
        """


rule sort_tsv_with_header:
    input:
        "{file}.{ext}",
    output:
        "{file}.sorted-col{colnum}.{ext}",
    log:
        "logs/{file}.sorted-col{colnum}.{ext}.log",
    conda:
        "../envs/coreutils.yaml"
    shell:
        """
        cat {input:q} | (read -r; printf "%s\n" "$REPLY"; LC_ALL=C sort -k {wildcards.colnum}b,{wildcards.colnum}) > {output:q} 2> {log:q}
        """


rule annotate_homology:
    input:
        homologous_regions="results/homology/blast-pf15-intergenic-regions-celegans-transcripts-{minlen}-bp-{minpctid}-pctid.sorted-col2.tsv",
        genes_transcripts="results/celegans-genes-transcipts.sorted-col2.tsv",
    output:
        report(
            "results/homology/blast-pf15-intergenic-regions-celegans-transcripts-{minlen}-bp-{minpctid}-pctid-genes.tsv",
            caption="../report/homology.rst",
            category="homology",
        ),
    log:
        "logs/homology/blast-pf15-intergenic-regions-celegans-transcripts-{minlen}-bp-{minpctid}-pctid-genes.log",
    conda:
        "../envs/coreutils.yaml"
    shell:
        """
        LC_ALL=C join -t $'\\t' --header -j2 {input.homologous_regions:q} {input.genes_transcripts:q} > {output:q} 2>> {log:q}
        """


rule homology_genelist:
    input:
        "results/homology/blast-pf15-intergenic-regions-celegans-transcripts-{minlen}-bp-{minpctid}-pctid-genes.tsv",
    output:
        report(
            "results/homology/blast-pf15-intergenic-regions-celegans-transcripts-{minlen}-bp-{minpctid}-pctid-genes.unique-genes.txt",
            caption="../report/homology-genelist.rst",
            category="homology",
        ),
    log:
        "logs/homology/blast-pf15-intergenic-regions-celegans-transcripts-{minlen}-bp-{minpctid}-pctid-genes.unique-genes.log",
    conda:
        "../envs/coreutils.yaml"
    shell:
        """
        awk 'NR>1 {{print $13}}' {input:q} | LC_ALL=C sort | uniq -c > {output:q} 2> {log:q}
        """


rule filter_homology_by_gene:
    input:
        homologous_regions="results/homology/blast-pf15-intergenic-regions-celegans-transcripts-{minlen}-bp-{minpctid}-pctid-genes.sorted-col13.tsv",
        gene_list="resources/{expression}-expressed-genes.sorted-col1.tsv",
    output:
        report(
            "results/homology/blast-pf15-intergenic-regions-celegans-transcripts-{minlen}-bp-{minpctid}-pctid-genes-expressed-{expression}.tsv",
            caption="../report/homology-expressed.rst",
            category="{expression}",
        ),
    log:
        "logs/homology/blast-pf15-intergenic-regions-celegans-transcripts-{minlen}-bp-{minpctid}-pctid-genes-expressed-{expression}.log",
    conda:
        "../envs/coreutils.yaml"
    shell:
        """
        LC_ALL=C join -t $'\\t' --header -1 1 -2 13 {input.gene_list:q} {input.homologous_regions:q} > {output:q} 2>> {log:q}
        """


rule homology_filter_by_gene_genelist:
    input:
        "results/homology/blast-pf15-intergenic-regions-celegans-transcripts-{minlen}-bp-{minpctid}-pctid-genes-expressed-{expression}.tsv",
    output:
        report(
            "results/homology/blast-pf15-intergenic-regions-celegans-transcripts-{minlen}-bp-{minpctid}-pctid-genes-expressed-{expression}.unique-genes.txt",
            caption="../report/homology-expressed-genelist.rst",
            category="{expression}",
        ),
    log:
        "logs/homology/blast-pf15-intergenic-regions-celegans-transcripts-{minlen}-bp-{minpctid}-pctid-genes-expressed-{expression}.unique-genes.log",
    conda:
        "../envs/coreutils.yaml"
    shell:
        """
        awk 'NR>1 {{print $1}}' {input:q} | LC_ALL=C sort | uniq -c > {output:q} 2> {log:q}
        """


rule homology_to_bed:
    input:
        "results/homology/blast-pf15-intergenic-regions-celegans-transcripts-{minlen}-bp-{minpctid}-pctid-genes-expressed-{expression}.tsv",
    output:
        "results/homology/blast-pf15-intergenic-regions-celegans-transcripts-{minlen}-bp-{minpctid}-pctid-genes-expressed-{expression}.bed",
    log:
        "results/homology/blast-pf15-intergenic-regions-celegans-transcripts-{minlen}-bp-{minpctid}-pctid-genes-expressed-{expression}-bed.log",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/blast-to-bed.py"
