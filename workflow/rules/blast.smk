rule blast_pf15_intergenic_celegans_transcripts:
    input:
        celegrans_transcripts="results/celegans-transcripts.fasta",
        pf15_intergenic_regions=f"{config['reference']['fasta']}.intergenic-regions.fasta",
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
       echo -e "qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore" > {output:q} 2> {log:q}
       awk '$3 >= {wildcards.minpctid} && $4 >= {wildcards.minlen} {{print $0}}' {input:q} >> {output:q} 2>> {log:q}
       """


# rule get_fasta:
#    input:
#        genome=config["reference"]["fasta"],
#        regions="results/diffexp_windows/{contrast}/significant_regions.gtf",
#    output:
#        "results/diffexp_windows/{contrast}/significant_regions.fasta",
#    log:
#        "logs/diffexp_windows/{contrast}/significant_regions.fasta.log",
#    conda:
#        "../envs/bedtools.yaml"
#    shell:
#        "bedtools getfasta -s -fi {input.genome:q} -bed {input.regions:q} -fo {output:q} 2> {log:q}"
#
#
# rule unzip_celegans_fasta:
#    input:
#        config["reference"]["celegans"],
#    output:
#        temp("results/celegans.fasta"),
#    log:
#        "logs/unzip_celegans.log",
#    conda:
#        "../envs/coreutils.yaml"
#    shell:
#        "gunzip --stdout --decompress --force {input:q} > {output:q} 2> {log:q}"
#
#
# rule blast_regions_celegans:
#    input:
#        genome="results/celegans.fasta",
#        regions="results/diffexp_windows/{contrast}/significant_regions.fasta",
#    output:
#        "results/diffexp_significant_regions/{contrast}/significant_regions.blast.txt",
#    log:
#        "logs/diffexp_significant_regions/{contrast}/significant_regions.blast.log",
#    conda:
#        "../envs/blast.yaml"
#    shell:
#        """
#        blastn -subject {input.genome:q} -query {input.regions:q} -task blastn-short -outfmt 6 -ungapped > {output:q} 2> {log:q}
#        """
#
#
# rule filter_blast:
#    input:
#        "results/diffexp_significant_regions/{contrast}/significant_regions.blast.txt",
#    output:
#        report(
#            "results/diffexp_significant_regions/{contrast}/significant_regions.blast.length-{length,[0-9]+}.txt",
#            caption="../report/filter_blast.rst",
#            category="Differential Expression (Significant Regions)",
#            subcategory="{contrast}",
#        ),
#    log:
#        "logs/diffexp_significant_regions/{contrast}/significant_regions.blast.length-{length}.log",
#    conda:
#        "../envs/coreutils.yaml"
#    shell:
#        """
#        echo -e "qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore" > {output:q} 2> {log:q}
#        awk '$4 >= {wildcards.length} {{print $0}}' {input:q} >> {output:q} 2>> {log:q}
#        """
#
#
# rule group_blast:
#    input:
#        blast_results="results/diffexp_significant_regions/{contrast}/significant_regions.blast.length-{length}.txt",
#        transcript_gene_map="resources/transcript_gene_map.tsv",
#    output:
#        "results/diffexp_significant_regions/{contrast}/significant_regions.blast.length-{length,[0-9]+}.grouped.txt",
#    log:
#        "logs/diffexp_significant_regions/{contrast}/significant_regions.blast.length-{length}.grouped.log",
#    conda:
#        "../envs/bedtools.yaml"
#    shell:
#        "join -1 2 -2 1 "
#        "<(tail -n +2 {input.blast_results:q} | sort -k2b,2) "
#        "<(sort -k1b,1 {input.transcript_gene_map:q}) "
#        "-a 1 -t $'\t' | "
#        "sort -k2b,2 | "
#        "bedtools groupby -opCols 13 -ops distinct -grp 2 | "
#        "sort -k1b,1 > {output:q} 2> {log:q}"
#
#
# rule annotate_diffexp_with_blast:
#    input:
#        diffexp_results="results/diffexp_significant_regions/{contrast}/results_significant_regions.tsv",
#        grouped_blast_results=expand(
#            "results/diffexp_significant_regions/{{contrast}}/significant_regions.blast.length-{length}.grouped.txt",
#            length=config["params"]["blast"]["min_length"],
#        ),
#    output:
#        report(
#            "results/diffexp_significant_regions/{contrast}/results_significant_regions_celegans_genes.tsv",
#            caption="../report/results_significant_regions_blast.rst",
#            category="Differential Expression (Significant Regions)",
#            subcategory="{contrast}",
#        ),
#    log:
#        "logs/diffexp_significant_regions/{contrast}/results_significant_regions_celegans_genes.log",
#    conda:
#        "../envs/coreutils.yaml"
#    shell:
#        'echo -e "region\tbaseMean\tlog2FoldChange\tlfcSE\tpvalue\tpadj\tc_elegans_gene_ids" > {output:q} 2> {log:q} && '
#        "join -a 1 -t $'\t' -j 1 "
#        "<(tail -n +2 {input.diffexp_results:q} | sort -k1b,1) "
#        "<(sort -k1b,1 {input.grouped_blast_results:q} ) >> "
#        "{output:q} 2>> {log:q}"