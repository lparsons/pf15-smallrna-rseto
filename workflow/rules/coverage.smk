rule pf15_genome_coverage:
    input:
        "results/aligned-{genome}/{sample}-{unit}/Aligned.sortedByCoord.out.bam",
    output:
        "results/aligned-{genome}-coverage/{sample}-{unit}-genomecov.bedgraph",
    log:
        "logs/aligned-{genome}-coverage/{sample}-{unit}-genomecov.log",
    params:
        "-bg",
    wrapper:
        "v3.10.2/bio/bedtools/genomecov"


rule pf15_genome_coverage_min:
    input:
        "results/aligned-{genome}-coverage/{sample}-{unit}-genomecov.bedgraph",
    output:
        "results/aligned-{genome}-coverage/{sample}-{unit}-genomecov-min{threshold}.bedgraph",
    log:
        "logs/aligned-{genome}-coverage/{sample}-{unit}-genomecov-min{threshold}.log",
    conda:
        "../envs/coreutils.yaml"
    shell:
        "awk '$4>={wildcards.threshold} {{print $0}}' {input:q} > {output:q} 2> {log:q}"


rule pf15_all_coverage:
    input:
        expand(
            "results/aligned-{genome}-coverage/{unit.sample_name}-{unit.unit_name}-genomecov-min{threshold}.bedgraph",
            unit=units.itertuples(),
            genome="pf15",
            threshold=config["params"]["min_coverage"],
        ),
    output:
        "results/aligned-{genome}-coverage/putative-smallrna-regions.bed",
    log:
        "logs/aligned-{genome}-coverage/putative-smallrna-regions.log",
    params:
        min_samples=len(samples),
    conda:
        "../envs/bedtools.yaml"
    shell:
        "bedtools multiinter -header -i {input:q} | awk '$4=={params.min_samples} {{printf (\"%s\\t%i\\t%i\\n\", $1, $2, $3)}}' > {output:q} 2> {log:q}"
