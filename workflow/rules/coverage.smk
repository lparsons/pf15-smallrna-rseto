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
