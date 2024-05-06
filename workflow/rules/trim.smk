rule rename_fastq:
    input:
        lambda wildcards: units.loc[
            (wildcards.sample_name, wildcards.unit_name), ["fq1"]
        ],
    output:
        fastq="results/fastq/{sample_name}-{unit_name}.fastq.gz",
    log:
        "logs/fastq/{sample_name}-{unit_name}.log",
    conda:
        "../envs/coreutils.yaml"
    shell:
        "ln -sr {input:q} {output:q} &> {log}"


rule cutadapt:
    input:
        fq1="results/fastq/{sample_name}-{unit_name}.fastq.gz",
    output:
        fastq="results/trimmed/{sample_name}-{unit_name}.fastq.gz",
        qc="results/trimmed/{sample_name}-{unit_name}.qc.txt",
    params:
        adapters=config["params"]["cutadapt"]["adapters"],
        extra=config["params"]["cutadapt"]["extra"],
    log:
        "logs/cutadapt/{sample_name}-{unit_name}.log",
    threads: 10
    wrapper:
        "v3.10.1/bio/cutadapt/se"
