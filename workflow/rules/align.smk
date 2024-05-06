rule align:
    input:
        fq1="results/trimmed/{sample}-{unit}.fastq.gz",
        idx="resources/{genome}-star-genome",
    output:
        # see STAR manual for additional output files
        aln="results/aligned-{genome}/{sample}-{unit}/Aligned.sortedByCoord.out.bam",
        log="results/aligned-{genome}/{sample}-{unit}/Log.out",
        unmapped="results/aligned-{genome}/{sample}-{unit}/unmapped.fastq.gz",
    log:
        "logs/aligned-{genome}/{sample}-{unit}.log",
    params:
        extra="--outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 10000000000 {}".format(
            config["params"]["star"]
        ),
    threads: 10
    resources:
        mem_mb=80000,
        time=120,
    wrapper:
        "v3.10.2/bio/star/align"


rule samtools_index:
    input:
        "results/aligned-{genome}/{sample}-{unit}/Aligned.sortedByCoord.out.bam",
    output:
        "results/aligned-{genome}/{sample}-{unit}/Aligned.sortedByCoord.out.bam.bai",
    log:
        "logs/aligned-{genome}/{sample}-{unit}-bamindex.log",
    params:
        "",  # optional params string
    # Samtools takes additional threads through its option -@
    threads: 4  # This value - 1 will be sent to -@
    wrapper:
        "v3.10.2/bio/samtools/index"


rule samtools_stats:
    input:
        "results/aligned-{genome}/{sample}-{unit}/Aligned.sortedByCoord.out.bam",
    output:
        "results/aligned-{genome}/{sample}-{unit}-samtools-stats.txt",
    params:
        extra="",  # Optional: extra arguments.
        region="",  # Optional: region string.
    log:
        "logs/aligned-{genome}/{sample}-{unit}-samtools-stats.log",
    wrapper:
        "v3.10.2/bio/samtools/stats"
