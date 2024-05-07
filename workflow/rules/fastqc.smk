rule fastqc:
    input:
        fastq="results/fastq/{sample}-{unit}.fastq.gz",
    output:
        html="results/qc/fastqc/{sample}-{unit}.html",
        zip="results/qc/fastqc/{sample}-{unit}_fastqc.zip",  # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
    params:
        "--quiet",
    log:
        "logs/fastqc/{sample}-{unit}.log",
    threads: 1
    wrapper:
        "v3.10.1/bio/fastqc"
