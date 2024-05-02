rule celegans_transcripts:
    input:
        storage(config["reference"]["celegans"]["transcripts_fasta"]),
    output:
        temp("results/celegans-transcripts.fasta"),
    log:
        "logs/unzip-celegans-transcripts-fasta.log",
    conda:
        "../envs/coreutils.yaml"
    shell:
        "gunzip --stdout --decompress --force {input:q} > {output:q} 2> {log:q}"


rule celegans_gtf:
    input:
        storage(config["reference"]["celegans"]["genome_gtf"]),
    output:
        temp("results/celegans-genes.gtf"),
    log:
        "logs/unzip-celegans-genome-gtf.log",
    conda:
        "../envs/coreutils.yaml"
    shell:
        "gunzip --stdout --decompress --force {input:q} > {output:q} 2> {log:q}"


rule gtf_to_gene_transcripts:
    input:
        "results/celegans-genes.gtf",
    output:
        "results/celegans-genes-transcipts.tsv",
    log:
        "logs/celegans-gtf-to-genes-transcipts.log",
    conda:
        "../envs/coreutils.yaml"
    shell:
        'gffread {input:q} --table "@geneid, transcript_id" -o {output:q} 2> {log:q}'


# rule star_index:
#     input:
#         fasta=config["reference"]["fasta"],
#     output:
#         directory("resources/star_genome"),
#     threads: 4
#     params:
#         extra=config["params"]["star_index"],
#     log:
#         "logs/star_index_genome.log",
#     cache: True
#     wrapper:
#         "v3.0.2/bio/star/index"
#
#
# rule gff_to_gtf:
#     input:
#         gff=config["reference"]["gff"],
#     output:
#         gtf=config["reference"]["gtf"],
#     conda:
#         "../envs/gffread.yaml"
#     log:
#         "logs/gffread_gff_to_gtf.log",
#     cache: True
#     shell:
#         "gffread -E -T --force-exons -o {output:q} {input:q} &> {log:q}"
#
#
rule samtools_faidx:
    input:
        fasta=config["reference"]["pf15"]["fasta"],
    output:
        fai=f"{config['reference']['pf15']['fasta']}.fai",
    params:
        "",  # optional params string
    log:
        "logs/samtools-faidx-genome.log",
    wrapper:
        "v3.0.2/bio/samtools/faidx"


# rule genmap_index:
#     input:
#         fasta=config["reference"]["fasta"],
#     output:
#         idx=f"{config['reference']['fasta']}.genmap.idx",
#     params:
#         config["params"]["genmap_index"],
#     log:
#         "logs/genmap_index.log",
#     conda:
#         "../envs/genmap.yaml"
#     shell:
#         """
#         genmap index --fasta-file {input.fasta:q} --index {output.idx:q} &> {log:q}
#         """
#
#
# rule mappability:
#     input:
#         idx=f"{config['reference']['fasta']}.genmap.idx",
#     output:
#         txt=expand(
#             "{fasta}.genmap.mappability.length-{{length}}.txt",
#             fasta=config["reference"]["fasta"],
#         ),
#         wig=expand(
#             "{fasta}.genmap.mappability.length-{{length}}.wig",
#             fasta=config["reference"]["fasta"],
#         ),
#         chromsizes=expand(
#             "{fasta}.genmap.mappability.length-{{length}}.chrom.sizes",
#             fasta=config["reference"]["fasta"],
#         ),
#         bedgraph=report(
#             expand(
#                 "{fasta}.genmap.mappability.length-{{length}}.bedgraph",
#                 fasta=config["reference"]["fasta"],
#             ),
#             caption="../report/mappability_bedgraph.rst",
#             category="Genome stats",
#         ),
#     params:
#         outputbase=expand(
#             "{fasta}.genmap.mappability.length-{{length}}",
#             fasta=config["reference"]["fasta"],
#         ),
#         other=config["params"]["genmap"],
#     threads: 10
#     log:
#         "logs/mappability-length-{length}.log",
#     conda:
#         "../envs/genmap.yaml"
#     shell:
#         "genmap map "
#         "--length {wildcards.length} "
#         " --index {input.idx:q} "
#         " --output {params.outputbase:q} "
#         " {params.other} --txt --wig --bedgraph "
#         " --threads {threads} "
#         " &> {log:q}"
#
#
# rule bedtools_makewindows:
#     input:
#         fai=f"{config['reference']['fasta']}.fai",
#     output:
#         windows=f"{config['reference']['fasta']}.windows.bed",
#     params:
#         window_size=config["params"]["genome_window_counts"]["window_size"],
#     conda:
#         "../envs/bedtools.yaml"
#     log:
#         "logs/bedtools_makewindows.log",
#     shell:
#         "bedtools makewindows "
#         "-g {input.fai:q} "
#         "-w {params.window_size} "
#         '| awk -F "\\t" \'OFS="\\t" {{print $0, $1 ":" $2 "-" $3 "(+)", "0", "+"; print $0, $1 ":" $2 "-" $3 "(-)", "0", "-"}}\' '
#         "> {output.windows:q} "
#         "2> {log:q}"
#
#
rule intergenic_regions:
    input:
        fai=f"{config['reference']['pf15']['fasta']}.fai",
        gtf=config["reference"]["pf15"]["gtf"],
    output:
        intergenic_regions=f"{config['reference']['pf15']['fasta']}.intergenic-regions.bed",
    conda:
        "../envs/bedtools.yaml"
    log:
        "logs/intergenic_regions.log",
    shell:
        "bedtools complement "
        "-i {input.gtf:q} "
        " -g {input.fai:q} "
        "> {output.intergenic_regions:q} "
        "2> {log:q}"


rule getfasta:
    input:
        fasta=f"{config['reference']['pf15']['fasta']}",
        fai=f"{config['reference']['pf15']['fasta']}.fai",
        bed=f"{config['reference']['pf15']['fasta']}.{{regions}}.bed",
    output:
        fasta=f"{config['reference']['pf15']['fasta']}.{{regions}}.fasta",
    conda:
        "../envs/bedtools.yaml"
    log:
        "logs/getfasta-{regions}.log",
    shell:
        "bedtools getfasta "
        "-fi {input.fasta:q} "
        " -bed {input.bed:q} "
        "> {output.fasta:q} "
        "2> {log:q}"


# rule bedtools_intergenic_windows:
#     input:
#         windows=f"{config['reference']['fasta']}.windows.bed",
#         gtf=config["reference"]["gtf"],
#     output:
#         intergenic_windows=f"{config['reference']['fasta']}.intergenic_windows.bed",
#     params:
#         window_size=config["params"]["genome_window_counts"]["window_size"],
#     conda:
#         "../envs/bedtools.yaml"
#     log:
#         "logs/bedtools_intergenic_windows.log",
#     shell:
#         "bedtools subtract "
#         "-A "
#         "-a {input.windows:q} "
#         "-b {input.gtf:q} "
#         "> {output.intergenic_windows:q} "
#         "2> {log:q}"
#
#
# rule bedtogenepred:
#     input:
#         bed=f"{config['reference']['fasta']}.intergenic_windows.bed",
#     output:
#         genepred=temp(f"{config['reference']['fasta']}.intergenic_windows.genepred"),
#     conda:
#         "../envs/ucsc-bedtogenepred.yaml"
#     log:
#         "logs/ucsc-bedtogenepred.log",
#     shell:
#         "bedToGenePred {input.bed:q} {output.genepred:q} &> {log}"
#
#
# rule genepredtogtf:
#     input:
#         genepred=f"{config['reference']['fasta']}.intergenic_windows.genepred",
#     output:
#         gtf=f"{config['reference']['fasta']}.intergenic_windows.gtf",
#     conda:
#         "../envs/ucsc-genepredtogtf.yaml"
#     log:
#         "logs/ucsc-genepredtogtf.log",
#     shell:
#         "genePredToGtf file {input.genepred:q} {output.gtf:q} &> {log}"
#
#
# rule transcript_to_gene_map:
#     input:
#         "results/celegans.fasta",
#     output:
#         "resources/transcript_gene_map.tsv",
#     conda:
#         "../envs/coreutils.yaml"
#     log:
#         "logs/transcript_to_gene_map.log",
#     shell:
#         """
#         awk -v OFS='\\t' '$1 ~ /\\>/ {{gsub(">","",$1); gsub("gene=","",$2); print $1, $2}}' {input:q} >{output:q} 2>{log:q}
#         """
