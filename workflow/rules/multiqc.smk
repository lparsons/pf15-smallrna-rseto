from os import path


rule multiqc:
    input:
        cutadapt=expand(
            "results/trimmed/{unit.sample_name}-{unit.unit_name}.qc.txt",
            unit=units.itertuples(),
        ),
        fastqc=expand(
            "results/qc/fastqc/{unit.sample_name}-{unit.unit_name}_fastqc.zip",
            unit=units.itertuples(),
        ),
        samtools_stats=expand(
            "results/aligned-pf15/{unit.sample_name}-{unit.unit_name}-samtools-stats.txt",
            unit=units.itertuples(),
        ),
    output:
        html=report(
            "results/qc/multiqc.html", caption="../report/multiqc.rst", category="QC"
        ),
        data=directory("results/qc/multiqc_data"),
    params:
        extra="--verbose",  # Optional: extra parameters for multiqc.
    log:
        "logs/multiqc.log",
    wrapper:
        "v3.10.2/bio/multiqc"


# rule multiqc:
#     input:
#         cutadapt=expand(
#             "results/trimmed/{unit.sample_name}-{unit.unit_name}.qc.txt",
#             unit=units.itertuples(),
#         ),
#         fastqc=expand(
#             "results/qc/fastqc/{unit.sample_name}-{unit.unit_name}_fastqc.zip",
#             unit=units.itertuples(),
#         ),
#         samtools_stats=expand(
#             "results/aligned-pf15/{unit.sample_name}-{unit.unit_name}-samtools-stats.txt",
#             unit=units.itertuples(),
#         ),
# #         counts_summary=expand(
# #             "results/counts/{unit.sample_name}-{unit.unit_name}-{strand}.txt.summary",
# #             unit=units.itertuples(),
# #             strand=("forward", "reverse"),
# #         ),
#     output:
#         html=report(
#             "results/qc/multiqc.html", caption="../report/multiqc.rst", category="QC"
#         ),
#     params:
#         output_dir=lambda wc, output: path.dirname(output.html),
#         output_name=lambda wc, output: path.basename(output.html),
#         input_dirs=lambda wc, input: set(path.dirname(fp) for fp in input),
#     log:
#         "logs/multiqc.log",
#     conda:
#         "../envs/multiqc.yaml"
#     shell:
#         "multiqc"
#         " --force"
#         " -o {params.output_dir:q}"
#         " -n {params.output_name:q}"
#         " {params.input_dirs}"
#         " &> {log:q}"
