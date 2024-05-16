import csv


def blast_to_bed(blast_input, bed_output):

    with open(blast_input) as blast_file:
        tsv_file = csv.DictReader(blast_file, delimiter="\t")

        with open(bed_output, "w", newline="") as bed_outfile:
            fieldnames = [
                "chrom",
                "start",
                "end",
                "WormBase Gene ID",
                "WormBase Transcript ID",
            ]
            writer = csv.DictWriter(bed_outfile, fieldnames=fieldnames, delimiter="\t")

            writer.writeheader()
            for row in tsv_file:
                region = row["qseqid"]
                (chrom, positions) = region.split(":")
                (qseq_start, qseq_end) = positions.split("-")
                qstart = int(row["qstart"])
                qend = int(row["qend"])
                offset = int(qseq_start) - 1
                start = qstart + offset
                end = qend + offset
                writer.writerow(
                    {
                        "chrom": chrom,
                        "start": start - 1,
                        "end": end,
                        "WormBase Gene ID": row["WormBase Gene ID"],
                        "WormBase Transcript ID": row["sseqid"],
                    }
                )


blast_to_bed(snakemake.input[0], snakemake.output[0])
