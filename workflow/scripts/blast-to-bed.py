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
                (start, end) = positions.split("-")
                writer.writerow(
                    {
                        "chrom": chrom,
                        "start": int(start),
                        "end": int(end),
                        "WormBase Gene ID": row["WormBase Gene ID"],
                        "WormBase Transcript ID": row["sseqid"],
                    }
                )


blast_to_bed(snakemake.input[0], snakemake.output[0])
