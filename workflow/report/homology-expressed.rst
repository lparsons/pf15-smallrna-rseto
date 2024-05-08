Homologous regions in genes expressed in {{ snakemake.wildcards.expression }}
(determined by list provided).

The columns are defined as follows:

* WormBase Gene ID: wormbase gene id for gene containing the transcript aligned
* sseqid: celegans transcript id
* qseqid: "name" of pf15 intergenic region (location in the genome)
* pident: percent identity of the match
* length: length of the match
* mismatch: number of mismatches
* gapopen: number of gap opens
* qstart: start of alignment in pf15 intergenic region
* qend: end of alignment in pf15 intergenic region
* sstart: start of alignment in celegans transcript
* send: end of alignment in celegans transcript
* evalue: expect value
* bitscore: bitscore
