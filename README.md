# Sitewise_comparision
The methods related to type identification for Mycoplasma pneumoniae.
# Quick Start (with examples)
~~~~~~~~~~~~~~
python Sitewise_comparision.py -r examples/GCA_900660465.1_50648_A01-3_genomic.fna -s examples/MP_P1.sites -q examples/HD0021418-B2.result.fastq -b examples/etoki.mapping.merged.bam
~~~~~~~~~~~~~~
# USAGE
~~~~~~~~~~~~~~
Usage: sitewise_comparison.py [OPTIONS]

Options:
  -r, --ref TEXT   reference sequence  [required]
  -s, --site TEXT  sites of concern, in format of <seq_name> <site>
                   [required]

  -q, --qry TEXT   query sequence  [required]
  -b, --bam TEXT   bam file specifying mapping results
  -o, --out TEXT   output file
  --minimap2 TEXT  default: /titan/softwares/bin/minimap2
  --samtools TEXT  default: /titan/softwares/bin/samtools
  --help           Show this message and exit.
~~~~~~~~~~~~~~
