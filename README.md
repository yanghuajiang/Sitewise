# Sitewise_comparision
The methods related to type identification for Mycoplasma pneumoniae.
# Quick Start (with examples)
~~~~~~~~~~~~~~
python Sitewise_comparision.py -r examples/GCA_900660465.1_50648_A01-3_genomic.fna -s examples/MP_P1.sites -q examples/HD0021418-B2.result.fastq -b examples/etoki.mapping.merged.bam
~~~~~~~~~~~~~~
# USAGE
~~~~~~~~~~~~~~
Usage: Sitewise_comparison.py [OPTIONS]

Options:
  -r, --ref TEXT   reference genome in fasta format  [required]
  -s, --site TEXT  type-specific SNPs in format of <seq_name> <site> <SNP>
                   [required]
  -q, --qry TEXT   query MAG in fastq format  [required]
  -b, --bam TEXT   bam file specifying mapping results  [required]
  --minimap2 TEXT  default: /users/softwares/bin/minimap2
  --samtools TEXT  default: /users/softwares/bin/samtools
  --help           Show this message and exit.
~~~~~~~~~~~~~~
