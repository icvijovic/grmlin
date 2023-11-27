
## About

This repository contains a python package that calls immunoglobulin V germline alleles from deep BCR sequencing data, using an approach described [here](https://biorxiv.org/cgi/content/short/2023.11.25.568681v1). If a database of known alleles is provided, it also annotates inferred germline sequences with respect to the sequences in that database.

## Planned updates
- Current algorithm doesn't rely on a database of known alleles, except to parse out V gene sequences from recombined VDJ sequences. Incorporating a VDJ sequence parsing %protocol would allow the user to skip this trivial reliance on an external database, and make it possible to provide a simple fasta file as input.
- Currently, various parameters based on which germline alleles are distinguished from convergent hypermutation events are passed from command line (default values calibrated for deep BCR RNA sequencing of healthy human peripheral blood). A future update should make it possible to self-consistently set these thresholds from the data.

## Requirements

required: see environment.yaml

optional: version 4 BLAST database of known alleles (for the annotation of inferred germlines)

## Installation

Clone repository, install dependencies and add grmlin to your path.

## Citing grmlin

When using grmlin in a publication, please cite the following paper:

"Reference-Free Germline Immunoglobulin Allele Discovery from B Cell Receptor Sequencing Data"

Ivana Cvijovic, Elizabeth R Jerison, Stephen R Quake
_bioRxiv_ 2023.11.25.568681; doi: https://doi.org/10.1101/2023.11.25.568681

## Usage
```
usage: grmlin [-h] [-annotate path_to_blast_db]
                   [-outdir output_directory]
                   [-trim_primers path_to_primer_fasta]
                   [--skip_preprocess] [--verbose]
                   [-max_log_v_evalue float_value]
                   [-max_log_j_evalue float_value]
                   [-min_v_sequence_length float_value]
                   [-min_j_sequence_length float_value]
                   [--allow_unproductive] [--allow_missing_cdr3]
                   [--allow_ns_in_sequence]
                   [-max_primer_error int_value] [--reverse_primers]
                   [-lineage_clustering_cutoff float_value]
                   [-hierarchical_clustering_method method]
                   [-min_germline_usage_fraction float_value]
                   [-min_germline_num_lineages float_value]
                   [-cloud_radius float_value]
                   [-min_usage_fraction_within_cloud float_value]
                            input_file_path

Constructs a database of Ig V germlines.

positional arguments:
  input_file_path       path to input AIRR-formatted table containing Ig
                        sequences

optional arguments:
  -h, --help            show this help message and exit
  -annotate path_to_blast_db
                        path to database containing known Ig V alleles; If
                        omitted, called germline genes will not be annotated
  -outdir output_directory
                        output directory (default: working directory)
  -trim_primers path_to_primer_fasta
                        path to fasta file containing primers to be trimmed;
                        If omitted, no sequences will be trimmed
  --skip_preprocess     if set, preprocessing steps will be skipped
  --verbose             verbose output
  -max_log_v_evalue float_value
  -max_log_j_evalue float_value
  -min_v_sequence_length float_value
  -min_j_sequence_length float_value
  --allow_unproductive
  --allow_missing_cdr3
  --allow_ns_in_sequence
  -max_primer_error int_value
  --reverse_primers
  -lineage_clustering_cutoff float_value
  -hierarchical_clustering_method method
  -min_germline_usage_fraction float_value
  -min_germline_num_lineages float_value
  -cloud_radius float_value
  -min_usage_fraction_within_cloud float_value
  ```
