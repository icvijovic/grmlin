
## About

This repository contains a python package that calls immunoglobulin V germline alleles from deep BCR sequencing data.

## Requirements

See environment.yaml

## Usage
```
usage: call_ig_germlines.py [-h] [-annotate ANNOTATE] [-outdir OUTDIR]
                            [-trim_primers PRIMER_FASTA_PATH]
                            [--skip_preprocess] [--verbose]
                            [-max_log_v_evalue MAX_LOG_V_EVALUE]
                            [-max_log_j_evalue MAX_LOG_J_EVALUE]
                            [--allow_unproductive] [--allow_missing_cdr3]
                            [--allow_ns_in_sequence]
                            [-min_v_sequence_length MIN_V_SEQUENCE_LENGTH]
                            [-min_j_sequence_length MIN_J_SEQUENCE_LENGTH]
                            [-max_primer_error MAX_PRIMER_ERROR]
                            [--reverse_primer]
                            [-lineage_clustering_cutoff LINEAGE_CLUSTERING_CUTOFF]
                            [-hierarchical_clustering_method HIERARCHICAL_CLUSTERING_METHOD]
                            [-min_germline_usage_fraction MIN_GERMLINE_USAGE_FRACTION]
                            [-min_germline_num_lineages MIN_GERMLINE_NUM_LINEAGES]
                            [-cloud_radius CLOUD_RADIUS]
                            [-min_usage_fraction_within_cloud MIN_USAGE_FRACTION_WITHIN_CLOUD]
                            input_file_path

Constructs a database of Ig V germlines.

positional arguments:
  input_file_path       path to input AIRR-formatted table containing Ig
                        sequences

optional arguments:
  -h, --help            show this help message and exit
  -annotate ANNOTATE    path to database containing known Ig V alleles if
                        omitted, called germline genes will not be annotated
  -outdir OUTDIR        output directory (default: working directory)
  -trim_primers PRIMER_FASTA_PATH
                        path to fasta file containing primers to be trimmed if
                        omitted, no sequences will be trimmed
  --skip_preprocess     if set, preprocessing steps will be skipped
  --verbose             verbose output
  -max_log_v_evalue MAX_LOG_V_EVALUE
  -max_log_j_evalue MAX_LOG_J_EVALUE
  --allow_unproductive
  --allow_missing_cdr3
  --allow_ns_in_sequence
  -min_v_sequence_length MIN_V_SEQUENCE_LENGTH
  -min_j_sequence_length MIN_J_SEQUENCE_LENGTH
  -max_primer_error MAX_PRIMER_ERROR
  --reverse_primer
  -lineage_clustering_cutoff LINEAGE_CLUSTERING_CUTOFF
  -hierarchical_clustering_method HIERARCHICAL_CLUSTERING_METHOD
  -min_germline_usage_fraction MIN_GERMLINE_USAGE_FRACTION
  -min_germline_num_lineages MIN_GERMLINE_NUM_LINEAGES
  -cloud_radius CLOUD_RADIUS
  -min_usage_fraction_within_cloud MIN_USAGE_FRACTION_WITHIN_CLOUD
  ```