
##About

This repository contains a python package that calls immunoglobulin V germline alleles from deep BCR sequencing data.

##Requirements

See environment.yaml

##Usage
```
usage: call_ig_germlines.py [-h] [-annotate ANNOTATE] [-outdir OUTDIR]
                            [-trim_primers PRIMER_FASTA_PATH]
                            [--skip_preprocess] [--verbose]
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
  ```