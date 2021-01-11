import numpy as np
import pandas as pd

import argparse

from germline_inference.clean_vdj import *
from germline_inference.germline_utils import *
import germline_inference.sequence_utils as sequence_utils


parser = argparse.ArgumentParser(description='Constructs a database of Ig V germlines.')

parser.add_argument('input_airr', metavar='input_file_path', type=str,
                    help='path to input AIRR-formatted table containing Ig sequences')

parser.add_argument('-annotate', type=str, required=False,
                    help='path to database containing known Ig V alleles\n'
                    'if omitted, called germline genes will not be annotated')

parser.add_argument('-outdir', type=str, default='./',
                    help='output directory '
                         '(default: working directory)')

parser.add_argument('-trim_primers', metavar='PRIMER_FASTA_PATH',
                    type=str, required=False,
                    help='path to fasta file containing primers to be trimmed\n'
                    'if omitted, no sequences will be trimmed')

parser.add_argument('--skip_preprocess', required=False, default=False,
                    dest='skip_preprocess', action='store_true',
                    help='if set, preprocessing steps will be skipped')

parser.add_argument('--verbose', required=False, default=False,
                    dest='verbose', action='store_true',
                    help='verbose output')

args = parser.parse_args()

# parse i/o paths
INPUT_FILENAME = args.input_airr
outdir = args.outdir

if args.trim_primers is None:
    trim_primers = False
else:
    trim_primers = True
    primer_dict = sequence_utils.fasta_to_dict(args.trim_primers)

skip_preprocess = args.skip_preprocess
verbose = args.verbose

sample_name = INPUT_FILENAME.split("/")[-1].split(".tsv.gz")[0]

if __name__ == '__main__':

    if args.annotate is None:
        print("Warning: Existing Ig V allele database not provided. "
              "Called germline alleles will not be annotated.")

    df = pd.read_csv(INPUT_FILENAME, sep='\t').dropna(how='all', axis=1)
    
    if not skip_preprocess:
        # ensure that all sequences pass QC criteria
        df = clean_vdj_dataframe(df,
                                 MAX_LOG_V_EVALUE=-60,
                                 MAX_LOG_J_EVALUE=-10,
                                 ALLOW_UNPRODUCTIVE=False,
                                 ALLOW_MISSING_CDR3=False,
                                 ALLOW_Ns_IN_SEQUENCE=False,
                                 MIN_V_SEQUENCE_LENGTH=160,
                                 MIN_J_SEQUENCE_LENGTH=20,
                                 verbose=verbose)

        # parse out templated v sequence
        df['v_sequence'] = df.apply(lambda x: x.sequence[int(x.v_sequence_start):
                                                         int(x.cdr3_start)],
                                    axis=1)

        v_templated_col = 'v_sequence'

        if trim_primers:
            # v primers assumed in fwd orientation
            if verbose:
                print(" Trimming primers...")
            df = trim_primer_sequence(df, v_templated_col, primer_dict,
                                      primer_orientation='fwd', max_error=2)
            v_templated_col = v_templated_col + "_trimmed"
            if verbose:
                print(" Done!")

        # extract v family info from v_call
        df['v_family'] = df.v_call.map(lambda x: x.split("-")[0].split("/")[0])

        # collapse truncated sequences
        v_seq_trunc_map = {}

        # for efficiency, only check whether v sequences within the same family
        # represent truncations of one another
        if verbose:
            print(" Collapsing sequences that represent truncations of one another...")
        for family in df.v_family.unique():
            family_truncation_map = construct_truncation_map(df[df.v_family == family][v_templated_col])
            v_seq_trunc_map.update(family_truncation_map)
        if verbose: 
            print(" Done!")

        df[v_templated_col+"_no_trunc"] = df[v_templated_col].map(v_seq_trunc_map)
        v_templated_col = v_templated_col + "_no_trunc"  

        # cache dataframe
        cache_dest = '{}/{}_preprocessed.tsv.gz'.format(outdir, sample_name)
        if verbose:
            print(" Caching cleaned AIRR dataframe, see {}".format(cache_dest))

        df.to_csv(cache_dest, sep = '\t')               
        
    else:
        v_templated_col = [x for x in df.columns if x.endswith("no_trunc")]
        if len(v_templated_col) != 1:
            raise KeyError("Failed to identify column containing cleaned sequences"
                           "Expected exactly one column of dataframe to end with string 'no_trunc'")
        v_templated_col = v_templated_col[0]

        # avoid renaming sample
        if sample_name.endswith("preprocessed"):
            sample_name = sample_name.split("_preprocessed")[0]

    # extract multilineage v sequence info
    multilin_vs, cdr3_dists_v = extract_multilineage_sequences(df,
                                                             v_templated_col,
                                                             retain_columns=['v_family'],
                                                             clustering_method='average',
                                                             clustering_cutoff=0.1)
    # cache multilineage sequences and cdr3 distribution
    multinlin_v_dest = "{}/{}_multilineage_vs.tsv.gz".format(outdir, sample_name)
    cdr3_dists_v_dest = "{}/{}_cdr3_distance_distribution.tsv.gz".format(outdir, sample_name)
   
    if verbose:
        print(" Caching information about all V sequences associated with"
              "multiple lineages, see:\n{}\n{}".format(multinlin_v_dest, cdr3_dists_v_dest))

    multilin_vs.to_csv(multinlin_v_dest, sep='\t')
    pd.Series(cdr3_dists_v).to_csv(cdr3_dists_v_dest, sep='\t')

    # call v germlines
    v_germline_df = call_germline_genes(multilin_vs,
                                        MIN_FRACTION=0.001,
                                        MIN_NUMBER=10,
                                        CLOUD_RADIUS=3,
                                        MIN_USAGE_FRAC_WITHIN_CLOUD=0.05,
                                        verbose=verbose)
    # annotate v germlines
    if args.annotate is None:
        pass
    else:
        v_germline_df = annotate_germline_calls(v_germline_df, args.annotate)
    
    if verbose:
      print(" Called a total of {n} germline alleles. ".format(n = v_germline_df.shape[0]))
      if args.annotate is not None:
        print(" {m} of these can be mapped to alleles in provided database.".format(
              m = v_germline_df.match.notna().sum()))

    # print germline_calls
    germline_dest = "{}/{}_v_germlines.tsv".format(outdir, sample_name)

    print(" Saving germline allele information to: {}".format(germline_dest))

    v_germline_df.to_csv(germline_dest, sep = '\t')
