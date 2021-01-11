import numpy as np
import pandas as pd

import argparse
import time

from clean_vdj import *
from germline_utils import *
import sequence_utils


parser = argparse.ArgumentParser(description='Calls Ig germlines.')
parser.add_argument('input_airr', metavar='input_file_path', type=str,
                    help='path to input airr file')
parser.add_argument('-outdir', type=str, default='./outs',
                    help='output directory '
                         '(default: ./outs)')
parser.add_argument('-trim_v_primers', metavar='v_primer_fasta_path',
                    type=str, required=False)
parser.add_argument('-trim_j_primers', metavar='j_primer_fasta_path',
                    type=str, required=False)
parser.add_argument('--from_cache', required=False, default=False,
                    dest='from_cache', action='store_true')

args = parser.parse_args()

# parse i/o paths
INPUT_FILENAME = args.input_airr
OUTDIR = args.outdir

if args.trim_v_primers is None:
    trim_v_primers = False
else:
    trim_v_primers = True
    v_primer_dict = sequence_utils.fasta_to_dict(args.trim_v_primers)

if args.trim_j_primers is None:
    trim_j_primers = False
else:
    trim_j_primers = True
    j_primer_dict = sequence_utils.fasta_to_dict(args.trim_j_primers)

from_cache = args.from_cache

if __name__ == '__main__':

    df = pd.read_csv(INPUT_FILENAME, sep='\t').dropna(how='all', axis=1)

    if not from_cache:
        # ensure that all sequences pass QC criteria
        df = clean_vdj_dataframe(df,
                                 MAX_LOG_V_EVALUE=-60,
                                 MAX_LOG_J_EVALUE=-10,
                                 ALLOW_UNPRODUCTIVE=False,
                                 ALLOW_MISSING_CDR3=False,
                                 ALLOW_Ns_IN_SEQUENCE=False,
                                 MIN_V_SEQUENCE_LENGTH=160,
                                 MIN_J_SEQUENCE_LENGTH=20,
                                 verbose=True)

        # parse out v and j templated sequences
        df['v_sequence'] = df.apply(lambda x: x.sequence[int(x.v_sequence_start):
                                                         int(x.cdr3_start)],
                                    axis=1)
        df['j_sequence'] = df.apply(lambda x: x.sequence[int(x.j_sequence_start):
                                                         int(x.j_sequence_end)],
                                    axis=1)

        v_templated_col, j_templated_col = 'v_sequence', 'j_sequence'

        if trim_v_primers:
            # v primers assumed in fwd orientation
            df = trim_primer_sequence(df, v_templated_col, v_primer_dict,
                                      primer_orientation='fwd', max_error=2)
            v_templated_col = v_templated_col + "_trimmed"

        if trim_j_primers:
            # j primers assumed in rev orientation
            df = trim_primer_sequence(df, j_templated_col, j_primer_dict,
                                      primer_orientation='rev', max_error=2)
            j_templated_col = j_templated_col + "_trimmed"

        # extract v family info from v_call
        df['v_family'] = df.v_call.map(lambda x: x.split("-")[0].split("/")[0])

        # collapse truncated sequences
        v_seq_trunc_map = {}

        # for efficiency, only check whether v sequences within the same family
        # represent truncations of one another
        print("constructing truncation map...")
        for family in df.v_family.unique():
            family_truncation_map = construct_truncation_map(df[df.v_family == family][v_templated_col])
            v_seq_trunc_map.update(family_truncation_map)
        df[v_templated_col+"_no_trunc"] = df[v_templated_col].map(v_seq_trunc_map)
        v_templated_col = v_templated_col + "_no_trunc"  

        # cache dataframe

        df = df.to_csv('germline_cache/{}_preprocessed.tsv.gz'.format(INPUT_FILENAME.split("/")[-1].split(".")[0]), sep = '\t')               
    else:
        v_templated_col = "v_sequence_no_trunc"
        j_templated_col = "j_sequence"
    print(df.columns)

    # extract multilineage v sequence info
    multilin_vs, cdr3_dists_v = extract_multilineage_sequences(df,
                                                             v_templated_col,
                                                             retain_columns=['v_family'],
                                                             clustering_method='average',
                                                             clustering_cutoff=0.1)
    # call v germlines
    v_germline_df = call_germline_genes(multilin_vs,
                                        MIN_FRACTION=0.001,
                                        MIN_NUMBER=10,
                                        CLOUD_RADIUS=3,
                                        MIN_USAGE_FRAC_WITHIN_CLOUD=0.05,
                                        verbose=True)
    # annotate v germlines
    IMGT_V_DB = '../../../Sticklebacks/data/test_human_data/imgt_dbs/imgt_human_ig_v.fasta'
    v_germline_df = annotate_germline_calls(v_germline_df, IMGT_V_DB)

    # to infer j germlines, focus only on subset of sequences that
    # have no v hypermutations
    germline_v_sequences = df[v_templated_col].map(lambda x: 
                                                   x in v_germline_df[v_templated_col].values)
    germline_subset = df[germline_v_sequences].copy()

    print(v_germline_df[['match','num_lineages','mutations']])

    # extract j group info from j_call
    germline_subset['j_group'] = germline_subset.j_call.map(lambda x:
                                  "".join([y.split("*")[0] for y in x.split(";")]))
    germline_subset['j_group'] = germline_subset.j_group.map(human_j_gene_groups)

    print("constructing j truncation map...")
    j_seq_trunc_map = {}
    for group in germline_subset.j_group.unique():
        group_subset = germline_subset[germline_subset.j_group == group]
        group_truncation_map = construct_truncation_map(group_subset[j_templated_col])
        j_seq_trunc_map.update(group_truncation_map)

    germline_subset[j_templated_col+"_no_trunc"] = germline_subset[j_templated_col].map(j_seq_trunc_map)
    j_templated_col = j_templated_col + "_no_trunc"

    multilin_js, cdr3_dists_j = extract_multilineage_sequences(germline_subset,
                                                             j_templated_col,
                                                             retain_columns=['j_group'],
                                                             clustering_method='average',
                                                             clustering_cutoff=0.1)

    # call j germlines
    j_germline_df = call_germline_genes(multilin_js,
                                        MIN_FRACTION=0.001,
                                        MIN_NUMBER=5,
                                        CLOUD_RADIUS=3,
                                        MIN_USAGE_FRAC_WITHIN_CLOUD=0.1,
                                        verbose=True)
     # annotate j germlines
    IMGT_J_DB = '../../../Sticklebacks/data/test_human_data/imgt_dbs/imgt_human_ig_j.fasta'
    j_germline_df = annotate_germline_calls(j_germline_df, IMGT_J_DB)


    print(j_germline_df[['match','num_lineages','mutations']].groupby(['match'])['mutations'].agg(";".join))

