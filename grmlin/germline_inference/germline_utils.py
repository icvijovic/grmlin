import pandas as pd
import numpy as np

from germline_inference.blast_utils import *
from germline_inference.cluster_vdj import *


def extract_multilineage_sequences(airr_df, sequence_column,
                                   retain_columns=[],
                                   clustering_method='average',
                                   clustering_cutoff=0.1):
    """Extracts sequences from airr dataframe that are associated
        with multiple putative VDJ recombination events (i.e. "lineages").

    Positional arguments:
        - airr_df: dataframe in AIRR format
        - sequence_column_name : column of airr_df with sequences

    Keyword arguments:
        - retain_columns: list of  columns with sequence information to retain
                          for further analysis (typically sequence family)
        - clustering_method: hierarchical clustering method to be used
                             in clustering sequences into lineages
        - clustering_cutoff: hierarchical clustering threshold

    Returns:
        - dataframe of unique sequences associated with multiple lineages,
          cdr3's and lineage counts associated with this dataframe.
        - array: all pairwise distances between cdr3's associated with more than
          one lineage
    """

    # compile cdr3 lengths if not already calculated
    if 'cdr3_length' in airr_df.columns:
        pass
    else:
        print("No cdr3_lenght column found in dataframe. Creating new column.")
        airr_df['cdr3_length'] = airr_df['cdr3'].str.len()

    groupby_cols = [sequence_column] + retain_columns + ['cdr3_length']
    unique_seq_df = airr_df.groupby(groupby_cols)['cdr3'].agg(list).to_frame().reset_index()

    unique_seq_df['num'] = unique_seq_df.cdr3.str.len()
    unique_seq_df['D'] = unique_seq_df.cdr3.apply(get_pairwise_distances)
    unique_seq_df['D_fractional'] = unique_seq_df.D / unique_seq_df.cdr3_length
    unique_seq_df['D_cutoff'] = unique_seq_df.cdr3_length * clustering_cutoff
    unique_seq_df['lineage_ids'] = unique_seq_df.apply(lambda x:
                                                       get_cluster_ids(x.D,
                                                                       x.D_cutoff,
                                                                       method=clustering_method),
                                                       axis=1)

    unique_seq_df['num_lineages'] = unique_seq_df['lineage_ids'].apply(np.unique).str.len()
    cumulative_lineages = np.zeros(unique_seq_df.shape[0])
    cumulative_lineages[1:] = np.cumsum(unique_seq_df['num_lineages'].values)[:-1]
    unique_seq_df['lineage_ids'] = unique_seq_df.lineage_ids + cumulative_lineages

    # compile distance distribution for sequences associated with multiple
    # cdr3's of the same length
    distance_distribution = unique_seq_df[unique_seq_df.num_lineages > 1].D_fractional.apply(list).sum()
    distance_distribution = np.asarray(distance_distribution, dtype=float)

    # finally, groupby sequence and additional info columns
    # retain list of cdr3's (separated by ",") associated with each lineage
    unique_seq_df.cdr3 = unique_seq_df.cdr3.apply(lambda x: np.asarray(x))

    unique_seq_df['lineage_cdr3s'] = unique_seq_df.apply(lambda x:
                                                         "|".join([",".join(x.cdr3[x.lineage_ids == y])
                                                                  for y in np.unique(x.lineage_ids)]),
                                                         axis=1)

    groupby_cols = retain_columns + [sequence_column]
    unique_seq_df = unique_seq_df.groupby(groupby_cols)[['lineage_cdr3s',
                                                         'num_lineages']].\
                                  agg({'lineage_cdr3s': "|".join,
                                       'num_lineages': sum})

    # retain subset that is associated with multiple lineages
    # and sort
    unique_seq_df = unique_seq_df[unique_seq_df.num_lineages > 1]
    unique_seq_df = unique_seq_df.sort_values(by='num_lineages',
                                              ascending=False).reset_index()

    return unique_seq_df, distance_distribution


def call_germline_genes(candidate_sequence_df,
                        MIN_FRACTION=0.001,
                        MIN_NUMBER=10,
                        CLOUD_RADIUS=3,
                        MIN_USAGE_FRAC_WITHIN_CLOUD=0.05,
                        verbose=False):
    """Calls germline genes from sequences in candidate_sequence_df """

    sequence_family_columns = [x for x in candidate_sequence_df.columns
                              if (('group' in x) or ('family' in x))]

    sequence_columns = [x for x in candidate_sequence_df.columns
                        if ('v_seq' in x)]

    if len(sequence_family_columns) != 1:
        raise KeyError("Cannot identify sequence family column in dataframe."
                    "Expected exactly one column with 'group' or 'family' in its name")
    else:
        sequence_family_column = sequence_family_columns[0]

    if len(sequence_columns) != 1:
        raise KeyError("Cannot identify sequence column in dataframe."
                    "Expected exactly one column with 'seq' in its name")
    else:
        sequence_column = sequence_columns[0]

    df = candidate_sequence_df[[sequence_family_column, sequence_column, 'num_lineages']]

    total_seqs = df.shape[0]
    total_lineages = candidate_sequence_df['num_lineages'].sum()

    CUTOFF = max(MIN_FRACTION * total_lineages, MIN_NUMBER)
    n_too_few_lineages = (df['num_lineages'] < CUTOFF).sum()

    if verbose:
        print(" {n_seq} candidate sequences are associated with a "
              "total of {n_total} lineages\n"
              "  discarding {n_small} sequences "
              "with fewer than {CUTOFF} lineages".format(
                n_seq = total_seqs,
                n_total=total_lineages,
                n_small=n_too_few_lineages,
                CUTOFF=CUTOFF))

    df = df[df['num_lineages'] > CUTOFF]
    accepted_sequences = []

    for family in df[sequence_family_column].unique():
        sequences = df[df[sequence_family_column] == family][sequence_column].values
        nlins = df[df[sequence_family_column] == family]['num_lineages'].values

        if len(sequences) > 0:
            D = get_pairwise_distances(sequences, squareform=True)
            # sequences are the cloud of another lineage if
            # - their Levenshtein distance is within CLOUD_RADIUS, and
            # - they are not associated with another, closer lineage
            within_cloud = ((D <= CLOUD_RADIUS)
                            * (D > 0)
                            * (D == np.min(D + 1000 * (D == 0))))

            # compile ratio of lineages associated with each sequence
            # and with potential donor sequence
            usage_ratio = ((within_cloud*nlins).T/nlins).T

            # if ratio is smaller than
            # the minimal fraction allowed for sequence
            # to be considered a separate germline,
            # then consider this sequence a repeated hypermutation event
            likely_convergent_hypermut = np.any((usage_ratio > 0)
                                                * (usage_ratio < MIN_USAGE_FRAC_WITHIN_CLOUD),
                                                axis=1)

            accepted_sequences.extend(sequences[~likely_convergent_hypermut])

    df = df[df[sequence_column].apply(lambda x: x in accepted_sequences)]

    if verbose:
        print("  discarding {n} sequences as they are in the vicinity"
              " of more commonly generated sequences".format(
                n=total_seqs - n_too_few_lineages - df.shape[0])
              )

    return df


def annotate_germline_calls(df, path_to_database):
    """Compares germline calls to alleles in database
       Calls mutations with respect to this database.

       Returns: dataframe with new columns corresponding to:
          - sequence_id of best match in database
          - mutations with respect to that match
          - number of mismatches with respect to that sequence
    """

    sequence_columns = [x for x in df.columns
                        if ('v_seq' in x)]

    if len(sequence_columns) != 1:
        raise KeyError("Cannot identify sequence column in dataframe."
                    "Expected exactly one column with 'seq' in its name")
    else:
        sequence_column = sequence_columns[0]

    # check if dataframe is uniquely indexed
    # this index will be passed to blastn
    if df.index.is_unique:
        pass
    else:
        print(" Waring! Dataframe is not uniquely indexed. Resetting index.")
        df = df.reset_index()

    df['sequence_id'] = df.index.astype(str) \
                      + "|n_lineages=" \
                      + df.num_lineages.astype(str)

    fasta_query_string = "\n".join((">" + df.sequence_id + "\n" + df[sequence_column]).values)

    blastn_out_df, blastn_err = pipe_to_blastn(fasta_query_string, path_to_database)
    blastn_matches = return_best_match(blastn_out_df)
    
    # ensure btop has not been accidentally cast as int
    blastn_matches['btop'] = blastn_matches['btop'].astype(str)
    blastn_matches['mutations'] = blastn_matches.apply(lambda x:
                                                       parse_btop(x.btop, x.sstart),
                                                       axis=1)

    blastn_matches = blastn_matches[['sequence_id', 'match', 'mismatch', 'mutations']]
    df = df.merge(blastn_matches, on='sequence_id', how='outer').drop('sequence_id', axis=1)

    return df


if __name__ == '__main__':

    print('Success!')
