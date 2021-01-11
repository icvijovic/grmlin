import pandas as pd
import numpy as np

import regex
import germline_inference.sequence_utils as sequence_utils


def clean_vdj_dataframe(airr_df,
                        MAX_LOG_V_EVALUE=-5,
                        MAX_LOG_J_EVALUE=-5,
                        ALLOW_UNPRODUCTIVE=False,
                        ALLOW_MISSING_CDR3=False,
                        ALLOW_Ns_IN_SEQUENCE=False,
                        MIN_V_SEQUENCE_LENGTH=160,
                        MIN_J_SEQUENCE_LENGTH=20,
                        verbose=False):
    """Clean AIRR dataframe to retain only high-quality VDJ sequences.

    Positional arguments:
    - airr_df: dataframe in AIRR format.

    Keyword arguments:
    - MAX_LOG_V_EVALUE: max natural log E-value for V gene alignment
    - MAX_LOG_J_EVALUE: max natural log E-value for J gene alignment
    - ALLOW_UNPRODUCITVE: retain sequences with an in-frame stop codon
    - ALLOW_MISSING_CDR3: retain sequences that do not contain a CDR3
    - ALLOW_Ns_IN_SEQUENCE: retain sequences that contain N's
    - MIN_V_SEQUENCE_LENGTH: shortest allowed V gene alignment
    - MIN_J_SEQUENCE_LENGTH: shortest allowed J gene alignment
    - verbose: print filtering log to stdout

    Returns:
    - filtered dataframe containing rows satisfying QC criteria

    """
    df = airr_df.copy()

    # verify that dataframe contains anticipated columns
    for col in ['sequence',
                'v_support',
                'j_support',
                'productive',
                'cdr3',
                'v_sequence_start',
                'v_sequence_end',
                'j_sequence_start',
                'j_sequence_end']:
        if not(col in df.columns):
            raise KeyError('The following columns were not found'
                           'in dataframe: {col}'.format(col=col))
    if verbose:
        print(" Input dataframe contains {n} sequences".format(n=df.shape[0]))

    good_v_map = np.log(df.v_support.astype(float)) < MAX_LOG_V_EVALUE
    good_j_map = np.log(df.j_support.astype(float)) < MAX_LOG_J_EVALUE
    good_v_len = (df.v_sequence_end - df.v_sequence_start + 1) > MIN_V_SEQUENCE_LENGTH
    good_j_len = (df.j_sequence_end - df.j_sequence_start + 1) > MIN_J_SEQUENCE_LENGTH

    if verbose:
        if not good_v_map.all():
            print("   {n} sequences discarded because of poor V gene support.".format(
                   n=((~good_v_map).sum())))
        else:
            print("   all sequences have good V gene support.")

        if not good_j_map.all():
            print("   {n} sequences discarded because of poor J gene support/".format(
                    n=((~good_j_map).sum())))
        else:
            print("   all sequences have good J gene support.")

        if not good_v_len.all():
            print("   {n} sequences discarded because their V gene alignment is too short.".format(
                    n=((~good_v_len).sum())))
        else:
            print("   all sequences have a long enough V gene alignment.")

        if not good_j_len.all():
            print("   {n} sequences discarded because their J gene alignment is too short.".format(
                    n=((~good_j_len).sum())))
        else:
            print("   all sequences have a long enough J gene alignment.")

    df = df[good_v_map & good_j_map & good_j_len & good_j_len]

    if ALLOW_UNPRODUCTIVE:
        pass
    else:
        productive = df.productive == 'T'
        if verbose:
            if not productive.all():
                print("   {n} sequences discarded because they are not productive.".format(
                      n=((~productive).sum())))
            else:
                print("   all sequences appear productive.")

        df = df[productive]

    if ALLOW_MISSING_CDR3:
        pass
    else:
        has_cdr3 = df.cdr3.notna()

        if verbose:
            if not has_cdr3.all():
                print("   {n} sequences discarded because they do not contain a CDR3.".format(
                    n=((~has_cdr3).sum())))
            else:
                print("   all sequences appear to have a CDR3.")
        df = df[has_cdr3]

    if ALLOW_Ns_IN_SEQUENCE:
        pass
    else:
        N_in_sequence = df.sequence.map(lambda x: "N" in x)

        if verbose:
            if N_in_sequence.any():
                print("   {n} sequences discarded because they contain an N.".format(
                     n=((N_in_sequence).sum())))
            else:
                print("   all sequences have unambiguous bases.")

        df = df[~N_in_sequence]

    if verbose:
        print(" Filtered dataframe contains {n} sequences.".format(
                n=df.shape[0]))

    return df


# TODO: ensure mapping is transitive
def construct_truncation_map(sequences, MAX_LEN_DIFF_OF_SHORTER_SEQ=6):
    """
    Map sequences that represent truncations of one another
    to the more common of each pair.

    Positional arguments:
    - sequences:  array-like, Iterable, or dict of sequences

    Keyword arguments:
    - MAX_LEN_DIFF_OF_SHORTER_SEQ : maximal number of nucleotides
                                    that can be removed in truncation map

    Returns:
    - dict of sequence : mapped_sequence pairs

    """
    if type(sequences) == pd.Series:
        pass
    else:
        sequences = pd.Series(sequences)

    # sort sequences by number of occurrences
    unique_sequences = sequences.value_counts().index.values

    sequence_map = {}

    for it1 in reversed(range(len(unique_sequences))):
        # look through all sequences more common than this one
        for it2 in range(it1+1):
            # if either is a truncation of the other,
            # then map less common sequence to more common one
            if len(unique_sequences[it1]) < len(unique_sequences[it2]):
                if (unique_sequences[it1] in unique_sequences[it2]):
                    sequence_map[unique_sequences[it1]] = unique_sequences[it2]
                    break
                else:
                    pass
            elif (unique_sequences[it2] in unique_sequences[it1]):
                # only map sequence to a shorter more common sequence if they
                # differ by fewer than MAX_LEN_DIFF_OF_SHORTER_SEQ nucleotides
                if len(unique_sequences[it1]) - len(unique_sequences[it2]) <= MAX_LEN_DIFF_OF_SHORTER_SEQ:
                    sequence_map[unique_sequences[it1]] = unique_sequences[it2]
                    break

    return sequence_map


# TODO: write rev primer implementation
def trim_primer_sequence(airr_df, seq_colname, primers,
                         primer_orientation='fwd', max_error=2):
    """ Trims primers from sequences found in specified column of dataframe.

    Positional arguments:
    - airr_df : dataframe containing sequences to be trimmed
    - seq_colname : column in which sequences are to be trimmed
    - primers : dict of primer sequences

    Keyword arguments:
    - primer_orientation : 'fwd' or 'rev' with respect to sequences in seq_colname
    - max_error : max Levenshtein distance for match with primer sequence

    Returns:
    -  copy of dataframe with two additional columns:
    - "{seq_colname}_trimmed": trimmed sequences
    - "{seq_colname}_primers": names of primers found

    """
    df = airr_df.copy()

    primer_lookup_dict = {v: k for k, v in primers.items()}

    if primer_orientation == 'fwd':
        sequences = df[seq_colname]
    elif primer_orientation == 'rev':
        sequences = df[seq_colname].apply(sequence_utils.rc)
    else:
        raise ValueError("Unknown primer orientation: {}. "
                         "Expected 'fwd' or 'rev'.".format(primer_orientation))

    _primer_parser = sequence_utils.make_searcher(primer_lookup_dict,
                                            max_levenshtein_dist=max_error)

    primer_search_results = sequences.apply(_primer_parser)

    def _trim(sequence, pspan):
        if pspan is (None, None):
            return sequence
        else:
            if primer_orientation == 'fwd':
                return sequence[pspan[1]:]
            else:
                return sequence[pspan[0]]
    
    df[seq_colname+"_primers"] = primer_search_results.map(lambda x: x[0])
    df["pspan"] = primer_search_results.map(lambda x: x[1])

    df[seq_colname+"_trimmed"] = df.apply(lambda row: _trim(row[seq_colname],
                                                            row.pspan),
                                          axis=1)
    df = df.drop(["pspan"], axis = 1)

    return df
