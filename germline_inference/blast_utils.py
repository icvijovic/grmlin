import regex
import numpy as np
import pandas as pd

from subprocess import Popen, PIPE
from io import BytesIO


DEFAULT_OUTPUT_OPTIONS = "qseqid sseqid pident length mismatch " \
                       + "gapopen qstart qend sstart send " \
                       + "evalue btop "


def pipe_to_blastn(query_string,
                   path_to_db,
                   output_options=DEFAULT_OUTPUT_OPTIONS):

    cmd = 'blastn -query <(echo -e \"{}\")'.format(query_string)
    cmd += ' -db ' + path_to_db
    cmd += ' -evalue 1e-05 '
    cmd += ' -outfmt "6 {out_options}"'.format(out_options=output_options)

    # print(cmd)
    p = Popen(cmd, stdout=PIPE, stderr=PIPE, shell=True, executable='/bin/bash')
    out, err = p.communicate()

    _col_names = {it: column for it, column in enumerate(output_options.split(" "))}
    _col_names.update({0: 'sequence_id', 1: 'match'})

    if len(out) > 0:
        out = BytesIO(out)
        out_df = pd.read_csv(out, sep='\t', header = None)
        out_df = out_df.rename(columns=_col_names)
    else:
       out_df =  pd.DataFrame({col: [] for col in _col_names.values()})

    return out_df, err


def return_best_match(blastn_out_df):

    if len(blastn_out_df) > 0:

        blastn_out_df['evalue'] = blastn_out_df['evalue'].astype(float)
        blastn_out_df['gapopen'] = blastn_out_df['gapopen'].astype(int)
        blastn_out_df['qstart'] = blastn_out_df['qstart'].astype(int)
        blastn_out_df['qend'] = blastn_out_df['qend'].astype(int)
        blastn_out_df['sstart'] = blastn_out_df['sstart'].astype(int)

        blastn_out_df['send'] = blastn_out_df['send'].astype(int)
        blastn_out_df['pident'] = blastn_out_df['pident'].astype(float)

        match_key = blastn_out_df.loc[blastn_out_df.groupby(['sequence_id'])['evalue'].transform('min') == blastn_out_df['evalue']]
        match_key = match_key.sort_values(['evalue', 'match']).drop_duplicates('sequence_id')

        return match_key

    else:
        print("Warning: Provided empty dataframe\n")
        return blastn_out_df


def parse_btop(btop, sstart=0, return_string = True):
    """ Parse btop """
    debug = False
    #######################
    # From the NCBI documentation:
    # a number represents this many matches,
    # two bases represent a mismatch and show query and reference base,
    # - base and gap or gap and base, show a gap in query or reference,
    # _<number>_ represents an insertion (gap in reference) of this number of bases,
    # %<number>% represents a deletion (gap in read) of this number of bases,
    ########################
    short_variant_strings = ['A', 'C', 'T', 'G', 'N']
    indel_strings = ['_', "%"]
    indel_string_type_dict = {"_": "+", "%": "-"}

    # in keeping with NCBI conventions, locations are 1-indexed
    sposition = sstart
    mutations = {}

    btop_pos = 0

    previous = -1

    while btop_pos + 1 < len(btop):

        if btop[btop_pos].isdigit():
            # this number indicates a match
            # read match length and move to next item
            digits = regex.split('(\d+)', btop[btop_pos:])[1]
            matches = int(digits)
            sposition += matches
            btop_pos += len(digits)
            if debug:
                print("match", matches)
            continue

        if btop[btop_pos] in short_variant_strings:
            # this base was not present in reference sequence
            # it's either a snp or short insertion
            if btop[btop_pos + 1] in short_variant_strings:
                # this is a snp
                derived = btop[btop_pos]
                ref = btop[btop_pos + 1]
                # record it and move to next item
                mut_string = "{ref}>{derived}".format(ref=ref, derived=derived)
                if mutations.get(sposition, None) is None:
                    mutations[sposition] = mut_string
                else:
                    mutations[sposition] = ','.join([mutations[sposition], mut_string])
                if debug:
                    print("snv", mut_string)
                btop_pos += 2
                sposition += 1

            elif btop[btop_pos + 1] == '-':
                # this is an insertion
                k = btop_pos + 1

                # find lenght of inserted sequence
                while True:
                    if (len(btop) <= k + 1):
                        break
                    if (btop[k+1].isdigit()):
                        break
                    if (len(btop) <= k + 2):
                        break
                    if (btop[k+2] == '-'):
                        k += 2
                    else:
                        break

                inserted_sequence = btop[btop_pos:k+1][::2]
                # record insertion and move to next item
                mut_string = "+{s}".format(s=inserted_sequence)
                if mutations.get(sposition, None) is None:
                    mutations[sposition] = mut_string
                else:
                    mutations[sposition] = ','.join([mutations[sposition], mut_string])
                if debug:
                    print(btop[btop_pos:][:20] + "...", end='   ')
                    print('ins', mut_string)
                btop_pos = k + 1

                # do not change sposition when documenting an insertion
            continue

        if btop[btop_pos] == '-':
            # this is a deletion from the reference sequence

            # find length of deleted sequence
            k = btop_pos
            while btop[k] == '-':
                k += 2
            deleted_sequence = btop[btop_pos:k+1][1::2]

            # record insertion and move to next item
            mut_string = "-{s}".format(s=deleted_sequence)
            if debug:
                print(btop[btop_pos:][:20] + "...", end='   ')
                print('del', mut_string)
            if mutations.get(sposition, None) is None:
                mutations[sposition] = mut_string
            else:
                mutations[sposition] = ','.join([mutations[sposition], mut_string])
            btop_pos = k
            sposition += len(deleted_sequence)
            continue

        if btop[btop_pos] in ['%', '_']:
            # this is a larger indel

            # find length of indel sequence
            k = btop_pos + 1
            while btop[k] != btop[btop_pos]:
                k += 1
            indel_length = int(btop[btop_pos:k])

            # record variant and move to next item
            mut_string = indel_string_type_dict[btop[btop_pos]] + str(indel_length)
            if mutations.get(sposition, None) is None:
                mutations[sposition] = mut_string
            else:
                mutations[sposition] = ','.join([mutations[sposition], mut_string])
            btop_pos = k
            if indel_string_type_dict[btop[btop_pos]] == '-':
                sposition += len(deleted_sequence)
            continue

    if return_string:
         return(",".join([str(x) + ":" + y for x, y in mutations.items()]))
    else:
        # returns dict rather than string:
        return(mutations)
