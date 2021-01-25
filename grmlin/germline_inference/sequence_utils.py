import regex

complements = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}


def rc(seq):
    """ returns reverse complement of sequence """
    return "".join(complements.get(base, base) for base in reversed(seq))


def OR(RE_LIST, ENHANCEDMATCH=False):
    """ regex capture one or none in list of REs """
    if ENHANCEDMATCH:
        return "(?e)(" + "|".join(["(?:"+RE+")" for RE in RE_LIST]) + ")"
    else:
        return "(" + "|".join(["(?:"+RE+")" for RE in RE_LIST]) + ")"


def make_searcher(options, max_levenshtein_dist=1):
    """"""
    error_string = "(?e){e<=%d}" % max_levenshtein_dist
    checkers = [(v, regex.compile("("+k+")"+error_string))
                for k, v in options.items()]

    def searcher(match):
        for (v, c) in checkers:
            if c.search(match):
                return v, c.search(match).span()
        else:
            return None, (None, None)
    return searcher


def make_corrector(options, max_levenshtein_dist=1):
    """ """
    error_string = "{e<=%d}" % max_levenshtein_dist
    checkers = [(v, regex.compile("("+k+")"+error_string))
                for k, v in options.items()]

    def corrector(match):
        for (v, c) in checkers:
            if c.match(match):
                return v

    return corrector

def read_fasta(filename):
    """ Reads multiline fasta from file, and returns sequences as seq_id, sequence tuples"""
    read_list = []

    header = None
    for line in open(filename,'r'):
        line = line.rstrip()
        
        if line.startswith(">"):
            if header is not None:
                read_list.append((header,new_seq))
            header = line[1:]
            new_seq = ""
        else:
            new_seq += line.upper()
    read_list.append((header,new_seq))
    return read_list

def fasta_to_dict(filename):
    """ Reads multiline fasta from file, and returns sequences as seq_id, sequence dict"""

    return {k: v for k, v in read_fasta(filename)}

