import numpy as np
import pandas as pd

from scipy.cluster.hierarchy import single, average, linkage, fcluster
from Levenshtein import distance


def get_pairwise_distances(seq_list, squareform=False):

    n = len(seq_list)

    distances = np.zeros((n, n))

    for i in range(0, n):
        for j in range(0, i):
            distances[i, j] = distance(seq_list[i], seq_list[j])
            distances[j, i] = distances[i, j]

    if squareform:
        return distances
    else:
        upper_triangular_indices = np.triu_indices(n, 1)
        condensed_distance_matrix = distances[upper_triangular_indices]

        return condensed_distance_matrix


def get_cluster_ids(condensed_distance_matrix, cutoff, method='average'):

    if len(condensed_distance_matrix) == 0:
        return np.asarray([0])
    else:
        if method == 'average':
            Z = average(condensed_distance_matrix)
        elif method == 'single':
            Z = single(condensed_distance_matrix)
        else:
            Z = linkage(condensed_distance_matrix, method=method)

        return fcluster(Z, cutoff, criterion='distance')
