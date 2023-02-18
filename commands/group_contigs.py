import logging
import pickle
from itertools import groupby
import click


def load_pickled_data(pickle_file):
    with open(pickle_file) as f:
        cosine_matrix = pickle.load(f)
        matrix_key = pickle.load(f)

    return cosine_matrix, matrix_key


def group_by_cosine(cosine_array, labels, threshold):
    def filter_by_cosine(array, dictionary):
        valid_contigs = []
        for n, cosine in enumerate(array):
            if cosine > threshold:
                valid_contigs.append(dictionary[n])

        return valid_contigs

    logging.info('grouping contigs...')

    normalised_keys = {n: k for n, k in enumerate(labels)}

    array_dict = {normalised_keys[n]: array for n, array in enumerate(cosine_array)}

    matched_contig = [filter_by_cosine(v, normalised_keys) for k, v in array_dict.items()]

    non_duplicate_contig_groups = list(matched_contig for matched_contig, _ in groupby(matched_contig))

    logging.info('writing data to pickle...')
    with open(f'grouped_contigs_({threshold}).pkl', 'wb') as f:
        pickle.dump(non_duplicate_contig_groups, f)


@click.command()
@click.option('-t', '--threshold', type=float, required=True, help="cosine similarity threshold for grouping")
@click.option('-i', '--input_data', type=str, required=True, help="pickle file generated with calculate_cosine")
def group_contigs(threshold, input_data):
    cosine, key = load_pickled_data(input_data)
    group_by_cosine(cosine, key, threshold)
