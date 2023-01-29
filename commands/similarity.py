import csv
import itertools
import pickle

import click
from sklearn.metrics.pairwise import cosine_similarity
from sklearn.feature_extraction.text import CountVectorizer


def extract_mapping_data(sam_file):
    contig_read_dictionary = {}
    print("extracting infromation from SAM file...")
    with open(sam_file) as file:
        for line in file:
            if line.startswith("@SQ"):
                seq_name = line.split("\t")[1][3:]
                contig_read_dictionary[seq_name] = []

            else:
                read_info = line.split("\t")
                query_name = read_info[0]
                reference_name = read_info[2]
                if reference_name == "*":
                    continue
                else:
                    contig_read_dictionary[reference_name].append(query_name)

    return contig_read_dictionary


def convert_reads_to_hexidemical(contig_read_dictionary, index):
    print("extracting infromation from index file...")
    with open(index, mode='r') as inp:
        index_reader = csv.reader(inp)
        ID_colour_dict = {rows[0]: rows[1] for rows in index_reader}

    print("converting Seq IDs to RGB values...")
    contig_rgb = {contig: list(map(lambda x: ID_colour_dict[x], contig_read_dictionary[contig])) for contig in
                  contig_read_dictionary}

    print("converting RGB values to hexidecimal values...")
    contig_hex_dictionary = {
        contig: list(map(lambda x: '#%02x%02x%02x' % tuple(map(int, x.split(","))), contig_rgb[contig])) for
        contig in contig_rgb}

    return contig_hex_dictionary


def calculate_cosine_similarity(contig_hexidemical_dictionary):
    print('converting profiles to text...')
    contig_adapter_profiles = [" ".join(contig_hexidemical_dictionary[contig]) for contig in
                               contig_hexidemical_dictionary]

    labels = contig_hexidemical_dictionary.keys()

    print('vectorizing profiles...')
    count_array = CountVectorizer(max_df=0.2)
    profile_count_array = count_array.fit_transform(contig_adapter_profiles)

    print('calculating cosine similarity...')
    cosine_array = cosine_similarity(profile_count_array, dense_output=False)

    cosine_dictionary = dict(cosine_array.todok())

    matched_contigs = {int(val): [] for val in range(1, len(labels) + 1)}

    for i in cosine_dictionary:
        key = str(i)[1:-1].split(",")
        reference_contig = int(key[0]) + 1
        matching_contig = int(key[1]) + 1
        cosine_value = cosine_dictionary[i]

        if cosine_value > 0.95:
            matched_contigs[reference_contig].append(matching_contig)

    matched_contig_groups = list(matched_contigs.values())

    matched_contig_groups.sort()
    non_duplicate_contig_groups = list(
        matched_contig_groups for matched_contig_groups, _ in itertools.groupby(matched_contig_groups))

    print(non_duplicate_contig_groups)

    with open('grouped_contigs.pkl', 'wb') as f:
        pickle.dump(non_duplicate_contig_groups, f)


@click.command()
@click.option('-i', '--samfile', type=str, required=True, help="SAM file")
@click.option('-x', '--index', type=str, required=True, help="Index file generated with colour mapper")
def calculate_similarity(samfile, index):
    contig_reads = extract_mapping_data(samfile)
    contig_hex = convert_reads_to_hexidemical(contig_reads, index)
    calculate_cosine_similarity(contig_hex)
