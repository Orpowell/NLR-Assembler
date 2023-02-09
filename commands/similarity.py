import csv
import pickle
import logging
from itertools import combinations

import click
from sklearn.metrics.pairwise import cosine_similarity
from sklearn.feature_extraction.text import CountVectorizer

logging.basicConfig(format='%(asctime)s - %(message)s', datefmt='%d-%b-%y %H:%M:%S', level=logging.INFO)


def extract_mapping_data(sam_file):
    contig_read_dictionary = {}
    logging.info("extracting infromation from SAM file...")
    with open(sam_file) as file:
        for line in file:
            if line.startswith("@SQ"):
                seq_name = line.split("\t")[1][3:]
                contig_read_dictionary[seq_name] = []

            elif line.startswith("@PG"):
                continue

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
    logging.info("extracting infromation from index file...")
    with open(index, mode='r') as inp:
        index_reader = csv.reader(inp)
        ID_colour_dict = {rows[0]: rows[1] for rows in index_reader}

    logging.info("converting Seq IDs to RGB values...")
    contig_rgb = {contig: list(map(lambda x: ID_colour_dict[x], contig_read_dictionary[contig])) for contig in
                  contig_read_dictionary}

    logging.info("converting RGB values to hexidecimal text values...")
    contig_hex_dictionary = {
        contig: " ".join(list(map(lambda x: '#%02x%02x%02x' % tuple(map(int, x.split(","))), contig_rgb[contig]))) for
        contig in contig_rgb}

    return contig_hex_dictionary


def calculate_cosine_similarity(contig_hexidecimal_dictionary):

    combos = combinations(contig_hexidecimal_dictionary, 2)

    cosine_dictionary = {}

    logging.info('calculating cosine similarity of contigs...')
    for t1, t2 in combos:
        text_list = [contig_hexidecimal_dictionary[t1], contig_hexidecimal_dictionary[t2]]
        cosine = CountVectorizer()
        transform = cosine.fit_transform(text_list)
        co_sim = cosine_similarity(transform)

        n = f"{t1} {t2}"
        m = f"{t2} {t1}"

        cosine_dictionary[n] = round(co_sim[0, 1], 2)
        cosine_dictionary[m] = round(co_sim[0, 1], 2)

    matched_contigs = {val: {val} for val in contig_hexidecimal_dictionary}

    logging.info('grouping contigs...')
    for i in cosine_dictionary:
        key = str(i).split()
        reference_contig = key[0]
        matching_contig = key[1]

        if cosine_dictionary[i] > 0.95:
            matched_contigs[reference_contig].add(matching_contig)

    logging.info('removing duplicates contig groups...')
    matched_contig_groups = list(matched_contigs.values())
    non_duplicate_contig_groups = [list(x) for x in set(tuple(x) for x in matched_contig_groups)]

    logging.info('writing data to pickle...')
    with open('grouped_contigs.pkl', 'wb') as f:
        pickle.dump(non_duplicate_contig_groups, f)


@click.command()
@click.option('-i', '--samfile', type=str, required=True, help="SAM file")
@click.option('-x', '--index', type=str, required=True, help="Index file generated with colour mapper")
def calculate_similarity(samfile, index):
    contig_reads = extract_mapping_data(samfile)
    contig_hex = convert_reads_to_hexidemical(contig_reads, index)
    calculate_cosine_similarity(contig_hex)
