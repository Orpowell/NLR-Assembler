import csv
import pickle
import logging
from itertools import combinations
import sys

import click
from sklearn.metrics.pairwise import cosine_similarity
from sklearn.feature_extraction.text import CountVectorizer

logging.basicConfig(stream=sys.stdout, format='%(asctime)s - %(message)s', datefmt='%d-%b-%y %H:%M:%S', level=logging.INFO)


def extract_mapping_data(sam_file):
    """
    The function takes a SAM file as an input and generates a dictionary of all the reads that have been mapped to a
    given refernce sequence in the SAM file as list.

    :param sam_file: reads mapped to contigs generated by assembly of reads
    :return: dictionary: {contig:[all mapped reads]}
    """
    contig_read_dictionary = {}
    logging.info("extracting infromation from SAM file...")
    with open(sam_file) as file:
        for line in file:

            # extract all reference sequence names from SAM file and add to dictionary with an empty list
            if line.startswith("@SQ"):
                seq_name = line.split("\t")[1][3:]
                contig_read_dictionary[seq_name] = []

            # skip the header added by samtools
            elif line.startswith("@PG"):
                continue

            # read line of SAM file and extract reference sequence and query sequence names then append query sequence
            # to a list in the dictionary using the reference sequence as a key
            else:
                read_info = line.split("\t")
                query_name = read_info[0]
                reference_name = read_info[2]
                # if a sequence is unmapped to a reference move onto the next read
                if reference_name == "*":
                    continue
                else:
                    contig_read_dictionary[reference_name].append(query_name)

    return contig_read_dictionary


def convert_reads_to_hexidemical(contig_read_dictionary, index):
    """
    The function takes a dictionary of contigs each with a list of reads mapped to that contig and converts the read
    names to a list of rgb values based on an index file. The list RGB values is then converted to a string of
    hexidecimal values that can later be analysed to compare the similarity of contigs.

    :param contig_read_dictionary: dictionary of reads mapped to each contig generated by extracting_mapping_data()
    :param index: csv file where each read has been assigned a colour based its adapter sequence
    :return: contig_hex_dicitonary: {contig: string of hexidecimal values assinged for each read}
    """
    logging.info("extracting infromation from index file...")
    with open(index, mode='r') as inp:
        index_reader = csv.reader(inp)
        ID_colour_dict = {rows[0]: rows[1] for rows in index_reader}

    logging.info("converting Seq IDs to RGB values...")
    contig_rgb = {contig: tuple(map(lambda x: ID_colour_dict[x], contig_read_dictionary[contig])) for contig in
                  contig_read_dictionary}

    logging.info("converting RGB values to hexidecimal text...")
    contig_hex_dictionary = {
        contig: " ".join(map(lambda x: '#%02x%02x%02x' % tuple(map(int, x.split(","))), contig_rgb[contig]))
        for contig in contig_rgb}

    return contig_hex_dictionary


def calculate_cosine_similarity(contig_hexidecimal_dictionary):
    """
    The function takes a dicitonary of contigs each with a string of hexidecimal values that represent all reads mapped
    to the contig. The cosine similarity of the hexidecimal string of all contigs are compared in a pairwise manner.
    Contigs with a cosine similarity above a given threshold (0.95) are grouped together into sets. Duplicate sets are
    removed and the non-duplicate list of group contig lists is saved as pickle.

    :param contig_hexidecimal_dictionary:
    :return: pickle file:
    """

    logging.info('calculating cosine similarity of contigs and grouping...')

    # generate list of all pairwise combinations of contigs
    combos = combinations(contig_hexidecimal_dictionary, 2)
    n_combos = len(tuple(combos))
    matched_contigs = {val: {val} for val in contig_hexidecimal_dictionary}

    i = 0
    for t1, t2 in combos:
        # extract text data for each contig and calculate cosine singularity
        text = (contig_hexidecimal_dictionary[t1], contig_hexidecimal_dictionary[t2])
        cosine = CountVectorizer()
        transform = cosine.fit_transform(text)
        cosine_matrix = cosine_similarity(transform)

        # If the cosine similairty of the two contigs is greater than 0.95 the
        if cosine_matrix[0, 1] > 0.95:
            matched_contigs[t1].add(t2)
            matched_contigs[t2].add(t1)

        # Track percentage of comparisons complete
        percentage_complete = i * 100 / n_combos
        if percentage_complete % 5 == 0:
            logging.info(f"{percentage_complete}% of contigs compared... ({i}/{n_combos})")
        i += 1

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
