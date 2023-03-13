import csv
import logging
import sys

from itertools import permutations
from collections import Counter

import click
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.feature_extraction.text import CountVectorizer
from sklearn.metrics.pairwise import cosine_similarity

logging.basicConfig(stream=sys.stdout, format='%(asctime)s - %(message)s', datefmt='%d-%b-%y %H:%M:%S',
                    level=logging.INFO)


def extract_mapping_data(sam_file, Blast_data):
    """
    Loads a SAM file and generates a dictionary of all the reads that have been mapped to a refernce sequence in the
    SAM file as a list of reads. The dictionary is then reduced to contain only contigs that have been annotated as
    NLRs using NLR annotator

    :param Blast_data:
    :param sam_file: reads mapped to contigs generated by assembly of reads
    :return: dictionary: {contig:[all mapped reads]}
    """
    contig_read_dictionary = {}
    logging.info("extracting reads from SAM file...")
    with open(sam_file) as file:
        for line in file:
            # extract all reference sequence names from SAM file and add to dictionary with an empty list
            if line.startswith("@SQ"):
                seq_name = line.split("\t")[1][3:]
                contig_read_dictionary[seq_name] = []

            # skip the header added by samtools
            elif line.startswith("@PG") or line.startswith("@HD"):
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

    logging.info("Loading BLAST data...")

    blast = pd.read_csv(Blast_data, header=None, sep="\t")

    blast_80 = blast[(blast[2] >= 80) & (blast[3] >= 80)]

    target_contigs = blast_80[1].unique()

    nlr_read_dictionary = {nlr: contig_read_dictionary[str(nlr)] for nlr in target_contigs}

    return nlr_read_dictionary


def convert_reads_to_hexidecimal(contig_read_dictionary, index):
    """
    The function takes a dictionary of contigs each with a list of reads mapped to that contig and converts the read
    names to a list of rgb values based on an index file. The list RGB values is then converted to a list of
    hexidecimal values that can later be analysed to compare the similarity of contigs.

    :param contig_read_dictionary: dictionary of reads mapped to each contig generated by extracting_mapping_data()
    :param index: csv file where each read has been assigned a colour based its adapter sequence
    :return: contig_hex_dicitonary: {contig: [hexidecimal values assinged for each read]}
    """
    logging.info("extracting information from index file...")
    with open(index, mode='r') as inp:
        index_reader = csv.reader(inp)
        ID_colour_dict = {rows[0]: rows[1] for rows in index_reader}

    logging.info("converting Seq IDs to RGB values...")
    contig_rgb = {contig: tuple(map(lambda x: ID_colour_dict[x], contig_read_dictionary[contig])) for contig in
                  contig_read_dictionary}

    logging.info("converting rgb values to hexidecimal...")
    contig_hex_dictionary = {
        contig: list(map(lambda x: '#%02x%02x%02x' % tuple(map(int, x.split(","))), contig_rgb[contig])) for
        contig in contig_rgb}

    return contig_hex_dictionary


def generate_cosine_matrix(contig_hex_dictionary):
    """
    Generates a text string for each contig using the hexidecimal profiles for each contig and then calculates the
    cosine similarity of all contigs in the dictionary. This is peformed in a pairwise manner and creates a matrix of all
    values.

    :param contig_hex_dictionary:
    :return: a 2D array of all cosine similarity values
    """
    logging.info("converting hexidecimal to strings...")
    contig_adapter_profiles = [" ".join(contig_hex_dictionary[contig]) for contig in
                               contig_hex_dictionary]

    logging.info("calculating cosine similarity...")
    count_array = CountVectorizer()
    profile_count_array = count_array.fit_transform(contig_adapter_profiles)
    cosine_array = cosine_similarity(profile_count_array)

    return cosine_array


def group_contigs(contig_hex, cosine_array, threshold):
    """
    generates a list of lists of all contigs that share a cosine similarity above the given threshold. The list is
    generated by grouping all contigs in a given row of the cosine matrix together. Each row in the matrix is assigned a
    number between 0 and n (max length of the matrix) as are the list of contigs from the dictionary used to generate the
    cosine matrix. Therefore, the two contigs compared at a given position in the matrix can be identified.

    EXAMPLE:

    contig dict:
    {0 : contig_1, 1 : contig_2, 2 : contig_3}

    cosine matrix:
    1.0 0.5 0.7
    0.5 1.0 0.8
    0.7 0.8 1.0

    @ threshold 0.7
    row 0: [0,2]
    row 1: [1,2]
    row 2: [0,1,2]

    convert positions to contig ids:
    row 0: [contig_1, contig_3]
    row 1: [contig_2, contig_3]
    row 2: [contig_1, contig_2, contig_3]

    The list of list is then filtered to remove a list that is a subset of another list. The count of each contig after
    this step in the list of lists should be 1. If this is not the case any list containing a contig that occurs more than
    once is removed from the list of list. These groups of contigs are often incorrect, massive and unhelpful. Any contig
    not found in the list of lists after this is added back to the list of lists as a list of itself. This ensures minimal
    data is lost during the filter process. Often complete NLRs are high promiscuous and have a high similairty to a range
    of partial and complete NLR contigs, thus the impact of this step is likely to be minimal. This list of lists is then
    returned as the output of the function

    :param contig_hex: the dictionary used to generate the cosine similarity matrix
    :param cosine_array: a matrix of the cosine similarity values generated by generate_cosine_matrix
    :param threshold: a float between 0 and 1.0
    :return: a list of Lists of contigs grouped together above the given threshold.
    """

    def get_key(key_value):
        for key, value in sudo_info.items():
            if key_value == value:
                return key

    def flatten(array):
        return [item for sublist in array for item in sublist]

    def get_docs(arr, docs_names, cosine_threshold):
        output_tuples = []
        for row in range(len(arr)):
            lst = [row + 1 + idx for idx, num in enumerate(arr[row, row + 1:]) if num >= cosine_threshold]
            for item in lst:
                output_tuples.append((docs_names[row], docs_names[item]))

        return output_tuples

    logging.info(f"-----{threshold}-----")
    logging.info(f"grouping contigs...")
    pairs = get_docs(cosine_array, list(contig_hex.keys()), threshold)
    grouping_dictionary = {key: [key] for key in contig_hex.keys()}

    [grouping_dictionary[pair[0]].append(pair[1]) for pair in pairs]

    groups = [v for k, v in grouping_dictionary.items() if len(v) > 1]

    logging.info(f"removing sub groups...")
    combos = permutations(groups, 2)
    sub_set = set()

    sudo_info = {k: v for k, v in enumerate(groups)}

    [sub_set.add(get_key(g1)) for g1, g2 in combos if set(g1).issubset(set(g2)) and g1 != g2]

    subs = [sudo_info[key] for key in sub_set]
    non_duplicate_group = [set(contig_group) for contig_group in groups if contig_group not in subs]

    logging.info(f"removing duplicate contigs...")
    bad_contigs = [k for k, v in Counter(flatten(non_duplicate_group)).items() if v > 1]

    for bc in bad_contigs:
        non_duplicate_group = [clean_group for clean_group in non_duplicate_group if bc not in clean_group]

    logging.info(f"re-adding contigs...")
    for val in list(contig_hex.keys()):
        if val not in flatten(non_duplicate_group):
            non_duplicate_group.append({val})

    return [list(contig) for contig in non_duplicate_group]


def find_optimal_grouping(cosine_threshold_dictionary):
    """
    Generates a bar chart that shows the number of lists with 2 or more elements in the list of lists of contigs at
    a range of cosine similarity thresholds. The key for the list with the greatest number of list with two or more
    elements is then returned. This allows for the cosine threshold with the greatest number of combined partial NLRs to
    be identified and selected for further analysis. In theory, more complete NLRs can be assembled from this list with
    the greatest number of elements greater than 2.

    :param cosine_threshold_dictionary: a dictionary containing the list of lists of grouped contigs at varying cosine
            thresholds
    :return: The key for list of lists in the dictionary with most number of lists with 2 or more elements
    """

    logging.info("Plotting of number of groups at cosine thresholds...")
    n_grouped_contigs = {str(k): len([contig for contig in v if len(contig) > 1]) for k, v in
                         cosine_threshold_dictionary.items()}
    x = n_grouped_contigs.keys()
    y = list(n_grouped_contigs.values())

    col = []
    optimal_threshold = max(y)
    for val in y:
        if val == optimal_threshold:
            col.append('red')
        else:
            col.append('blue')

    plt.bar(*zip(*n_grouped_contigs.items()), color=col)
    plt.xlabel('Cosine Similarity Threshold')
    plt.ylabel('# Grouped Contigs')

    for i in range(len(x)):
        plt.text(i, y[i] + 2.5, y[i], ha='center')

    plt.savefig('cosine_threshold.png', dpi=300)

    optimal_threshold = max(n_grouped_contigs, key=n_grouped_contigs.get)
    logging.info("barchart was saved as cosine_threshold.png")

    logging.info(f"Optimal cosine threshold identified... ({optimal_threshold})")

    return optimal_threshold


def write_grouped_contig_fasta(assembly, grouped_contigs):
    """
    converts lists of grouped contigs into a dictionary where each key is the contig ID and the value is concatenated
    sequence of all contigs in the group. For groups with more than 1 element, contigs spaced by 1000 N's. All sequence
    information comes from the orignal assembly.

    :param assembly: fasta file containing the sequences for all contigs in the original assembly
    :param grouped_contigs: list of grouped contigs generated by filter_by_strand
    :return: a fasta file containing the sequence all NLR contigs (new and old) generated by the pipeline.
    """

    print('writing contigs to fasta file...')
    contig_sequence = {}
    with open(assembly) as file:
        for line in file:
            if line.startswith('>'):
                contig_sequence[line[1:-1]] = next(file)[:-1]

    spacer = 'N' * 1000

    grouped_contigs = [list(map(str, contig)) for contig in grouped_contigs]

    seqs = {">" + "_".join(contig): f"{spacer}".join(list(map(contig_sequence.get, contig))) for contig in
            grouped_contigs}

    with open('grouped_assemblies.fa', 'w') as file:
        for k, v in seqs.items():
            file.write(k + "\n")
            file.write(v + "\n")


@click.command()
@click.option('-i', '--samfile', type=str, required=True, help="SAM file")
@click.option('-b', '--blast', type=str, required=True, help="blast file")
@click.option('-x', '--index', type=str, required=True, help="Index file generated with colour mapper")
@click.option('-a', '--assembly', type=str, required=True, help="assembly fasta")
def group(samfile, blast, index, assembly):
    nlr_contig_reads = extract_mapping_data(samfile, blast)
    contig_hex = convert_reads_to_hexidecimal(nlr_contig_reads, index)
    cosine_matrix = generate_cosine_matrix(contig_hex)

    threshold_dictionary = {threshold: group_contigs(contig_hex, cosine_matrix, threshold) for threshold in
                            [0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1]}

    optimal_threshold = find_optimal_grouping(threshold_dictionary)

    best_contig_grouping = threshold_dictionary[float(optimal_threshold)]

    write_grouped_contig_fasta(assembly, best_contig_grouping)
