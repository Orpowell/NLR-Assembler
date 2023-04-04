import csv
import logging
import sys
import numpy as np
from collections import Counter

import click
import pandas as pd
from sklearn.feature_extraction.text import TfidfVectorizer
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

    logging.info("Removing poor quality contigs...")

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
    :return: a 2D nlr_dict of all cosine similarity values
    """
    logging.info("converting hexidecimal to strings...")
    contig_adapter_profiles = [" ".join(contig_hex_dictionary[contig]) for contig in
                               contig_hex_dictionary]

    logging.info("calculating cosine similarity...")
    tfidf_array = TfidfVectorizer()
    profile_count_array = tfidf_array.fit_transform(contig_adapter_profiles)
    cosine_array = cosine_similarity(profile_count_array)
    tf = pd.DataFrame(cosine_array, index=contig_hex_dictionary.keys(), columns=contig_hex_dictionary.keys())
    return tf


def group_contigs(cosine_dataframe):
    logging.info("Grouping contigs...")

    def dynamic_grouping(df_row, dataframe):
        test = dataframe[df_row].nlargest(n=6)
        threshold = test.iloc[1]

        if threshold < 0.5:
            contig_group = [test.index[0]]
            return sorted(contig_group)

        else:
            contig_group = list(test[test > (test.iloc[1] - 0.1)].index)
            return sorted(contig_group)

    contig_groups = []
    i = 0
    points = np.percentile(np.arange(len(cosine_dataframe)), [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100]).astype(int)
    for row in cosine_dataframe.columns.unique():

        if i in points:
            logging.info(f"{round((i * 100) / len(cosine_dataframe))}% complete")

        target_contig_group = dynamic_grouping(row, cosine_dataframe)
        target_contig_group.sort()
        i += 1
        if target_contig_group not in contig_groups:
            contig_groups.append(target_contig_group)

    return contig_groups


def merge_contig_groups(contig_groups):
    def merge_overlapping_lists(lists):
        merged_lists = {}
        i = 0
        points = np.percentile(np.arange(len(lists)), [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100]).astype(int)
        for lst in lists:
            if i in points:
                logging.info(f"{round((i * 100) / len(lists))}% complete")
            overlap_found = False
            for key in merged_lists:
                if any(elem in lst for elem in merged_lists[key]):
                    merged_lists[key].extend(lst)
                    overlap_found = True
                    break
            if not overlap_found:
                merged_lists[tuple(lst)] = lst
            i += 1
        return list(map(set, list(merged_lists.values())))

    logging.info("Merging overlapping contigs...")
    merged_contig_groups = merge_overlapping_lists(contig_groups)
    merged_list = list(map(sorted, map(list, merged_contig_groups)))
    merged_list.sort(key=len, reverse=True)

    return merged_list


def remove_subset_contigs(merged_contig_groups):
    def flatten(array):
        return [item for sublist in array for item in sublist]

    logging.info("removing subset contigs...")
    nr_list = list(map(sorted, map(list, merged_contig_groups)))
    nr_list.sort(key=len, reverse=True)
    count = Counter(flatten(nr_list))
    dupes = [k for k in count if count[k] > 1]

    i = 0
    points = np.percentile(np.arange(len(dupes)), [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100]).astype(int)
    for contig in dupes:
        if i in points:
            logging.info(f"{round((i * 100) / len(dupes))}% complete")

        i += 1
        query_group = [cgroup for cgroup in nr_list if contig in cgroup]

        if len(query_group) < 2:
            continue
        remove_group = min(query_group, key=len)
        nr_list.remove(remove_group)

    return nr_list


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
    logging.info("----- running NLR-Assembler group -----")
    nlr_contig_reads = extract_mapping_data(samfile, blast)
    contig_hex = convert_reads_to_hexidecimal(nlr_contig_reads, index)
    cosine_matrix = generate_cosine_matrix(contig_hex)
    raw_contig_grouping = group_contigs(cosine_matrix)
    merged_contig_grouping = merge_contig_groups(raw_contig_grouping)
    final_contig_grouping = remove_subset_contigs(merged_contig_grouping)
    write_grouped_contig_fasta(assembly, final_contig_grouping)
