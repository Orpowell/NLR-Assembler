import logging
import pickle
import sys

import click
import pandas as pd

logging.basicConfig(stream=sys.stdout, format='%(asctime)s - %(message)s', datefmt='%d-%b-%y %H:%M:%S',
                    level=logging.INFO)


def load_pickled_data(pickle_file):
    """
    Load pickled list of lists
    :param pickle_file: pickle file generate by group_by_cosine containing list of grouped contigs
    :return: unpickled list of grouped contigs
    """
    logging.info('loading contig grouping data...')
    with open(pickle_file, 'rb') as f:
        contig_grouping = pickle.load(f)

    return contig_grouping


def load_annotation_data(NLR_annotation):
    """
    Load NLR annotations as two dictionary with the strand and annotation for each contig
    E.g {contig_1 : CC-NBARC-LRR} and {contig_1 : +}
    :param NLR_annotation: Path to file for NLR annotator output for assembly
    :return: Dicitonaries for the strand and annotation of each contig individually
    """
    logging.info('loading NLR annotations')
    contig_strand_dictionary = {}
    contig_annotation_dictionary = {}
    with open(NLR_annotation) as annotations:
        for line in annotations:
            line_data = line.split("\t")
            contig_strand_dictionary[line_data[0]] = line_data[5]
            contig_annotation_dictionary[line_data[0]] = line_data[2]

    return contig_annotation_dictionary, contig_strand_dictionary


def load_mapping_data(contig_list, minimap2_alignment):

    logging.info("loading contig mapping data...")
    keys = [item for sublist in contig_list for item in sublist]
    with open(minimap2_alignment) as file:

        mapping_data = {}
        for line in file:
            if line.startswith('@'):
                continue

            else:
                row = line.split()
                mapping_data[row[0]] = [row[2], row[3], len(row[9])]

        nlr_data = {k: mapping_data[k] for k in keys}

    return nlr_data


def filter_by_strand(strand_info, contig_groups):
    """
    Sorts a list of grouped contigs by the strand each contig is on using information on the contigs
    generated by NLR annotator. This takes place in 4 steps:

        1. The list is split into list with a length of 1, A, and a length greater than 1, B
        2. B is divided into list where all contigs have the same strand, C and mixed strands, D
        3. Each list in D is split into positive and negative strand lists and appended to a new list, E
        4. Lists A,C and E are combined to generate a list of contigs grouped by strand

    :param strand_info: Dictionary with the strand of each contig generated by load_annotation_data
    :param contig_groups: list of lists of grouped contigs loaded from pickle with load_pickled_data
    :return: list of grouped contigs that have been sorted into subgroups (new contig lists) by strand
    """
    logging.info('sorting contigs by strand...')
    singletons = []
    grouped = []  #
    same_stranded_groups = []
    mix_stranded_groups = []  #
    seperated_stranded_groups = []

    for val in contig_groups:
        if len(val) == 1:
            singletons.append(val)
        else:
            grouped.append(val)

    for val in grouped:
        strand_group = [strand_info[n] for n in val]

        if len(set(strand_group)) == 1:
            same_stranded_groups.append(val)

        else:
            mix_stranded_groups.append(val)

    for val in mix_stranded_groups:
        strand_group = [strand_info[n] for n in val]
        temp_pos = []
        temp_neg = []

        for n, strand in enumerate(strand_group):
            if strand == "+":
                temp_pos.append(val[n])
            else:
                temp_neg.append(val[n])

        seperated_stranded_groups.append(temp_pos)
        seperated_stranded_groups.append(temp_neg)

    return seperated_stranded_groups + singletons + same_stranded_groups


def annotate_grouped_contigs(sorted_contigs, annotation_dictionary, error_mapping):
    """
    Creates new NLR annotations for groups of contigs and counts the occurence of each NLR annotation in the grouped contigs
    E.g ['CC-NBARC', 'NBARC', 'NBARC-LRR'] would be classified as CC-NBARC-LRR

    :param error_mapping:
    :param sorted_contigs: list of contig groups organised by strand using filter_by_strand
    :param annotation_dictionary: A dictionary with the NLR annotation of each ech contig generated by load_annotation_data
    :return: A dictionary with a total count of each NLR type found in the list of grouped contigs
    """

    def congregate_annotations(contig_group):
        contig_group.sort()

        if len(contig_group) == 1:
            return contig_group[0]

        elif contig_group == ['NBARC', 'NBARC-LRR']:
            return 'NBARC-LRR'

        elif contig_group == ['NBARC', 'NBARC', 'NBARC-LRR']:
            return 'NBARC-LRR'

        elif contig_group == ['CC-NBARC', 'NBARC']:
            return 'CC-NBARC'

        elif contig_group == ['CC-NBARC', 'NBARC', 'NBARC']:
            return 'CC-NBARC'

        elif contig_group == ['NBARC', 'NBARC']:
            return 'NBARC'

        elif contig_group == ['CC-NBARC', 'NBARC', 'NBARC-LRR']:
            return 'CC-NBARC-LRR'

        elif contig_group == ['CC-NBARC', 'NBARC-LRR']:
            return 'CC-NBARC-LRR'

        else:
            return 'CLUSTER'

    def process_mapping_data(mapping_dict, contig):
        raw_data = list(map(mapping_dict.get, contig))
        chromosome = set([val[0] for val in raw_data])

        if chromosome == {'ChrUnknown'}:
            return 0

        elif len(chromosome) > 1:
            return 1

        else:
            positions = [int(val[1]) for val in raw_data]
            start = min(positions)
            length = sum([int(val[2]) for val in raw_data]) + ((len(raw_data) - 1) * 1000)

            position = f'{" ".join(chromosome)}:{start}-{start + length}'

            return position

    if error_mapping is not None:
        logging.info('annotating contigs (includes contig mapping data)...')
        mapping = {"_".join(contig): [" ".join([annotation_dictionary[n] for n in contig]),
                                      congregate_annotations([annotation_dictionary[n] for n in contig]),
                                      process_mapping_data(error_mapping, contig)] for contig in sorted_contigs}
        data = pd.DataFrame(mapping).transpose()
        data.columns = ['grouped_annotations', 'overall_annotation', 'position']
        data.to_csv('assembly_annotations.txt', sep='\t')

        total_count = len(data)
        mismatch_count = len(data[data.position == 1])
        unknown_count = len(data[data.chromosome == 0])

        logging.info(f"Unmappable contigs: {100 * unknown_count / total_count:.3}% ({unknown_count} of {total_count})")
        logging.info(f"Mismatched contigs: {100 * mismatch_count / total_count:.3}% ({mismatch_count} of {total_count})")
        logging.info(f"Total contig errors: {100 * (unknown_count + mismatch_count) / total_count:.3}% {(unknown_count + mismatch_count)} of {total_count})")

    else:
        logging.info('annotating contigs (without contig mapping data)...')
        annotations = {"_".join(contig): [" ".join([annotation_dictionary[n] for n in contig]),
                                          congregate_annotations([annotation_dictionary[n] for n in contig])] for contig
                       in
                       sorted_contigs}

        data = pd.DataFrame(annotations).transpose()
        data.columns = ['grouped_annotations', 'overall_annotation']
        data.to_csv('assembly_annotations.txt', sep='\t')

    return data


def generate_assembly_grouping_statistics(assembly_nlr_data, grouping_nlr_data):
    """
    laods the NLR annotations for the grouped contigs and the original assembly and compares the difference in number of
    each annotation. Changes are measures as the difference in count and percentage change respective to the original
    assembly.

    :param assembly_nlr_data: path to NLR annotations for the original assembly
    :param grouping_nlr_data: dicitonary of NLR annotation counts for the newly grouped contigs
    :return: a csv file containing the comparison of the original assembly and pipelined assembly annotations
    """
    logging.info('calculating assembly statistics')
    assembly = pd.read_csv(assembly_nlr_data, sep="\t", header=None)
    comparison_counts = {'assembly': assembly[2].value_counts().to_dict(),
                         'cosine_grouping': grouping_nlr_data.overall_annotation.value_counts().to_dict()}
    data = pd.DataFrame(comparison_counts)
    data.drop([index for index in data.index if 'TIR' in index], inplace=True)

    data['relative_change'] = data.cosine_grouping - data.assembly

    data['percentage_change'] = 100 * (data.cosine_grouping / data.assembly) - 100
    data.fillna(0).astype(int)

    data.to_csv('assembly_statistics.txt', sep='\t')


def write_grouped_contig_fasta(assembly, grouped_contigs):
    """
    converts lists of grouped contigs into a dictionary where each key is the contig ID and the value is concatenated
    sequence of all contigs in the group. For groups with more than 1 element, contigs spaced by 1000 N's. All sequence
    information comes from the orignal assembly.

    :param assembly: fasta file containing the sequences for all contigs in the original assembly
    :param grouped_contigs: list of grouped contigs generated by filter_by_strand
    :return: a fasta file containing the sequence all NLR contigs (new and old) generated by the pipeline.
    """

    logging.info('writing contigs to fasta file...')
    contig_sequence = {}
    with open(assembly) as file:
        for line in file:
            if line.startswith('>'):
                contig_sequence[line[1:-1]] = next(file)[:-1]

    spacer = 'N' * 1000
    seqs = {">" + "_".join(contig): f"{spacer}".join(list(map(contig_sequence.get, contig))) for contig in
            grouped_contigs}

    with open('grouped_assemblies.fa', 'w') as file:
        for k, v in seqs.items():
            file.write(k + "\n")
            file.write(v + "\n")


@click.command()
@click.option('-i', '--input_data', type=str, required=True, help="pickle file generated with calculate_cosine")
@click.option('-n', '--nlr', type=str, required=True, help="NLR annotator file")
@click.option('-c', '--assembly_contigs', type=str, required=True, help="assembly contigs fasta file")
@click.option('-e', '--error_mapping_data', type=str, required=False, default=None, help="assembly contigs fasta file")
def annotate(input_data, nlr, assembly_contigs, error_mapping_data):
    contig_annotations, contig_strands = load_annotation_data(nlr)
    contigs = load_pickled_data(input_data)
    strand_filtered_contigs = filter_by_strand(contig_strands, contigs)

    if error_mapping_data is not None:
        mapping_data = load_mapping_data(strand_filtered_contigs, error_mapping_data)

    else:
        mapping_data = None

    write_grouped_contig_fasta(assembly_contigs, strand_filtered_contigs)
    annotated_grouped_contigs = annotate_grouped_contigs(strand_filtered_contigs, contig_annotations,
                                                         mapping_data)
    generate_assembly_grouping_statistics(nlr, annotated_grouped_contigs)
