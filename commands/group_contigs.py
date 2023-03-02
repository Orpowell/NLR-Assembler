import logging
import pickle
import click
import pandas as pd


def load_pickled_data(pickle_file):
    logging.info('loading contig grouping data...')
    with open(pickle_file) as f:
        contig_grouping = pickle.load(f)

    return contig_grouping


def load_annotation_data(NLR_annotation):
    logging.info('loading NLR annotations')
    contig_strand_dictionary = {}
    contig_annotation_dictionary = {}
    with open(NLR_annotation) as annotations:
        for line in annotations:
            line_data = line.split("\t")
            contig_strand_dictionary[line_data[0]] = line_data[5]
            contig_annotation_dictionary[line_data[0]] = line_data[2]

    return contig_annotation_dictionary, contig_strand_dictionary


def filter_by_strand(strand_info, contig_groups):
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


def annotate_grouped_contigs(sorted_contigs, annotation_dictionary):
    logging.info('annotating contigs...')
    annotations = sorted([sorted([annotation_dictionary[n] for n in val]) for val in sorted_contigs], key=len,
                         reverse=True)

    nlr_count = {k: 0 for k in set([item for sublist in annotations for item in sublist])}
    nlr_count['CLUSTER'] = 0

    for val in annotations:
        if len(val) == 1:
            nlr_count[val[0]] += 1

        elif val == ['NBARC', 'NBARC-LRR']:
            nlr_count['NBARC-LRR'] += 1

        elif val == ['CC-NBARC', 'NBARC']:
            nlr_count['CC-NBARC'] += 1

        elif val == ['NBARC', 'NBARC']:
            nlr_count['NBARC'] += 1

        elif val == ['CC-NBARC', 'NBARC', 'NBARC-LRR']:
            nlr_count['CC-NBARC-LRR'] += 1

        elif val == ['CC-NBARC', 'NBARC-LRR']:
            nlr_count['CC-NBARC-LRR'] += 1

        else:
            nlr_count['CLUSTER'] += 1

            for annot in val:
                nlr_count[annot] += 1

    return nlr_count


def generate_assembly_grouping_statistics(assembly_nlr_data, grouping_nlr_data):
    logging.info('calculating assembly statistics')
    assembly = pd.read_csv(assembly_nlr_data, sep="\t", header=None)
    comparison_counts = {'assembly': assembly[2].value_counts().to_dict(), 'cosine_grouping': grouping_nlr_data}
    data = pd.DataFrame(comparison_counts)
    data.drop([index for index in data.index if 'TIR' in index], inplace=True)

    data['relative_change'] = data.cosine_grouping - data.assembly

    data['percentage_change'] = 100 * (data.cosine_grouping / data.assembly) - 100
    data.fillna(0).astype(int)

    data.to_csv('assembly_statistics.txt', sep='\t')


def write_grouped_contig_fasta(assembly, grouped_contigs):
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
def group_contigs(input_data, nlr, assembly_contigs):
    contig_annotations, contig_strands = load_annotation_data(nlr)
    contigs = load_pickled_data(input_data)
    strand_filtered_contigs = filter_by_strand(contig_strands, contigs)
    annotated_grouped_contigs = annotate_grouped_contigs(strand_filtered_contigs, contig_annotations)
    generate_assembly_grouping_statistics(nlr, annotated_grouped_contigs)
    write_grouped_contig_fasta(assembly_contigs, annotated_grouped_contigs)
