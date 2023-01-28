import seaborn as sns
from sklearn.metrics.pairwise import cosine_similarity
from sklearn.feature_extraction.text import CountVectorizer

import csv
import click


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

    print('plotting cosine similarity...')
    cosine_plot = sns.heatmap(cosine_array.toarray(), cmap="crest", annot=True, xticklabels=labels, yticklabels=labels)

    print('saving plot...')
    fig = cosine_plot.get_figure()
    fig.savefig("cosine_similarity.png", bbox_inches='tight')


@click.command()
@click.option('-i', '--samfile', type=str, required=True, help="SAM file")
@click.option('-x', '--index', type=str, required=True, help="Index file generated with colour mapper")
def calculate_similarity(samfile, index):
    contig_reads = extract_mapping_data(samfile)
    contig_hex = convert_reads_to_hexidemical(contig_reads, index)
    calculate_cosine_similarity(contig_hex)
