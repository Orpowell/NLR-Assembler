import click
import logging
import numpy as np
import matplotlib.pyplot as plt


class NLR:
    def __init__(self, name, seq_length):
        self.name = name
        self.coverage_array = np.zeros(seq_length)
        self.coverage = 0

    def add_coverage(self, start, end):
        for n in list(range(start - 1, end, 1)):
            self.coverage_array[n] += 1

    def get_coverage(self):
        return (np.count_nonzero(self.coverage_array) * 100) / len(self.coverage_array)


def load_NLR_data(fasta_file):
    NLR_dict = {}
    logging.info("loading NLR sequences...")
    with open(fasta_file) as file:
        for line in file:
            if line.startswith(">"):
                NLR_dict[line[1:]] = len(next(file))

    return NLR_dict


def load_BLAST_data(blast_table):
    BLAST_list = []
    logging.info("loading BLAST data...")
    with open(blast_table) as file:
        for n, line in enumerate(file):
            row = line.split()
            positions = [int(row[8]), int(row[9])]
            BLAST_list.append([row[1], min(positions), max(positions)])

    return BLAST_list


def calculate_NLR_coverage(nlr_class_dict, blast_array):
    logging.info('calculating NLR coverage...')
    [nlr_class_dict[val[0]].add_coverage(val[1], val[2]) for val in blast_array]
    coverage_array = np.asarray([nlr.get_coverage() for nlr in nlr_class_dict.values()])
    return coverage_array, coverage_array.mean()


def plot_coverage_histogram(coverage_array):
    logging.info("plotting histogram of coverage...")
    ax = plt.figure().gca()
    ax.yaxis.set_major_locator(plt.MaxNLocator(integer=True))
    plt.hist(coverage_array)
    plt.xlabel("Percentage Coverage %")
    plt.ylabel("Frequency")
    plt.xlim(0, 101)
    plt.savefig("coverage_histogram.png", dpi=300)
    plt.show()


@click.command()
@click.option('-b', '--blast', type=str, required=True, help="SAM file")
@click.option('-n', '--nlr', type=str, required=True, help="NLR annotator file")
def coverage(nlr, blast):
    nlr_data = load_NLR_data(nlr)
    nlr_dict = {k: NLR(k, v) for k, v in nlr_data.items()}
    blast_data = load_BLAST_data(blast)
    coverage_data, coverage_mean = calculate_NLR_coverage(nlr_dict, blast_data)
    plot_coverage_histogram(coverage_data)
    logging.info(f"Average NLR Covereage: {coverage_mean()}")
