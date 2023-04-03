import logging

import pandas as pd
import click


def load_blast_data(blast_data):
    logging.info("loading BLAST data...")
    blast = pd.read_csv(blast_data, sep="\t", header=None)
    blast = blast[blast[1] != "ChrUnknown"]
    blast[13] = (blast[7] - blast[6]) * 100 / blast[12]

    return blast


def load_grouped_contigs(contig_data):
    logging.info("loading assembly data...")
    with open(contig_data) as file:
        grouped_contigs = [line[1:-1] for line in file if line.startswith(">") and "_" in line]

    return grouped_contigs


class contig:
    def __init__(self, name, data):
        self.name = name
        self.sub_contigs = list(map(int, name.split("_")))
        self.blast = data[data[0].isin(self.sub_contigs)]
        self.chromosome = None
        self.genome_coverage = None
        self.contig_count = 0
        self.contig_total = len(self.sub_contigs)
        self.genome_coordinates = 'XXX'
        logging.info(f"Initializing {self.name}")

    def calculate_coverage(self):
        logging.info(f"Calculating coverage {self.name}")
        cond = self.blast[8] > self.blast[9]
        self.blast.loc[cond, [8, 9]] = self.blast.loc[cond, [9, 8]].values
        best_chromosome = self.blast.sort_values(by=[11, 13], ascending=False).drop_duplicates(0)[
            1].value_counts().idxmax()
        filtered = self.blast[self.blast[1] == best_chromosome]
        final = filtered.sort_values(by=[11, 13], ascending=False).drop_duplicates(0)

        start = int(min(final[8]))
        end = int(max(final[9]))

        self.chromosome = best_chromosome
        self.genome_coverage = end - start
        self.contig_count = len(final)
        self.genome_coordinates = f"{best_chromosome}:{start}-{end}"

    def get_data(self):
        return self.name, self.contig_count, self.contig_total, self.chromosome, self.genome_coverage, self.genome_coordinates


@click.command()
@click.option('-b', '--blast_path', type=str, required=True, help="BLAST file")
@click.option('-c', '--contig_path', type=str, required=True, help="Contig Assembly")
def query_coverage(blast_path, contig_path):
    logging.info("----- running NLR-Assembler query-coverage -----")
    blast_data = load_blast_data(blast_path)
    group_data = load_grouped_contigs(contig_path)

    contig_dict = {contig_group: contig(contig_group, blast_data) for contig_group in group_data}
    [contig_dict[k].calculate_coverage() for k in contig_dict.keys()]

    for k in contig_dict.keys():
        try:
            contig_dict[k].calculate_coverage()
        except ValueError:
            logging.error(f"Could not calculate coverage for {k}")

    data_matrix = [contig_dict[k].get_data() for k in contig_dict.keys()]
    summary = pd.DataFrame(data_matrix)
    logging.info("Saving query coverage data to query_coverage.txt ...")
    summary.to_csv("query_coverage.txt", header=None, sep="\t")
    logging.info(
        f"Percentage contigs covering 60 Kb or less: {len(summary[summary[5] < 60000][5]) / len(summary) * 100:.4}% ({len(summary[summary[5] < 60000][5])} of {len(summary)})")
    logging.info(
        f"Percentage contigs covering 1 Mb or less: {len(summary[summary[5] < 100000][5]) / len(summary) * 100:.4}% ({len(summary[summary[5] < 100000][5])} of {len(summary)})")

