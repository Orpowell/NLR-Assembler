import click
import logging
import numpy as np
import pandas as pd


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
                header = line.split()
                NLR_dict[header[0][1:]] = len(next(file))

    return {k: NLR(k, v) for k, v in NLR_dict.items()}


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
    return len(coverage_array), coverage_array.mean()


def determine_assembly_coverage(nlr, blast):
    nlr_data = load_NLR_data(nlr)
    blast_data = load_BLAST_data(blast)
    coverage_data, coverage_mean = calculate_NLR_coverage(nlr_data, blast_data)
    total_contigs = pd.read_csv(blast, header=None, sep="\t")[0].unique()
    logging.info(f"Average NLR Covereage: {coverage_mean}")
    logging.info(f"Total Contigs: {total_contigs}")
    return [coverage_mean, total_contigs]


@click.command()
@click.option('-b', '--draft_assemlby_blast', type=str, required=True, help="SAM file")
@click.option('-c', '--final_assemlby_blast', type=str, required=True, help="SAM file")
@click.option('-n', '--nlr_annotation_path', type=str, required=True, help="NLR annotator file")
def nlr_coverage(draft_assemlby_blast, final_assemlby_blast, nlr_annotation_path):
    logging.info("----- running NLR-Assembler nlr-coverage -----")
    draft_assembly_stats = determine_assembly_coverage(nlr_annotation_path, draft_assemlby_blast)
    final_assembly_stats = determine_assembly_coverage(nlr_annotation_path, final_assemlby_blast)

    cc = pd.DataFrame([draft_assembly_stats, final_assembly_stats], index=["draft", "final"], columns=["coverage of NLRs (%)", "contigs"]).transpose()
    cc["PD"] = [cc.final.iloc[0]-cc.draft.iloc[0], (cc[['draft', 'final']].pct_change(axis=1)['final'][1] * 100)]
    cc.to_csv("NLR_coverage.txt", sep="\t")
