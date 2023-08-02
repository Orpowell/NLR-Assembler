import logging
import pandas as pd
import click


def load_blast_data(blast_data):
    """
    Load blast data into a dataframe, removes all alignments with unknown chromosomes and calculates the percentage 
    query coverage.
    
    :param blast_data: raw blastn data
    :return blast: a DataFrame of the cleaned blast data 
    
    """
    
    logging.info("loading BLAST data...")
    blast = pd.read_csv(blast_data, sep="\t", header=None) # load data
    blast = blast[blast[1] != "ChrUnknown"] # Remove unknown chromosomes
    blast[13] = (blast[7] - blast[6]) * 100 / blast[12] # calculate percentage coverage

    return blast


def load_grouped_contigs(contig_data):
    """
    Load the final assembly and filter for the contig groups
    
    :param contig_data: final assembly produced by NLR-Assembler
    :return grouped_contigs: a list of all contig groups in the final assembly
    """
    
    logging.info("loading assembly data...")
    with open(contig_data) as file:
        grouped_contigs = [line[1:-1] for line in file if line.startswith(">") and "_" in line]

    return grouped_contigs


class contig:
    '''
    contig holds all the blast data for a set of contigs grouped together in the final assembly. The class is initialised 
    using the name of the grouped contig and raw blast data. The raw blast data is then filtered down to contigs only in 
    the contig group.
    
    Attributes
    ----------
    name: name of the contig group in the final assembly
    sub_contigs: list of all constituent contigs in the contig group
    blast: blast data filtered for all contigs in sub_contigs
    chromosome: most likely cnadidate chromosome for the contigs
    genome_coverage: genomic range covered by all contigs on the chromosome
    contig_count: number of contigs in sub_contigs with a blast hit on chromosome
    contig_total: total number of constituent contigs
    genome_coordinates: genomic coordinates for the genome coverage for use by IGV
    
    Methods
    -------
    calculate_coverage()
        calculates the region covered by the grouped contigs
    '''
    
    def __init__(self, name, data):
        logging.info(f"Initializing {self.name}")
        self.name = name
        self.sub_contigs = list(map(int, name.split("_")))
        self.blast = data[data[0].isin(self.sub_contigs)]
        self.chromosome = None
        self.genome_coverage = None
        self.contig_count = 0
        self.contig_total = len(self.sub_contigs)
        self.genome_coordinates = 'XXX'
        

    def calculate_coverage(self):
        logging.info(f"Calculating coverage {self.name}")
        
        # Reorganise blast hit coordinates to ensure all hit are in the same direction
        cond = self.blast[8] > self.blast[9]
        self.blast.loc[cond, [8, 9]] = self.blast.loc[cond, [9, 8]].values
        
        # Identify the highest likelihood chromosome for the contig group and filter the blast data for the chromosome
        best_chromosome = self.blast.sort_values(by=[11, 13], ascending=False).drop_duplicates(0)[
            1].value_counts().idxmax()
        filtered = self.blast[self.blast[1] == best_chromosome]
        final = filtered.sort_values(by=[11, 13], ascending=False).drop_duplicates(0)

        # Identify min/max of genomic range on target chromosome
        start = int(min(final[8]))
        end = int(max(final[9]))

        #update class atrributes
        self.chromosome = best_chromosome
        self.genome_coverage = end - start
        self.contig_count = len(final)
        self.genome_coordinates = f"{best_chromosome}:{start}-{end}"

    
    def get_data(self):
        """
        return all necessary data from class after coverage has been calculated
        
        return: all class attributes
        """
        return self.name, self.contig_count, self.contig_total, self.chromosome, self.genome_coverage, self.genome_coordinates


@click.command()
@click.option('-b', '--blast', type=str, required=True, help="BLAST file")
@click.option('-a', '--assembly', type=str, required=True, help="Assembly")
def contig_coverage(blast, assembly):
    '''
    Calculates the region of the genome covered by contigs grouped together in the final assembly
    
    :param assembly: path to the final assembly
    :param blast: path to the raw blast data
    :return None: A csv file is saved with the output file
    '''
    logging.info("----- running NLR-Assembler query-coverage -----")
    blast_data = load_blast_data(blast)
    group_data = load_grouped_contigs(assembly)

    # generate a dictionary of all contig groups 
    contig_dict = {contig_group: contig(contig_group, blast_data) for contig_group in group_data}

    # Generate the final output csv containing data for all contig groups
    for k in contig_dict.keys():
        try:
            contig_dict[k].calculate_coverage()
        except ValueError:
            logging.error(f"Could not calculate coverage for {k}")

    data_matrix = [contig_dict[k].get_data() for k in contig_dict.keys()]
    summary = pd.DataFrame(data_matrix)
    
    # Generate and log summary statistics
    logging.info(
        f"Percentage contigs covering 60 Kb or less: {len(summary[summary[5] < 60000][5]) / len(summary) * 100:.4}% ({len(summary[summary[5] < 60000][5])} of {len(summary)})")
    logging.info(
        f"Percentage contigs covering 1 Mb or less: {len(summary[summary[5] < 100000][5]) / len(summary) * 100:.4}% ({len(summary[summary[5] < 100000][5])} of {len(summary)})")

    logging.info("Saving query coverage data to query_coverage.txt ...")
    summary.to_csv("query_coverage.txt", sep="\t")