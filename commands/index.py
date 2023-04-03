import os
import shutil
import tempfile
import re
import random
import csv
import multiprocessing

from functools import partial
from itertools import zip_longest
from multiprocessing import get_context

import click

random.seed(22)  # Ensures adapter sequence : colour dictionary generated is the same each time


# splits a fastq file into a given batch size
def split_fastq(fastq, batch_size, temporary_directory):
    """
    :param fastq: fastq file
    :param batch_size: number of sequences per subfile
    :param temporary_directory: location of the subdirectory
    :return: a list of all subfiles created, use to generate pool for multiprocessing
    """

    def grouper(n, iterable, fillvalue=None):
        """Collect data into fixed-length chunks or blocks"""
        args = [iter(iterable)] * n
        return zip_longest(fillvalue=fillvalue, *args)

    file_list = []

    with open(fastq) as file:
        for i, g in enumerate(grouper(batch_size * 4, file, fillvalue=None)):
            with tempfile.NamedTemporaryFile('w', delete=False) as fout:
                for j, line in enumerate(g, 1):  # count number of lines in group
                    if line is None:
                        j -= 1  # don't count this line
                        break
                    fout.write(line)
            file_name = f'{temporary_directory}/fastq_{i}.fastq'
            os.rename(fout.name, file_name)
            file_list.append(file_name)

    return file_list


def get_adapter_info(adapter_file):
    """
    :param adapter_file: text file with list of adapters (1 per line)
    :return adapter_sequences: a list of all adapters, used to correct sequencing errors
    :return seq_to_color_dict: dictionary with each adapter assigned to a random colour, used to assign colour to correct reads
    """

    print(f'Loading adapter sequence data from {adapter_file}')
    with open(adapter_file) as file:
        raw_adapters = file.readlines()
        adapter_sequences = [seq[:-1] for seq in raw_adapters]  # remove newline character from each sequence
        seq_to_color_dict = {adapter: f'{random.randint(0, 255)},{random.randint(0, 255)},{random.randint(0, 255)}' for
                             adapter in adapter_sequences}

        return adapter_sequences, seq_to_color_dict


def index_subfile(adapter_array, adapter_dict, fastq):
    """
    Assign each sequence ID in a fastq file a colour based on a dictionary provided
    :param adapter_array: a list of all possible adapters available
    :param adapter_dict: a dictionary with each adapter assigned a random colour
    :param fastq: a fastq file containing linked-reads with the same adapters used for other inputs
    :return: writes a CSV file with each sequence ID assigned a colour based on adapter sequence
    """

    def adapter_correcter(raw_adapter, search_set):
        """
        Corrects a given adapter if it contains a sequencing error
        :param raw_adapter: 16 character string (adapter) extracted from a raw fastq read
        :param search_set: a list of all availavle adapters
        :return: If adapter cannot be corrected None  otherwise the corrected adapter is returned
        """
        # If there's an ambigous base (N) in the adapter sequence search for the true sequence
        if 'N' in raw_adapter:
            search_string = raw_adapter.replace('N', '[ATCG]')  # convert string to regex
            true_adapter = re.search(fr'{search_string}', str(search_set))  # search for real adapter

            # If no adapter is found return true
            if true_adapter:
                return true_adapter.group()

        else:
            return raw_adapter

    # Opens a fastq file and creates a dictionary of sequence id: adapter (First 16 bases of each read)
    print(f'extracting sequences from {fastq}...')
    seq_dict = {}
    with open(fastq) as file:
        for line in file:
            if re.match('@', line) is not None:
                seq_id = line.split()[0][1:]
                seq_dict[seq_id] = next(file)[:16]

    # Correct sequencing errors in the previously generated sequence id: adapter dictionary
    print(f'Correcting sequencing data from {fastq}...')
    seq_true_dict = {seq_id: adapter_correcter(seq, adapter_array) for seq_id, seq in seq_dict.items()}

    # Assigns colour to each sequence ID based on the adapter if no adapter is found sequences are coloured black
    print(f'Assigning colours to sequence ids from {fastq}...')
    seq_id_colour_dict = {}
    for seq_id, seq in seq_true_dict.items():
        if seq is None:
            seq_id_colour_dict[seq_id] = '0,0,0'
        else:

            try:
                seq_id_colour_dict[seq_id] = adapter_dict[seq]

            except KeyError:
                seq_id_colour_dict[seq_id] = '0,0,0'

    # Name output file based on the input fastq file
    outfile = f'{fastq[:-6]}_index.csv'

    # Write the sequence ID : colour dictionary to a csv file
    with open(outfile, 'w') as csvfile:
        writer = csv.writer(csvfile)
        for key, value in seq_id_colour_dict.items():
            writer.writerow([key, value])


def compile_csv_data(directory):
    """
    Merges all intermediate CSV files in a directory
    :param directory: a directory containing intermediate csv files
    :return: CSV file
    """

    file_list = os.listdir(directory)  # List all files in a directory
    results = re.findall(r'fastq_[0-9]+_index.csv', str(file_list))  # filter list for all intermediate csv files
    results.sort()  # Order files numerically, ensures output file is consistent

    # Merge all intermediate csv files by loading each as a dictionary and writing it to the final output
    print(f'Compiling indexed sequences and writing to: {os.getcwd()}/multiprocess_index.csv')
    with open('multiprocess_index.csv', 'w') as output:
        writer = csv.writer(output)
        for file in results:
            with open(f'{directory}/{file}', 'r') as infile:
                subfile_dict = {rows[0]: rows[1] for rows in csv.reader(infile)}
                for key, value in subfile_dict.items():
                    writer.writerow([key, value])


@click.command()
@click.option('-f', '--fastq', type=str, required=True, help="Fastq file")
@click.option('-a', '--adapters', type=str, required=True, help="Adapter sequence list")
@click.option('-c', '--cores', type=int, required=False, help="Number of cores", default=1)
@click.option('-s', '--split', type=int, required=False, help="Maximum number of sequences per subfile")
def index(fastq, adapters, cores, split=1000000):
    """
    Click command to genrate a colour index of all reads in a fastq file based on adapter sequence of the reads
    :param fastq: linked-read fastq file to be processed
    :param adapters: list of all adapters that could be in linked-read fastq file
    :param cores: number of cores to use for multiprocessing
    :param split: number of sequences per subfile generated
    :return: csv file with each sequence ID assigned a colour based on adapter sequence
    """
    print("----- running NLR-Assembler index -----")
    # Ensure only available cores are used
    if cores > multiprocessing.cpu_count():
        print(f'Too many cores specified {cores}! {multiprocessing.cpu_count()} cores will be used')
        cores = multiprocessing.cpu_count()

    temp_dir = tempfile.mkdtemp()  # create a temporary directory
    print(f'temporary directory cretaed at: {temp_dir}')
    sub_fastq = split_fastq(fastq, split, temp_dir)  # split fastq file into intermediates for multiprocessing
    adapter_list, adapter_colour_dict = get_adapter_info(adapters)  # extract adapter sequence information
    pool = get_context("spawn").Pool(
        processes=cores)  # generate a pool for multiprocessing with specified number of cores
    pool.map(partial(index_subfile, adapter_list, adapter_colour_dict),
             sub_fastq)  # multiprocess fastq intermediate files
    compile_csv_data(temp_dir)  # merge all intermediate csv files into a single output
    shutil.rmtree(temp_dir)  # delete the temporary directory

# example command
# python3 main.py index -f test_data.fastq -a 4M-with-alts-february-2016.txt -c 1 -s 1000
