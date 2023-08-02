#!/usr/local/bin/python3
import click
import logging
import sys
from commands.group import group
from commands.index import index
from commands.nlr_coverage import nlr_coverage
from commands.contig_coverage import contig_coverage

logging.basicConfig(stream=sys.stdout, format='%(asctime)s - %(message)s', datefmt='%d-%b-%y %H:%M:%S',
                    level=logging.INFO)

@click.group(help="CLI tool to manage full development cycle of projects")
def cli():
    pass


cli.add_command(index)
cli.add_command(group)
cli.add_command(nlr_coverage)
cli.add_command(contig_coverage)


if __name__ == '__main__':
    cli()
