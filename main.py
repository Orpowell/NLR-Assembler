#!/usr/local/bin/python3
import click
from commands import calculate_cosine
from commands import group_contigs


@click.group(help="CLI tool to manage full development cycle of projects")
def cli():
    pass


cli.add_command(calculate_cosine.calculate_similarity)
cli.add_command(group_contigs.group_contigs)

if __name__ == '__main__':
    cli()
