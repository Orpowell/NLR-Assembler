#!/usr/local/bin/python3
import click
from commands.group_by_cosine import calculate_similarity
from commands.group_contigs import group_contigs


@click.group(help="CLI tool to manage full development cycle of projects")
def cli():
    pass


cli.add_command(calculate_similarity)
cli.add_command(group_contigs)

if __name__ == '__main__':
    cli()
