#!/usr/local/bin/python3
import click
from commands import similarity


@click.group(help="CLI tool to manage full development cycle of projects")
def cli():
    pass


cli.add_command(similarity.calculate_similarity)

if __name__ == '__main__':
    cli()
