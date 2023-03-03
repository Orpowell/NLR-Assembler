#!/usr/local/bin/python3
import click
from commands.group import group
from commands.annotate import annotate
from commands.index import index


@click.group(help="CLI tool to manage full development cycle of projects")
def cli():
    pass


cli.add_command(index)
cli.add_command(group)
cli.add_command(annotate)

if __name__ == '__main__':
    cli()
