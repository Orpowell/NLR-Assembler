#!/usr/local/bin/python3
import click
from commands.group import group
from commands.annotate import annotate
from commands.index import index
from commands.coverage import coverage


@click.group(help="CLI tool to manage full development cycle of projects")
def cli():
    pass


cli.add_command(index)
cli.add_command(group)
cli.add_command(annotate)
cli.add_command(coverage)

if __name__ == '__main__':
    cli()
