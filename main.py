#!/usr/local/bin/python3
import click
from commands.group import group
from commands.index import index
from commands.nlr_coverage import nlr_coverage
from commands.query_coverage import query_coverage


@click.group(help="CLI tool to manage full development cycle of projects")
def cli():
    pass


cli.add_command(index)
cli.add_command(group)
cli.add_command(nlr_coverage)
cli.add_command(query_coverage)


if __name__ == '__main__':
    cli()
