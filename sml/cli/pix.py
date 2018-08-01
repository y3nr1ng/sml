import logging

import click
import yaml

@click.command()
@click.argument('settings', type=click.Path(exists=True))
def main(settings):
    """Initiate the analysis pipline by parameters in SETTINGS."""
    with open(settings, 'r') as fd:
        settings = yaml.safe_load(fd)

    #DEBUG
    from pprint import pprint
    pprint(settings)

if __name__ == '__main__':
    main()
