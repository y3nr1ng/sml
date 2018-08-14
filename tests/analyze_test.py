import logging

import click
import imageio
import pandas as pd

from sml.analysis.analyzer import Analyzer

# set logger
handler = logging.StreamHandler()
formatter = logging.Formatter(
    '%(levelname).1s %(asctime)s [%(name)s] %(message)s', '%H:%M:%S'
)
handler.setFormatter(formatter)
logging.basicConfig(level=logging.DEBUG, handlers=[handler])
logger = logging.getLogger(__name__)

@click.command()
@click.argument('path', type=click.Path(exists=True))
def main(path):
    frames = imageio.volread(path)
    print(data.shape)

    #TODO fetch data to analyzer
    analyzer = Analyzer()

    analyzer.process(frames)

if __name__ == '__main__':
    main()
