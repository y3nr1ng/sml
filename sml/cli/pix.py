import logging
import click
import yaml
import imageio

from sml.analysis.analyzer import Analyzer

@click.command()
@click.argument('settings', type=click.Path(exists=True))
def main(settings):
    # Initiate the analysis pipline by parameters in SETTINGS.
    with open(settings, 'r') as fd:
        settings = yaml.safe_load(fd)

    # Set logger.
    handler = logging.StreamHandler()
    formatter = logging.Formatter(
        '%(levelname).1s %(asctime)s [%(name)s] %(message)s', '%H:%M:%S'
    )
    handler.setFormatter(formatter)
    log_level = logging.WARNING if settings['VERBOSE'] == 0 else logging.DEBUG
    logging.basicConfig(level=log_level, handlers=[handler])
    logger = logging.getLogger(__name__)

    # Start the image analysis.
    frames = imageio.volread(settings['INPF_RAW_IMAGE'])

    analyzer = Analyzer(settings)
    analyzer.process(frames)

    #DEBUG
    #from pprint import pprint
    #pprint(settings)
    print(frames.shape)

if __name__ == '__main__':
    main()
