import logging

import click
import pandas as pd

from sml.calibrate.fitting import generate_lookup_function

@click.command()
@click.argument('raw_input', type=click.Path(exists=True))
@click.argument('output')
@click.option('-v', '--verbose', is_flag=True)
def main(raw_input, output, verbose):
    # set logger
    handler = logging.StreamHandler()
    formatter = logging.Formatter(
        '%(levelname).1s %(asctime)s [%(name)s] %(message)s', '%H:%M:%S'
    )
    handler.setFormatter(formatter)
    log_level = logging.WARNING if not verbose else logging.DEBUG
    logging.basicConfig(level=log_level, handlers=[handler])
    logger = logging.getLogger(__name__)

    data = pd.read_csv(raw_input)

    z, x, y = data['z'].values, data['x'].values, data['y'].values
    f, fx, fy, z0 = generate_lookup_function(z, x, y, model='huang', tol=1e-5)

    data['xo'] = fx(z)
    data['yo'] = fy(z)
    data['zcal'] = f(x, y)
    data['z0'] = data['z']-z0

    data.to_csv('result.csv', index=False, float_format='%.4f')

if __name__ == '__main__':
    main()
