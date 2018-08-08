import logging

import click
import pandas as pd

from sml.calibrate.fitting import fit_curve

@click.command()
@click.argument('raw_input', type=click.Path(exists=True))
@click.argument('output')
@click.option('-v', '--verbose', is_flag=True)
def main(raw_input, output, verbose):
    data = pd.read_csv(raw_input)
    p = fit_curve(data['z'].values, data['x'].values)
    print(p)
    
if __name__ == '__main__':
    main()
