import logging
import click
import pandas as pd
import numpy as np

from sml.calibrate.fitting import generate_lookup_function

#---------------------------------------------------------------------------
#
#       Output data format:
#
#       - z:    z-coordinate of raw data.
#       - x:    spot width in x direction of raw data.
#       - y:    spot width in y direction of raw data.
#       - xo:   calculated spot width in x direction.
#       - yo:   calculated spot width in y direction.
#       - zcal: calculated & center shifted z-coordinate.
#       - z0:   center shifted z-coordinate of raw data.
#
#---------------------------------------------------------------------------

@click.command()
@click.argument('raw_input', type=click.Path(exists=True))
@click.argument('output')
@click.option('-p', '--parafout', default='NULL',
              help='output filename of fitting parameters')
@click.option('-v', '--verbose', is_flag=True, help='verbose output')
def main(raw_input, output, parafout, verbose):
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
    f = generate_lookup_function(z, x, y, method='huang', model='huang',
                                 tol=1e-5)

    data['xo']   = f.fw(z)
    data['yo']   = f.fh(z)
    data['zcal'] = f(x, y)
    data['z0']   = data['z']-f.z0
    data.to_csv(output, index=False, float_format='%.4f')

    if (parafout != 'NULL'):
        para_name = ('w0', 'a', 'b', 'c', 'd')
        func_name = ('x', 'y')
        para_data = np.concatenate((f.fw.arguments, f.fh.arguments),
                                   axis=0).reshape(2,5)
        paras = pd.DataFrame(para_data, index=func_name, columns=para_name)
        paras.to_csv(parafout, index=True, float_format='%.8E')

if __name__ == '__main__':
    main()
