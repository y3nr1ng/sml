import click

@click.command()
@click.argument('raw_input', type=click.Path(exists=True),
                help="A list of ellipical axis lengths for X and Y, with respect to z position.")
@click.argument('output',
                help="Fitting results of the calibration curve for both X and Y axes.")
@click.option('-v', '--verbose', is_flag=True)
def main(raw_input, output, verbose):
    print(verbose)

if __name__ == '__main__':
    main()
