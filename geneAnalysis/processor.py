import click
from analysis import perGeneAnalysis, setDictionaries


@click.command()
@click.option('--parametersfile', type=click.Path(exists=True), required=True)
@click.option('--taskid', type=int, required=True)

def process(parametersfile, taskid):
    details = open(parametersfile).readlines()
    geneName = details[taskid - 1]
    geneName = ''.join(geneName.split())
    tablelookup, treelookup = setDictionaries()
    perGeneAnalysis(geneName, tablelookup, treelookup)


if __name__ == '__main__':
    process()


