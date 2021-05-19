import click
from guide import geneCombinationAnalysis, setDictionaries


@click.command()
@click.option('--parametersfile', type=click.Path(exists=True), required=True)
@click.option('--taskid', type=int, required=True)

def process(parametersfile, taskid):
    details = open(parametersfile).readlines()
    geneTuple = details[taskid - 1 + 8000]
    #geneTuple = ''.join(geneTuple.split())
    
    tablelookup, treelookup = setDictionaries()
    geneCombinationAnalysis(geneTuple, tablelookup, treelookup)


if __name__ == '__main__':
    process()


