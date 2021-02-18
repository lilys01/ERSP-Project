"""Extract marker genes from Woltka-generated per-gene feature table.

Usage: parition(input.biom, gene_to_marker.map, output_dir)

Adapation of /projects/wol/woltka/extract_markers.py 

Output:
    One BIOM table per marker gene.
"""

from sys import argv
from os import makedirs
from os.path import join
from biom import Table, load_table
from biom.util import biom_open


def partition(biomTable, map, outputDir):

    # load input BIOM table
    table = biomTable
    samples = table.ids(axis='sample')
    genby = table.generated_by

    # read gene-to-marker mapping
    marker_dict = {}
    with open(map, 'r') as f:
        g = None
        for line in f:
            line = line.rstrip()
            if line.startswith('>'):
                g = line[1:]
                marker_dict[g] = {}
            else:
                idx, marker = line.split('\t')
                marker_dict[g][idx] = marker

    # extract marker genes from table
    tables = {}
    for gene, row in zip(table.ids(axis='observation'),
                         table.iter_data(dense=True, axis='observation')):
        g, _, idx = gene.partition('_')
        mdic = marker_dict[g]
        if idx in mdic:
            marker = mdic[idx]
            if marker not in tables:
                tables[marker] = [g], [row]
            else:
                tables[marker][0].append(g)
                tables[marker][1].append(row)

    # write per-marker tables
    outdir = outputDir
    makedirs(outdir, exist_ok=True)
    for marker, (gs, rows) in tables.items():
        with biom_open(join(outdir, f'{marker}.biom'), 'w') as f:
           # Table(rows, gs, samples).to_hdf5(f,'foobar') 
#replace line with import feature table (ar = qiime2.Artifact.import_data('FeatureTable[Frequency]', T) where T= biom.table(). then save (ar.save(name))
            ar = qiime2.Artifact.import_data('FeatureTable[Frequency]', f)
            ar.save('{marker}.qzv')

if __name__ == '__main__':
    main()
