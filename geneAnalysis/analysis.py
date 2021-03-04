import qiime2
from qiime2.plugins import diversity, feature_table, emperor, phylogeny, taxa
from qiime2 import Artifact, Metadata
import biom
import os
import glob
import pandas as pd
from partitioning import partition

#importing data
fTable = Artifact.load('genes-ftable.qza')
tree = Artifact.load('../tree.qza')
mdata = Metadata.load('../sample-metadata.txt')

'''
#filtering feature table
filtResult = feature_table.actions.filter_samples(min_frequency=100000, metadata=mdata, where="[env_material]='feces'",table=fTable)
filtTable = filtResult.filtered_table
#filtTable.save('filtered-gene-table.qza')

#paritioning 
btable = filtTable.view(biom.Table)
partition(btable,'/projects/wol/release/markers/phylophlan.intree.map', 'gene-tables') 
'''
#set dictionaries with gene-tree and gene-table pairs for easy look up

tablelookup = {}
for f in glob.glob('gene-tables/*.qza'):
 basename = os.path.basename(f)
 no_extension = os.path.splitext(basename)[0]
 tablelookup[no_extension] = f

treelookup = {}
for t in glob.glob('trees/*.nwk'):
 base = os.path.basename(t)
 no_ext = os.path.splitext(base)[0]
 treelookup[no_ext] = t 

'''for key in tablelookup.keys():'''
key='p0097'
newtable = Artifact.load(tablelookup.get(key))
newtree = qiime2.Artifact.import_data('Phylogeny[Rooted]',treelookup.get(key))

distMatrix = diversity.pipelines.beta_phylogenetic(table=newtable, phylogeny=newtree, metric='unweighted_unifrac')
dm = distMatrix.distance_matrix

PCoA = diversity.actions.pcoa(distance_matrix=dm)
pcoaResults = PCoA.pcoa 

'''
viz = emperor.actions.plot(pcoaResults, mdata)
viz.visualization.save('p0097_emperor.qzv')

summary = feature_table.actions.summarize(table=newtable, sample_metadata=mdata)
summary.visualization.save('p0097_summary.qzv')
'''

'''Add a new category column to our metadata'''


'''df = pd.DataFrame('../sample-metadata.txt')'''
df = mdata.to_dataframe()

'''bTable = biom.load_table('gene-tables/p0097.biom')'''
bTable = qiime2.Artifact.load('gene-tables/p0097.qza').view(biom.Table)

newCategory = pd.Series(bTable.sum('sample'), index=bTable.ids())
df['p0097-read-counts'] = newCategory

new_md = qiime2.Metadata(df)

viz = emperor.actions.plot(pcoaResults, new_md)
viz.visualization.save('p0097_emperor.qzv')

summary = feature_table.actions.summarize(table=newtable, sample_metadata=new_md)
summary.visualization.save('p0097_summary.qzv')





 
