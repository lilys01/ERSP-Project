import qiime2
from qiime2.plugins import feature_table, emperor, phylogeny, taxa, empress
from qiime2 import Artifact, Metadata
import biom
import os
from partitioning import partition

#importing data
fTable = Artifact.load('genes-ftable.qza')
tree = Artifact.load('../tree.qza')
mdata = Metadata.load('../sample-metadata.txt')

#filtering feature table
filtResult = feature_table.actions.filter_samples(min_frequency=100000, metadata=mdata, where="[env_material]='feces'",table=fTable)
filtTable = filtResult.filtered_table
#filtTable.save('filtered-gene-table.qza')

#paritioning 
btable = filtTable.view(biom.Table)
partition(btable,'/projects/wol/release/markers/phylophlan.intree.map', 'gene-tables') 


