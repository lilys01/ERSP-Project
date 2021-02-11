import qiime2
from qiime2.plugins import biodiversity, feature_table, emperor, phylogeny, taxa, empress
from qiime2 import Artifact, Metadata

#importing data
fTable = Artifact.load('genes-ftable.qza')
tree = Artifact.load('../tree.qza')
mdata = Metadata.load('../sample-metadata.txt')

#filtering feature table
