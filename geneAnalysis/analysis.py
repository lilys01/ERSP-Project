import qiime2
from qiime2.plugins import diversity, feature_table, emperor, phylogeny, taxa
from qiime2.plugins.feature_table.methods import rarefy
from qiime2 import Artifact, Metadata
import biom
import os
import glob
import pandas as pd
import numpy as np 
from partitioning import partition

#importing data
fTable = Artifact.load('genes-ftable.qza')
tree = Artifact.load('../tree.qza')
mdata = Metadata.load('../sample-metadata.txt')

'''
Filters and partitions tables
'''
def gentables():
    
    #filtering feature table
    filtResult = feature_table.actions.filter_samples(min_frequency=100000, metadata=mdata, where="[env_material]='feces'",table=fTable)
    filtTable = filtResult.filtered_table

    filtTable.save('filtered-gene-table.qza')

    #paritioning 
    btable = filtTable.view(biom.Table)
    partition(btable,'/projects/wol/release/markers/phylophlan.intree.map', 'gene-tables') 



#set dictionaries with gene-tree and gene-table pairs for easy look up
def setDictionaries():
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


'''

compute unweughted unifrac on each gene with it's table and tree
compute pcoa and make an emperor plot for each
add a new column for gene read counts to the metadata

 '''

def perGeneAnalysis(key):
   #print(key) 
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

   newtable = Artifact.load(tablelookup.get(key))
   newtree = qiime2.Artifact.import_data('Phylogeny[Rooted]',treelookup.get(key))

   bTable = qiime2.Artifact.load(tablelookup.get(key)).view(biom.Table)

   #rarefraction on gene table
   #quartile,bins = pd.qcut(bTable.sum('sample'), 4, retbins=True)
   quartile = pd.qcut(bTable.sum('sample'),4,labels=False)
   rarefyFt, = feature_table.methods.rarefy(newtable, int(quartile[3]))

   #unifrac
   distMatrix = diversity.pipelines.beta_phylogenetic(table=rarefyFt, phylogeny=newtree, metric='unweighted_unifrac')
   dm = distMatrix.distance_matrix
   dm.save('gene-dms/' + key + '-dist-matrix.qza')

   PCoA = diversity.actions.pcoa(distance_matrix=dm)
   pcoaResults = PCoA.pcoa 


   #add column of gene counts to metadata

   df = mdata.to_dataframe()

   newCategory = pd.Series(bTable.sum('sample'), index=bTable.ids())

   columnName = key + '-read-counts'

   df[columnName] = newCategory

   new_md = qiime2.Metadata(df)

   viz = emperor.actions.plot(pcoaResults, new_md)
   viz.visualization.save('gene-emp-plots/' + key + '-emperor.qzv')
  

def summarize():
   summary = feature_table.actions.summarize(table=newtable, sample_metadata=mdata)
   summary.visualization.save('p0097_summary.qzv')




