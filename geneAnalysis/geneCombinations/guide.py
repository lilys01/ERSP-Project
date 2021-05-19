import qiime2
from qiime2.plugins import  diversity, feature_table, emperor, phylogeny, taxa
from qiime2.plugins.diversity.visualizers import mantel
from qiime2.plugins.feature_table.methods import rarefy
from qiime2 import Artifact, Metadata
import biom, os, glob, skbio, csv
import pandas as pd
import numpy as np
import re 

import unifrac
from unifrac import _methods
from unifrac._meta import CONSOLIDATIONS
#from partitioning import partition

#importing data
fTable = Artifact.load('../genes-ftable.qza')
speciesTree = Artifact.load('../../tree.qza')
mdata = Metadata.load('../../sample-metadata.txt')

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



'''
compute unweighted unifrac on each gene with it's table and tree
compute pcoa and make an emperor plot for each
add a new column for gene read counts to the metadata
group covid+ metadata together
perform mantel test on each gene against species tree
 '''

def geneCombinationAnalysis(keysString, tablelookup, treelookup): 
   
   tables_list =[]
   trees_list = []
   rarefyFts = []

   # create a list from the passed string tuple
   keysList = re.findall(r"'(.*?)'",keysString) 

   # create a passable string in format 'p0abc-p0def-...-p0xyz' 
   name = '-'.join(map(str,keysList))
    
   # create list of trees and tables of the genes in the tuples
   for key in keysList:
      tables_list.append("rarefiedTables/" + key + ".biom")
      trees_list.append("../trees/" + key + ".nwk")
   
   tables_tuple = tuple(tables_list)
   trees_tuple = tuple(trees_list)   
   print(tables_tuple)
   #unifrac
   dMatrix = unifrac.meta(tables=tables_tuple, phylogenies=trees_tuple, method ='unweighted')
   dMatrix.write('/panfs/panfs1.ucsd.edu/panscratch/jvsantan/ERSP-Project/geneCombinations/dms/' + name + '-dm.qza')   
   dm = qiime2.Artifact.import_data('DistanceMatrix', dMatrix)   

   PCoA = diversity.actions.pcoa(distance_matrix=dm)
   pcoaResults = PCoA.pcoa 

   df = mdata.to_dataframe()

   #add new metadata grouping together all positive related covid metadata
   positives = []
   for sample, row in df.iterrows():
     resp = row['covid_suspected_positive'] 
     resp = str(resp)
     resp = resp.split()[0]
     if (resp == 'Yes,'):
       positives.append('Yes')
     else:
       positives.append('No')
   df['covid_positive'] = positives

   #make new dataframe metadata object
   new_md = qiime2.Metadata(df)

   viz = emperor.actions.plot(pcoaResults, new_md)
   viz.visualization.save('/panfs/panfs1.ucsd.edu/panscratch/jvsantan/ERSP-Project/geneCombinations/emp-plots/' + name + '-emperor.qzv')

   
   #permanova on each gene table with new metadata column we made
   permanovaSkbio(new_md, dm, 'covid_positive', name)


   #mantel test on each gene against species tree (tree var imported above)
   speciesDm = qiime2.Artifact.load('../../taxonomicAnalysis/taxonomic-dist-matrix.qza')
   #mantelTestQiime(dm, speciesDm, label=True, name1=key, name2='species', intersectIds=True)  
   mantelTestSkbio(dm, speciesDm, name)

   

'''
creates two dictionaries to map gene names to their partitioned gene tables and their respective trees
'''
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
   return tablelookup, treelookup

'''
generates lists of unique tuples
'''
def uniqueTuples():
   #import list of genes in p0XZY format accounting for the 19 missing ones
   genes = open("genesIds.txt").read()
   genes_list = genes.split("\n")
   genes_list.remove('')

   #initialize a list
   list =[]

   #generate unique tuples
   for k in range(2,11):
      for i in range(10):
         list.append(sample(genes_list,k))

   return list

'''
save a summary of a feature table to help determine filtering
'''
def summarize(newtable, mdata, key):
   summary = feature_table.actions.summarize(table=newtable, sample_metadata=mdata)
   summary.visualization.save(key + '_summary.qzv')




'''
uses permanova qiime2 plugin to generate a permanova visualization
categories - a list of strings representing the columns in metadata
mdata - qiime2 metadata
dm - qiime2 artifact, distance matrix
geneName - string, the gene ID
'''
def permanovaVis(categories, mdata, dm, geneName):
   for dataType in categories: 
      data = mdata.get_column(dataType)
      results = diversity.visualizers.beta_group_significance(distance_matrix=dm, metadata=data, pairwise=True)
      results.visualization.save('pnova-results/'+ dataType + '_' + geneName + '_pnova')



'''
computes permanova using skbio and outputs a csv file with the results
mdata - qiime2 metadata object
dm - qiime2 artifact, distance matrix
columnName - string, column in metadata
geneName - string, gene ID

'''
def permanovaSkbio(mdata, dm, columnName, geneName):
   dm = dm.view(skbio.DistanceMatrix)
   mdata = mdata.to_dataframe()
   pnovaResult = skbio.stats.distance.permanova(dm, mdata, column=columnName)
   #pnovaResult.to_csv('pnova-results/'+ geneName + '-' + columnName + '-pnova')
   tStat = pnovaResult.get(key='test statistic') 
   pValue = pnovaResult.get(key='p-value')
   nGroups = pnovaResult.get(key='number of groups')
   with open('pnovaResultsGC.csv','a',newline='') as file:
     writer = csv.writer(file)
     writer.writerow([geneName, tStat, pValue, nGroups])




'''
qiime2 mantel test on two distance matrices, 
dm1,dm2 - qiime2 distance matrix
label - bool, needs to be true to apply labels to visualization
name1,name2 - the labels for visualization
intersectIds - bool, will filter out unmatched ids between dm's
'''
def mantelTestQiime(dm1, dm2, label, name1, name2, intersectIds):
   #dm1 = qiime2.Artifact.load(dm1)
   #dm2 = qiime2.Artifact.load(dm2)
   name1 = str(name1) 
   if(label == True):
      mant = mantel(dm1, dm2, label1=name1,label2=name2, intersect_ids=intersectIds)
      mant, = mant
      mant.save(name1+'-'+name2+'-mantel.qzv')
   else:
      mant = mantel(dm1,dm2)
      mant, = mant
      mant.save('mantel.qzv')




'''
mantel test using skbio
'''
def mantelTestSkbio(dm1, dm2, geneName):
   dm1 = dm1.view(skbio.DistanceMatrix)
   dm2 = dm2.view(skbio.DistanceMatrix)  
   coeff, p_value, n = skbio.stats.distance.mantel(dm1, dm2, method='pearson', permutations=999, alternative='two-sided', strict=False) 
   with open('mantelResultsGC.csv','a',newline='') as file:
     writer = csv.writer(file)
     writer.writerow([geneName, coeff, p_value, n])

'''
computes unifrac and saves file 
params:
fttable - feature table qiime artifact 
tree - corresponding phylogeny qiime artifact
metric - unweighted_unifrac or weighted_unifrac 

'''
def our_unifrac(fttable, tree, type):
   distMatrix = diversity.pipelines.beta_phylogenetic(table=fttable, phylogeny=tree, metric=type)
   dm = distMatrix.distance_matrix
   dm.save('gene-dms/' + key + '-dist-matrix.qza')
   return dm 








