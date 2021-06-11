import qiime2
from qiime2.plugins import diversity, feature_table, emperor, phylogeny, taxa
from qiime2.plugins.diversity.visualizers import mantel
from qiime2.plugins.feature_table.methods import rarefy
from qiime2 import Artifact, Metadata
import biom, os, glob, skbio, csv
import pandas as pd
import numpy as np 
from partitioning import partition

#importing data
fTable = Artifact.load('genes-ftable.qza')
speciesTree = Artifact.load('../tree.qza')
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



'''
compute unweighted unifrac on each gene with it's table and tree
compute pcoa and make an emperor plot for each
add a new column for gene read counts to the metadata
group covid+ metadata together
perform mantel test on each gene against species tree
 '''

def perGeneAnalysis(key, tablelookup, treelookup): 
   
   newtable = Artifact.load(tablelookup.get(key))
   newtree = qiime2.Artifact.import_data('Phylogeny[Rooted]',treelookup.get(key))

   bTable = qiime2.Artifact.load(tablelookup.get(key)).view(biom.Table)

   #rarefraction on gene table
   quartile = pd.qcut(bTable.sum('sample'),4)
   rarefyFt, = feature_table.methods.rarefy(newtable, int(quartile.categories[-1].left))

   #unifrac
   distMatrix = diversity.pipelines.beta_phylogenetic(table=rarefyFt, phylogeny=newtree, metric='unweighted_unifrac')
   dm = distMatrix.distance_matrix
   dm.save('gene-dms/' + key + '-dist-matrix-unfilt.qza')

   PCoA = diversity.actions.pcoa(distance_matrix=dm)
   pcoaResults = PCoA.pcoa 


   #add column of gene counts to metadata

   df = mdata.to_dataframe()

   newCategory = pd.Series(bTable.sum('sample'), index=bTable.ids())

   columnName = key + '-read-counts'

   df[columnName] = newCategory
   

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
   viz.visualization.save('gene-emp-plots/' + key + '-emperor-unfilt.qzv')

   
   #permanova on each gene table with new metadata column we made
   #permanovaSkbio(new_md, dm, 'covid_positive', key)


   #mantel test on each gene against species tree (tree var imported above)
   #speciesDm = qiime2.Artifact.load('../taxonomicAnalysis/taxonomic-dist-matrix.qza')
   #mantelTestQiime(dm, speciesDm, label=True, name1=key, name2='species', intersectIds=True)  
   #mantelTestSkbio(dm, speciesDm, key)



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
      results.visualization.save(dataType + '_' + geneName + '_pnova')



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
   with open('pnovaResults.csv','a',newline='') as file:
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
   with open('mantelResults.csv','a',newline='') as file:
     writer = csv.writer(file)
     writer.writerow([geneName, coeff, p_value, n])

'''
computes unifrac and saves file 
params:
fttable - feature table qiime artifact 
tree - corresponding phylogeny qiime artifact
metric - unweighted_unifrac or weighted_unifrac 

'''
def unifrac(fttable, tree, type):
   distMatrix = diversity.pipelines.beta_phylogenetic(table=fttable, phylogeny=tree, metric=type)
   dm = distMatrix.distance_matrix
   dm.save('gene-dms/' + key + '-dist-matrix.qza')
   return dm 








