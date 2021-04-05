import qiime2
from qiime2.plugins import diversity, feature_table, emperor, phylogeny
from qiime2 import Artifact
#setting data variables
fTable = Artifact.load('../feature-table.qza')
tree = Artifact.load('../tree.qza')
metadata = qiime2.Metadata.load('../sample-metadata.txt')

#Distance Matrix
distMatrix = diversity.pipelines.beta_phylogenetic(table=fTable, phylogeny=tree, metric='unweighted_unifrac')
dm = distMatrix.distance_matrix
dm.save('taxonomic-dist-matrix.qza')


#obtaining PCoA
PCoA = diversity.actions.pcoa(distance_matrix=dm)
pcoaResults = PCoA.pcoa

#emperor visulization
viz = emperor.actions.plot(pcoaResults, metadata)
viz.visualization.save('emperor_plot.qzv')

#permanova
columns = ['age_cat','antibiotic_history','types_of_plants']
for dataType in columns: 
	data = metadata.get_column(dataType)
	results = diversity.visualizers.beta_group_significance(distance_matrix=dm, metadata=data, pairwise=True)
	results.visualization.save(dataType)


	
