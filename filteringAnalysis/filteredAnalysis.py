import qiime2
from qiime2.plugins import diversity, feature_table, emperor, phylogeny, taxa, empress
from qiime2 import Artifact

#setting data variables
fTable = Artifact.load('../feature-table.qza')
tree = Artifact.load('../tree.qza')
mdata = qiime2.Metadata.load('../sample-metadata.txt')

#filtering feature table
filtResult = feature_table.actions.filter_samples(min_frequency=1000, metadata=mdata, where="[env_material]='feces'",table=fTable)
filtTable = filtResult.filtered_table

#summarize table
summary = feature_table.actions.summarize(table=filtTable, sample_metadata=mdata)
summary.visualization.save('filtTable-summary.qzv')

#create taxa bar plots
taxonomyInfo = Artifact.load('taxonomy.qza')
plots = taxa.actions.barplot(table=filtTable, metadata=mdata, taxonomy=taxonomyInfo)
plots.visualization.save('taxa-bar-plots.qzv')

#Distance Matrix
distMatrix = diversity.pipelines.beta_phylogenetic(table=filtTable, phylogeny=tree, metric='unweighted_unifrac')
dm = distMatrix.distance_matrix

#obtaining PCoA
PCoA = diversity.actions.pcoa(distance_matrix=dm)
pcoaResults = PCoA.pcoa

#emperor visulization
viz = emperor.actions.plot(pcoaResults, mdata)
viz.visualization.save('emperor_plot.qzv')

#permanova
columns = ['age_cat','antibiotic_history','types_of_plants']
for dataType in columns: 
	data = mdata.get_column(dataType)
	results = diversity.visualizers.beta_group_significance(distance_matrix=dm, metadata=data, pairwise=True)
	results.visualization.save(dataType)

#create empire plot
fmetadata = taxonomyInfo.view(qiime2.Metadata)
empire = empress.actions.community_plot(tree=tree,sample_metadata=mdata, feature_metadata=fmetadata, feature_table=filtTable, pcoa=pcoaResults)
empire.visualization.save('empire_plot.qzv')


