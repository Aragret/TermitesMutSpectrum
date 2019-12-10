library(ape)

tree = read.tree('Alina/other_projects/TermitesMutSpectrum/AllTermites/data/cockroaches11_19/forAlign/results/phylogeny/iqtree_mito.treefile')

is.rooted(tree)

rooted.tree = root(tree, 'Locusta_migratoria_tibetensis_1', resolve.root = TRUE)

is.rooted(rooted.tree)

write.tree(rooted.tree, 'Alina/other_projects/TermitesMutSpectrum/AllTermites/results/rootedTreeCockroaches.newick')
