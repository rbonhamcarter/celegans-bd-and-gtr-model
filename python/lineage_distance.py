''' 
This script imports all trees in the outputs folder and plots a histrogram showing
the distribution of all pairwise lineage distances in the set of trees.

Trees are expected to be in .nex format.
'''
import os
import dendropy
import numpy as np
import matplotlib.pyplot as plt

OUTPUT_DIR = "./output/"
TRUE_TREE_FILEPATH = None    #  EDIT IF TRUE TREE FILE EXISTS


def parse_trees_paths_from_directory(directory):
    source_files = []
    for item in os.listdir(directory):
        if "tree" in item:
            path = os.path.join(directory, item)
            source_files.append(path)

    return source_files

def get_lineage_dist_histogram_from_sim_trees(source_files, taxon_namespace):
    # Compute the lineage distances for each tree
    lineage_distances = []
    for path in source_files:
        tree = dendropy.Tree.get(path=path, 
                                schema="nexus",
                                taxon_namespace=taxon_namespace,
                                store_tree_weights=True,
                                rooting="force-rooted")
        pdm = tree.phylogenetic_distance_matrix()
        lineage_distances.extend(pdm.distances())

    # Compute the histogram
    counts, bins = np.histogram(lineage_distances, bins=100, density=False)
    num_trees = len(source_files)
    normed_counts = [i/num_trees for i in counts] # Normalize by number of trees

    return normed_counts, bins

def get_lineage_dist_histogram_from_true_tree(true_tree_path):
    # Compute the lineage distances
    lineage_distances = []
    tree = dendropy.Tree.get(path=true_tree_path, 
                            schema="newick",
                            store_tree_weights=True,
                            rooting="force-rooted")
    pdm = tree.phylogenetic_distance_matrix()
    lineage_distances = pdm.distances()
    total_time = tree.max_distance_from_root()

    # Compute the histogram
    counts, bins = np.histogram(lineage_distances, bins=100, density=False)

    return counts, bins

def main():
    source_files = parse_trees_paths_from_directory(OUTPUT_DIR)
    taxon_namespace = dendropy.TaxonNamespace()

    # Set up the plot
    if TRUE_TREE_FILEPATH == None:
        fig, ax = plt.subplots(1, 1, figsize=(8, 4))    
        ax.set_xlabel("Pairwise Lineage Distance")
        ax.set_ylabel("Average Frequency per Tree")

        # Get and plot the histrogram for simulated trees
        sim_normed_counts, sim_bins = get_lineage_dist_histogram_from_sim_trees(source_files, taxon_namespace)
        sim_width = np.max(sim_bins)/np.size(sim_bins) 
        ax.bar(sim_bins[:-1], sim_normed_counts, width=sim_width, alpha=0.5, color='steelblue')

    else:
        fig, ax = plt.subplots(1, 2, figsize=(8, 4))    
        ax[0].set_xlabel("Pairwise Lineage Distance")
        ax[0].set_ylabel("Average Frequency per Tree")
        ax[1].set_xlabel("Pairwise Lineage Distance")
        ax[1].set_xlabel("Frequency")

        # Get and plot the two histrograms (one for simulated trees, one for true tree)
        sim_normed_counts, sim_bins = get_lineage_dist_histogram_from_sim_trees(source_files, taxon_namespace)
        sim_width = np.max(sim_bins)/np.size(sim_bins) 
        ax[0].bar(sim_bins[:-1], sim_normed_counts, width=sim_width, alpha=0.5, color='steelblue')

        true_counts, true_bins = get_lineage_dist_histogram_from_true_tree(TRUE_TREE_FILEPATH)
        true_width = np.max(true_bins)/np.size(true_bins) 
        ax[1].bar(true_bins[:-1], true_counts, width=true_width, alpha=0.5, color='steelblue')
    
    plt.show()


if __name__ == "__main__":
    main()