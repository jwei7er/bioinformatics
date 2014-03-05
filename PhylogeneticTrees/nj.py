# Written by Jordan Weiler
# March 2, 2014
#
# Python script for Phylogenetic Inference using Neighbor-Joining Algorithm
#  -Input:  reads file of pairwise distances
#  -Output: prints Newick tree representation for phylogenetic tree

from tree import Tree
import sys

# Global variables
species = []
pairwise_matrix = None

def ComputeNodeAvgDist():
    """
    Compute the average distance for each node to every other node
    """
    node_distances = []
    for i in range(len(species)):
        # Sum distance to other nodes / (#nodes - 2)
        sum_node = 0
        for j in range(len(species)):
            if j == i: continue

            if i < j: sum_node += pairwise_matrix[i][j]
            else:     sum_node += pairwise_matrix[j][i]
        
        if len(species)-2 <= 0:
            node_distances.append(sum_node)
        else:
            node_distances.append(sum_node / (len(species)-2))
    
    return node_distances

def FindLowestPair(species_distances):
    """
    Find lowest pairwise value in 2d pairwise matrix for NJ
    """
    min_calc, min_val, min_i, min_j = -1, -1, -1, -1

    # Spin through 2d matrix and find lowest pairwise distance minus average distances
    for i in range(len(species)-1):
        for j in range(i+1, len(species)):
            if pairwise_matrix[i][j] - species_distances[i] - species_distances[j] < min_calc or min_calc == -1:
                min_calc = pairwise_matrix[i][j] - species_distances[i] - species_distances[j]
                min_val = pairwise_matrix[i][j]
                min_i, min_j = i, j

    return min_val, min_i, min_j

def CreateJoinedDistances(a, b, min_val):
    """
    Create the combined distance for A and B with every other node for NJ
    """
    joined_values = []
    for i in range(len(species)):
        if i == a or i == b: continue
        
        # Get val for a species
        if i < a: val_a = pairwise_matrix[i][a]
        else:     val_a = pairwise_matrix[a][i]
        
        # Get val for b species
        if i < b: val_b = pairwise_matrix[i][b]
        else:     val_b = pairwise_matrix[b][i]
        
        # Average values
        joined_values.append((val_a + val_b - min_val)/2.0)

    return joined_values

def CombineSpecies(new_values, old_a, old_b):
    """
    Remove species A and B and combine them into a new AB species
    """
    # Add combined species to species list
    species.append("(" + species[old_a] + "),(" + species[old_b] + ")")
    
    # Remove joined species from species list  
    species.pop(old_b)
    species.pop(old_a)
    
    # Remove joined species from pairwise matrix (rows and columns)
    for row in pairwise_matrix:
        row.pop(old_b)
        row.pop(old_a)
    pairwise_matrix.pop(old_b)
    pairwise_matrix.pop(old_a)
    
    # Add joined species into pairwise matrix (joined values into column and entire blank bottom row)
    i = 0
    for row in pairwise_matrix:
        row.append(new_values[i])
        i+=1
    pairwise_matrix.append([-1 for i in range(len(species))])

def NJ_Algorithm():
    """
    Neighbor-Joining algorithm to join species and create phylogenetic tree
    """
    phylo_tree = None
    tree_dict = dict()

    # Loop through species until down to one pair
    for k in range(len(species)-1):
        # Get average distance for every node
        node_avg = ComputeNodeAvgDist()
        
        # Find lowest pair of species to join from average distance species list
        min_val, min_a, min_b = FindLowestPair(node_avg)
        
        if species[min_a] in tree_dict:
            # If the first joined species is already in the tree dictionary, retrieve it
            child1 = tree_dict[species[min_a]]
        else:
            # Else create new species tree
            child1 = Tree(species[min_a])
        # v_i = D_ij/2 + (u_i - u_j)/2
        child1.value = min_val/2.0 + (node_avg[min_a] - node_avg[min_b])/2.0
        
        if species[min_b] in tree_dict:
            # If the second joined species is already in the tree dictionary, retrieve it
            child2 = tree_dict[species[min_b]]
        else:
            # Else create new species tree
            child2 = Tree(species[min_b])
        # v_j = D_ij/2 + (u_j - u_i)/2
        child2.value = min_val/2.0 + (node_avg[min_b] - node_avg[min_a])/2.0
        
        # Create new tree with both A and B species joined 
        key = "(" + species[min_a] + "),(" + species[min_b] + ")"
        phylo_tree = Tree(key)
        phylo_tree.value = min_val/2.0

        # Put lower value to left of tree
        if child1.value <= child2.value:
            phylo_tree.AddLeftChild(child1, child1.value)
            phylo_tree.AddRightChild(child2, child2.value)
        else:
            phylo_tree.AddLeftChild(child2, child2.value)
            phylo_tree.AddRightChild(child1, child1.value)
        
        #phylo_tree.printNodeInfo()
        tree_dict[key] = phylo_tree
   
        # Create joined values to put back into 2d matrix 
        joinedDist = CreateJoinedDistances(min_a, min_b, min_val)
 
        # Remove A and B individually and insert AB joined
        CombineSpecies(joinedDist, min_a, min_b)

    return phylo_tree

if __name__ == "__main__":
    """
    The main function called when nj.py is run from the command line:

    > python nj.py distancefile
    """
    
    # Read in pairwise distance file
    if len(sys.argv) < 2:
        print ""
        print "ERROR! Pairwise distance file required as input"
        print ""
        sys.exit(1)
    
    # Spin through input file and store pairwise data to build 2d matrix later
    pairwise_data = []
    with open(sys.argv[1], "r") as fin:
        for line in fin:
            # Skip blank lines or lines starting with #
            if line[0] == "#" or line.strip() == "": continue
            
            sl = line.split()
            if sl[0] not in species:
                species.append(sl[0])
            if sl[1] not in species:
                species.append(sl[1])
            pairwise_data.append(sl)
    #print species
    
    # Initialize 2d pairwise matrix to size of species list
    pairwise_matrix = [[-1 for j in range(len(species))] for i in range(len(species))]
    
    # Spin through pairwise data and build 2d matrix
    for line in pairwise_data:
        sp1 = species.index(line[0])
        sp2 = species.index(line[1])
        if sp1 < sp2:
            i, j = sp1, sp2
            key = line[0] + "," + line[1]
        else:
            i, j = sp2, sp1
            key = line[1] + "," + line[0]
        pairwise_matrix[i][j] = float(line[2])
    #print pairwise_matrix
    
    # Run NJ algorithm against distance matrix
    phylo_tree = NJ_Algorithm()
    
    # Final output in Newick notation
    # i.e. (((A:1, E:1):2, C:3):1, (D:1.5,B:1.5):2.5)
    print phylo_tree.BuildNewick()
  

