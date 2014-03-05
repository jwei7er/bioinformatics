# Written by Jordan Weiler
# March 2, 2014
#
# Python script for Phylogenetic Inference using UPGMA Algorithm
#  -Input:  reads file of pairwise distances
#  -Output: prints Newick tree representation for phylogenetic tree

from tree import Tree
import sys

# global variables
species = []
pairwise_matrix = None

def FindLowestPair():
    """
    Find lowest pairwise value in 2d pairwise matrix
    """
    min_val, min_i, min_j = -1, -1, -1

    # Spin through 2d matrix and find lowest pairwise distance
    for i in range(len(species)-1):
        for j in range(i+1, len(species)):
            if pairwise_matrix[i][j] < min_val or min_val == -1:
                min_val = pairwise_matrix[i][j]
                min_i, min_j = i, j

    return min_val, min_i, min_j

def CreateJoinedDistances(a, b):
    """
    Create the combined distance for A and B with every other node
    """
    joined_values = []
    for i in range(len(species)):
        if i == a or i == b: continue
    
        # Get val for i species
        if i < a: val_a = pairwise_matrix[i][a]
        else:     val_a = pairwise_matrix[a][i]
    
        # Get val for j species
        if i < b: val_b = pairwise_matrix[i][b]
        else:     val_b = pairwise_matrix[b][i]
    
        # Take the average of the values
        joined_values.append((val_a + val_b)/2.0)

    return joined_values

def CombineSpecies(new_values, old_a, old_b):
    """
    Remove species A and B and combine them into a new AB species
    """
    # Add combined species into species list
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

def UPGMA_Algorithm():
    """
    UPGMA algorithm to join species and create phylogenetic tree
    """
    phylo_tree = None
    tree_dict = dict()

    # Loop through species until down to one pair
    for k in range(len(species)-1):
        min_val, min_i, min_j = FindLowestPair()
    
        if species[min_i] in tree_dict:
            # If the first joined species is already in the tree dictionary, retrieve it
            child1 = tree_dict[species[min_i]]
        else:
            # Else create new species tree
            child1 = Tree(species[min_i])
            child1.value = min_val/2.0
    
        if species[min_j] in tree_dict:
            # If the second joined species is already in the tree dictionary, retrieve it
            child2 = tree_dict[species[min_j]]
        else:
            # Else create new species tree
            child2 = Tree(species[min_j])
            child2.value = min_val/2.0
     
        # Create new tree with both A and B species joined 
        key = "(" + species[min_i] + "),(" + species[min_j] + ")"
        phylo_tree = Tree(key)
        phylo_tree.value = min_val/2.0

        # Put lower value to left of tree
        if child1.value <= child2.value:
            phylo_tree.AddLeftChild(  child1, (phylo_tree.value - child1.value) )
            phylo_tree.AddRightChild( child2, (phylo_tree.value - child2.value) )
        else:
            phylo_tree.AddLeftChild(  child2, (phylo_tree.value - child2.value) )
            phylo_tree.AddRightChild( child1, (phylo_tree.value - child1.value) )
        
        #phylo_tree.printNodeInfo()
        tree_dict[key] = phylo_tree
     
        # Create joined values to put back into 2d matrix 
        joined_values = CreateJoinedDistances(min_i, min_j)
        
        # Remove A and B individually and insert AB joined
        CombineSpecies(joined_values, min_i, min_j)

    return phylo_tree

if __name__ == "__main__":
    """
    The main function called when upgma.py is run from the command line:

    > python upgma.py distancefile
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
            i,j = sp1, sp2
            key = line[0] + "," + line[1]
        else:
            i,j = sp2, sp1
            key = line[1] + "," + line[0]
        pairwise_matrix[i][j] = float(line[2])
    #print pairwise_matrix
 
    # Run UPGMA algorithm against distance matrix 
    phylo_tree = UPGMA_Algorithm()
    
    # Final output in Newick notation
    # i.e. (((A:1, E:1):2, C:3):1, (D:1.5,B:1.5):2.5)
    print phylo_tree.BuildNewick()
  
