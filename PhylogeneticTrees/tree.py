# Written by Jordan Weiler
# March 2, 2014
#
# Tree class used to build phylogenetic tree

class Tree:
    def __init__(self, name):
        self.name = name
        self.value = 0.0
        self.left = None
        self.left_value = 0.0
        self.right = None
        self.right_value = 0.0
        self.parent = None
      
    def AddLeftChild(self, tree, value):
        """
        Add left child to tree
        """
        self.left = tree
        self.left_value = value
        tree.parent = self

    def AddRightChild(self, tree, value):
        """
        Add right child to tree
        """
        self.right = tree
        self.right_value = value
        tree.parent = self

    def PrintNodeInfo(self):
        """
        Print node information for debugging
        """
        print "Name:", self.name
        print "Value:", self.value
        print "Left:", self.left
        print "Left Value:", self.left_value
        print "Right:", self.right
        print "Right Value:", self.right_value

    def BuildNewick(self):
        """
        Recursively build Newick representation of tree
        i.e. (((A:1, E:1):2, C:3):1, (D:1.5,B:1.5):2.5)
        """
        if self.left is None and self.right is None:
            return self.name + ":" + str(self.value)
        elif self.left.left is None and self.right.right is None:
            return "(" + self.left.BuildNewick() + ", " + self.right.BuildNewick() + ")"
        elif self.left.left is None:
            return "(" + self.left.BuildNewick() + ", " + self.right.BuildNewick() + ":" + str(self.right_value) + ")"
        elif self.right.right is None:
            return "(" + self.left.BuildNewick() + ":" + str(self.left_value) + ", " + self.right.BuildNewick() + ")"
        else:
            return "(" + self.left.BuildNewick() + ":" + str(self.left_value) + ", " + self.right.BuildNewick() + ":" + str(self.right_value) + ")"

