### This README file describes the upgma.py and nj.py scripts
=========================================================

###### Purpose

  The purpose of these scripts is to use different distance methods for phylogenetic inference. A file containing pairwise distances is read in and the output is a Newick representation of the phylogenetic tree.

  **upgma.py** uses the UPGMA (Unweighted Pair Group Method with
  Arithmetic Mean) distance method.

  **nj.py** uses the Neighbor-Joining method.

  **tree.py** contains a tree class used in both **upgma.py** and **nj.py**.


###### Additional Features

Several distance files were used to test **upgma.py** and **nj.py**
* *sarich.txt* which was given to us
* *isosceles.txt* which was taken from the lecture notes
* *honeybee.txt* which was taken from V. Makarenkov


**How to Run the Code**

    python upgma.py distancefile

  or

    python nj.py distancefile

  NOTE: The distance file is assumed to be space delimited with 
        float values for the distances.

