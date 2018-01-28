import csv

from Bio import Phylo
from Bio.Phylo.Consensus import adam_consensus, strict_consensus, majority_consensus
from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor

constructor = DistanceTreeConstructor()

# reading the csv file of distance matrix

Global = open('csv/Global.csv', "rt", encoding="utf8")

# names to be displayed in on the tree

names = ['Zaire', 'Sudan', 'Reston', 'TaiForest', 'Bundibugyo']

read = csv.reader(Global)
n = 0
matrix = []
# filling the matrix of distances
for row in read:
     sub = []
     for i in range(n):
          sub.append(float(row[i]))
     sub.append(0.0)
     matrix.append(sub)
     n += 1     # the matrix has to be lower triangular
dm = DistanceMatrix(names=names, matrix=matrix)

NJtreeNP = constructor.nj(dm)

NJtreeNP.ladderize()   # Flip branches so deeper clades are displayed at top
Phylo.draw(NJtreeNP)
