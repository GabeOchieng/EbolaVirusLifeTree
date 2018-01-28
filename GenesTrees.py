import csv

from Bio import Phylo
from Bio.Phylo.Consensus import adam_consensus, strict_consensus, majority_consensus
from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor
from Bio import Phylo


constructor = DistanceTreeConstructor()

#    reading the csv files of Genes

NP = open('csv/NP.csv', "rt", encoding="utf8")
VP35 = open('csv/VP35.csv', "rt", encoding="utf8")
VP40 = open('csv/VP40.csv', "rt", encoding="utf8")
GP = open('csv/GP.csv', "rt", encoding="utf8")
VP30 = open('csv/VP30.csv', "rt", encoding="utf8")
VP24 = open('csv/VP24.csv', "rt", encoding="utf8")
L = open('csv/L.csv', "rt", encoding="utf8")

#    names to be displayed in on the trees

names = ['Zaire', 'Sudan', 'Reston', 'TaiForest', 'Bundibugyo']

read = csv.reader(NP)
n = 0
matrixNP = []
#    filling the matrix of distances of NP gene
for row in read:
     sub = []
     for i in range(n):
          sub.append(int(float(row[i])))
     sub.append(0.0)
     matrixNP.append(sub)
     n += 1     # the matrix has to be lower triangular
dm = DistanceMatrix(names=names, matrix=matrixNP)
print(dm)
NJtreeNP = constructor.nj(dm)
UPGMAtreeNP = constructor.upgma(dm)

Phylo.draw(NJtreeNP)

Phylo.draw(UPGMAtreeNP)



read = csv.reader(VP35)
n = 0
matrixVP35 = []
#    filling the matrix of distances of VP35 gene
for row in read:
     sub = []
     for i in range(n):
          sub.append(int(float(row[i])))
     sub.append(0.0)
     matrixVP35.append(sub)
     n += 1
dmVP35 = DistanceMatrix(names=names, matrix=matrixVP35)

NJtreeVP35 = constructor.nj(dmVP35)
UPGMAtreeVP35 = constructor.upgma(dmVP35)

Phylo.draw(NJtreeVP35)

Phylo.draw(UPGMAtreeVP35)



read = csv.reader(VP40)
n = 0
matrixVP40 = []
# filling the matrix of distances of VP40 gene
for row in read:
     sub = []
     for i in range(n):
          sub.append(int(float(row[i])))
     sub.append(0.0)
     matrixVP40.append(sub)
     n += 1
dmVP40 = DistanceMatrix(names=names, matrix=matrixVP40)

NJtreeVP40 = constructor.nj(dmVP40)
UPGMAtreeVP40 = constructor.upgma(dmVP40)

Phylo.draw(NJtreeVP40)

Phylo.draw(UPGMAtreeVP40)



read = csv.reader(GP)
n = 0
matrixGP = []
# filling the matrix of distances of GP gene
for row in read:
     sub = []
     for i in range(n):
          sub.append(int(float(row[i])))
     sub.append(0.0)
     matrixGP.append(sub)
     n += 1
dmGP = DistanceMatrix(names=names, matrix=matrixGP)

NJtreeGP = constructor.nj(dmGP)
UPGMAtreeGP = constructor.upgma(dmGP)

Phylo.draw(NJtreeGP)

Phylo.draw(UPGMAtreeGP)



read = csv.reader(VP30)
n = 0
matrixVP30 = []
# filling the matrix of distances of VP30 gene
for row in read:
     sub = []
     for i in range(n):
          sub.append(int(float(row[i])))
     sub.append(0.0)
     matrixVP30.append(sub)
     n += 1
dmVP30 = DistanceMatrix(names=names, matrix=matrixVP30)

NJtreeVP30 = constructor.nj(dmVP30)
UPGMAtreeVP30 = constructor.upgma(dmVP30)

Phylo.draw(NJtreeVP30)

Phylo.draw(UPGMAtreeVP30)



read = csv.reader(VP24)
n = 0
matrixVP24 = []
# filling the matrix of distances of VP24 gene
for row in read:
     sub = []
     for i in range(n):
          sub.append(int(float(row[i])))
     sub.append(0.0)
     matrixVP24.append(sub)
     n += 1
dmVP24 = DistanceMatrix(names=names, matrix=matrixVP24)

NJtreeVP24 = constructor.nj(dmVP24)
UPGMAtreeVP24 = constructor.upgma(dmVP24)

Phylo.draw(NJtreeVP24)

Phylo.draw(UPGMAtreeVP24)



read = csv.reader(L)
n = 0
matrixL = []
# filling the matrix of distances of L gene
for row in read:
     sub = []
     for i in range(n):
          sub.append(int(float(row[i])))
     sub.append(0.0)
     matrixL.append(sub)
     n += 1
dmL = DistanceMatrix(names=names, matrix=matrixL)

NJtreeL = constructor.nj(dmL)
UPGMAtreeL = constructor.upgma(dmL)

Phylo.draw(NJtreeL)

Phylo.draw(UPGMAtreeL)

#    list the trees built by NJ algorithm
treesNJ = [NJtreeVP35, NJtreeGP, NJtreeL, NJtreeNP, NJtreeVP24, NJtreeVP30, NJtreeVP40]
#    list the trees built by UPGMA algorithm
treesUPGMA = [UPGMAtreeVP35, UPGMAtreeGP, UPGMAtreeL, UPGMAtreeNP, UPGMAtreeVP24, UPGMAtreeVP30, UPGMAtreeVP40]

#    build the consensus trees
consensusUPGMA = majority_consensus(treesUPGMA, cutoff=0.5)
consensusNJ = majority_consensus(treesNJ, cutoff=0.5)

Phylo.draw(consensusNJ)
Phylo.draw(consensusUPGMA)

