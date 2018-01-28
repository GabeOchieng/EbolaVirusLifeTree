import edlib
import numpy as np
from Bio import pairwise2

text_file = open("genes/zairNP.txt", "r")
zaire1 = text_file.read()

text_file = open("genes/sudanNP.txt", "r")
sudan1 = text_file.read()

text_file = open("genes/restonNP.txt", "r")
reston1 = text_file.read()

text_file = open("genes/taiForestNP.txt", "r")
taiForest1 = text_file.read()

text_file = open("genes/bundibugyoNP.txt", "r")
bundibugyo1 = text_file.read()

a = edlib.align(zaire1, sudan1)
zaire_sudan1 = a["editDistance"]
a = edlib.align(zaire1, reston1)
zaire_restone1 = a["editDistance"]
a = edlib.align(zaire1, taiForest1)
zaire_taiforest1 = a["editDistance"]
a = edlib.align(zaire1, bundibugyo1)
zaire_bundibugyo1 = a["editDistance"]

a = edlib.align(sudan1, reston1)
sudan_restone1 = a["editDistance"]
a = edlib.align(sudan1, taiForest1)
sudan_taiforest1 = a["editDistance"]
a = edlib.align(sudan1, bundibugyo1)
sudan_bundibugyo1 = a["editDistance"]

a = edlib.align(reston1, taiForest1)
reston_taiforest1 = a["editDistance"]
a = edlib.align(reston1, bundibugyo1)
reston_bundibugyo1 = a["editDistance"]

a = edlib.align(taiForest1, bundibugyo1)
taiForest_bundibugyo1 = a["editDistance"]



text_file = open("genes/zairVP35.txt", "r")
zaire2 = text_file.read()

text_file = open("genes/sudanVP35.txt", "r")
sudan2 = text_file.read()

text_file = open("genes/restonVP35.txt", "r")
reston2 = text_file.read()

text_file = open("genes/taiForestVP35.txt", "r")
taiForest2 = text_file.read()

text_file = open("genes/bundibugyoVP35.txt", "r")
bundibugyo2 = text_file.read()

a = edlib.align(zaire2, sudan2)
zaire_sudan2 = a["editDistance"]
a = edlib.align(zaire2, reston2)
zaire_restone2 = a["editDistance"]
a = edlib.align(zaire2, taiForest2)
zaire_taiforest2 = a["editDistance"]
a = edlib.align(zaire2, bundibugyo2)
zaire_bundibugyo2 = a["editDistance"]

a = edlib.align(sudan2, reston2)
sudan_restone2 = a["editDistance"]
a = edlib.align(sudan2, taiForest2)
sudan_taiforest2 = a["editDistance"]
a = edlib.align(sudan2, bundibugyo2)
sudan_bundibugyo2 = a["editDistance"]

a = edlib.align(reston2, taiForest2)
reston_taiforest2 = a["editDistance"]
a = edlib.align(reston2, bundibugyo2)
reston_bundibugyo2 = a["editDistance"]

a = edlib.align(taiForest2, bundibugyo2)
taiForest_bundibugyo2 = a["editDistance"]



text_file = open("genes/zairVP40.txt", "r")
zaire3 = text_file.read()

text_file = open("genes/sudanVP40.txt", "r")
sudan3 = text_file.read()

text_file = open("genes/restonVP40.txt", "r")
reston3 = text_file.read()

text_file = open("genes/taiForestVP40.txt", "r")
taiForest3 = text_file.read()

text_file = open("genes/bundibugyoVP40.txt", "r")
bundibugyo3 = text_file.read()

a = edlib.align(zaire3, sudan3)
zaire_sudan3 = a["editDistance"]
a = edlib.align(zaire3, reston3)
zaire_restone3 = a["editDistance"]
a = edlib.align(zaire3, taiForest3)
zaire_taiforest3 = a["editDistance"]
a = edlib.align(zaire3, bundibugyo3)
zaire_bundibugyo3 = a["editDistance"]

a = edlib.align(sudan3, reston3)
sudan_restone3 = a["editDistance"]
a = edlib.align(sudan3, taiForest3)
sudan_taiforest3 = a["editDistance"]
a = edlib.align(sudan3, bundibugyo3)
sudan_bundibugyo3 = a["editDistance"]

a = edlib.align(reston3, taiForest3)
reston_taiforest3 = a["editDistance"]
a = edlib.align(reston3, bundibugyo3)
reston_bundibugyo3 = a["editDistance"]

a = edlib.align(taiForest3, bundibugyo3)
taiForest_bundibugyo3 = a["editDistance"]


text_file = open("genes/zairGP.txt", "r")
zaire4 = text_file.read()

text_file = open("genes/sudanGP.txt", "r")
sudan4 = text_file.read()

text_file = open("genes/restonGP.txt", "r")
reston4 = text_file.read()

text_file = open("genes/taiForestGP.txt", "r")
taiForest4 = text_file.read()

text_file = open("genes/bundibugyoGP.txt", "r")
bundibugyo4 = text_file.read()

a = edlib.align(zaire4, sudan4)
zaire_sudan4 = a["editDistance"]
a = edlib.align(zaire4, reston4)
zaire_restone4 = a["editDistance"]
a = edlib.align(zaire4, taiForest4)
zaire_taiforest4 = a["editDistance"]
a = edlib.align(zaire4, bundibugyo4)
zaire_bundibugyo4 = a["editDistance"]

a = edlib.align(sudan4, reston4)
sudan_restone4 = a["editDistance"]
a = edlib.align(sudan4, taiForest4)
sudan_taiforest4 = a["editDistance"]
a = edlib.align(sudan4, bundibugyo4)
sudan_bundibugyo4 = a["editDistance"]

a = edlib.align(reston4, taiForest4)
reston_taiforest4 = a["editDistance"]
a = edlib.align(reston4, bundibugyo4)
reston_bundibugyo4 = a["editDistance"]

a = edlib.align(taiForest4, bundibugyo4)
taiForest_bundibugyo4 = a["editDistance"]


text_file = open("genes/zairVP30.txt", "r")
zaire5 = text_file.read()

text_file = open("genes/sudanVP30.txt", "r")
sudan5 = text_file.read()

text_file = open("genes/restonVP30.txt", "r")
reston5 = text_file.read()

text_file = open("genes/taiForestVP30.txt", "r")
taiForest5 = text_file.read()

text_file = open("genes/bundibugyoVP30.txt", "r")
bundibugyo5 = text_file.read()

a = edlib.align(zaire5, sudan5)
zaire_sudan5 = a["editDistance"]
a = edlib.align(zaire5, reston5)
zaire_restone5 = a["editDistance"]
a = edlib.align(zaire5, taiForest5)
zaire_taiforest5 = a["editDistance"]
a = edlib.align(zaire5, bundibugyo5)
zaire_bundibugyo5 = a["editDistance"]

a = edlib.align(sudan5, reston5)
sudan_restone5 = a["editDistance"]
a = edlib.align(sudan5, taiForest5)
sudan_taiforest5 = a["editDistance"]
a = edlib.align(sudan5, bundibugyo5)
sudan_bundibugyo5 = a["editDistance"]

a = edlib.align(reston5, taiForest5)
reston_taiforest5 = a["editDistance"]
a = edlib.align(reston5, bundibugyo5)
reston_bundibugyo5 = a["editDistance"]

a = edlib.align(taiForest5, bundibugyo5)
taiForest_bundibugyo5 = a["editDistance"]


text_file = open("genes/zairVP24.txt", "r")
zaire6 = text_file.read()

text_file = open("genes/sudanVP24.txt", "r")
sudan6 = text_file.read()

text_file = open("genes/restonVP24.txt", "r")
reston6 = text_file.read()

text_file = open("genes/taiForestVP24.txt", "r")
taiForest6 = text_file.read()

text_file = open("genes/bundibugyoVP24.txt", "r")
bundibugyo6 = text_file.read()

a = edlib.align(zaire6, sudan6)
zaire_sudan6 = a["editDistance"]
a = edlib.align(zaire6, reston6)
zaire_restone6 = a["editDistance"]
a = edlib.align(zaire6, taiForest6)
zaire_taiforest6 = a["editDistance"]
a = edlib.align(zaire6, bundibugyo6)
zaire_bundibugyo6 = a["editDistance"]

a = edlib.align(sudan6, reston6)
sudan_restone6 = a["editDistance"]
a = edlib.align(sudan6, taiForest6)
sudan_taiforest6 = a["editDistance"]
a = edlib.align(sudan6, bundibugyo6)
sudan_bundibugyo6 = a["editDistance"]

a = edlib.align(reston6, taiForest6)
reston_taiforest6 = a["editDistance"]
a = edlib.align(reston6, bundibugyo6)
reston_bundibugyo6 = a["editDistance"]

a = edlib.align(taiForest6, bundibugyo6)
taiForest_bundibugyo6 = a["editDistance"]


text_file = open("genes/zairL.txt", "r")
zaire7 = text_file.read()

text_file = open("genes/sudanL.txt", "r")
sudan7 = text_file.read()

text_file = open("genes/restonL.txt", "r")
reston7 = text_file.read()

text_file = open("genes/taiForestL.txt", "r")
taiForest7 = text_file.read()

text_file = open("genes/bundibugyoL.txt", "r")
bundibugyo7 = text_file.read()

a = edlib.align(zaire7, sudan7)
zaire_sudan7 = a["editDistance"]
a = edlib.align(zaire7, reston7)
zaire_restone7 = a["editDistance"]
a = edlib.align(zaire7, taiForest7)
zaire_taiforest7 = a["editDistance"]
a = edlib.align(zaire7, bundibugyo7)
zaire_bundibugyo7 = a["editDistance"]

a = edlib.align(sudan7, reston7)
sudan_restone7 = a["editDistance"]
a = edlib.align(sudan7, taiForest7)
sudan_taiforest7 = a["editDistance"]
a = edlib.align(sudan7, bundibugyo7)
sudan_bundibugyo7 = a["editDistance"]

a = edlib.align(reston7, taiForest7)
reston_taiforest7 = a["editDistance"]
a = edlib.align(reston7, bundibugyo7)
reston_bundibugyo7 = a["editDistance"]

a = edlib.align(taiForest7, bundibugyo7)
taiForest_bundibugyo7 = a["editDistance"]


matrix1 = np.asarray([[0, zaire_sudan1, zaire_restone1, zaire_taiforest1, zaire_bundibugyo1], # zaire
                      [zaire_sudan1, 0, sudan_restone1, sudan_taiforest1, sudan_bundibugyo1], # sudan
                      [zaire_restone1, sudan_restone1, 0, reston_taiforest1, reston_bundibugyo1], # reston
                      [zaire_taiforest1, sudan_taiforest1, reston_taiforest1, 0, taiForest_bundibugyo1], # tiforest
                      [zaire_bundibugyo1, sudan_bundibugyo1, reston_bundibugyo1, taiForest_bundibugyo1, 0]  # bundibugyo
                     ])

np.savetxt("csv/NP1.csv", matrix1, delimiter=",")


matrix2 = np.asarray([[0, zaire_sudan2, zaire_restone2, zaire_taiforest2, zaire_bundibugyo2], # zaire
                      [zaire_sudan2, 0, sudan_restone2, sudan_taiforest2, sudan_bundibugyo2], # sudan
                      [zaire_restone2, sudan_restone2, 0, reston_taiforest2, reston_bundibugyo2], # reston
                      [zaire_taiforest2, sudan_taiforest2, reston_taiforest2, 0, taiForest_bundibugyo2], # tiforest
                      [zaire_bundibugyo2, sudan_bundibugyo2, reston_bundibugyo2, taiForest_bundibugyo2, 0]  # bundibugyo
                     ])

np.savetxt("csv/VP351.csv", matrix2, delimiter=",")


matrix3 = np.asarray([[0, zaire_sudan3, zaire_restone3, zaire_taiforest3, zaire_bundibugyo3], # zaire
                      [zaire_sudan3, 0, sudan_restone3, sudan_taiforest3, sudan_bundibugyo3], # sudan
                      [zaire_restone3, sudan_restone3, 0, reston_taiforest3, reston_bundibugyo3], # reston
                      [zaire_taiforest3, sudan_taiforest3, reston_taiforest3, 0, taiForest_bundibugyo3], # tiforest
                      [zaire_bundibugyo3, sudan_bundibugyo3, reston_bundibugyo3, taiForest_bundibugyo3, 0]  # bundibugyo
                     ])

np.savetxt("csv/VP401.csv", matrix3, delimiter=",")


matrix4 = np.asarray([[0, zaire_sudan4, zaire_restone4, zaire_taiforest4, zaire_bundibugyo4], # zaire
                      [zaire_sudan4, 0, sudan_restone4, sudan_taiforest4, sudan_bundibugyo4], # sudan
                      [zaire_restone4, sudan_restone4, 0, reston_taiforest4, reston_bundibugyo4], # reston
                      [zaire_taiforest4, sudan_taiforest4, reston_taiforest4, 0, taiForest_bundibugyo4], # tiforest
                      [zaire_bundibugyo4, sudan_bundibugyo4, reston_bundibugyo4, taiForest_bundibugyo4, 0]  # bundibugyo
                     ])

np.savetxt("csv/GP1.csv", matrix4, delimiter=",")


matrix5 = np.asarray([[0, zaire_sudan5, zaire_restone5, zaire_taiforest5, zaire_bundibugyo5], # zaire
                      [zaire_sudan5, 0, sudan_restone5, sudan_taiforest5, sudan_bundibugyo5], # sudan
                      [zaire_restone5, sudan_restone5, 0, reston_taiforest5, reston_bundibugyo5], # reston
                      [zaire_taiforest5, sudan_taiforest5, reston_taiforest5, 0, taiForest_bundibugyo5], # tiforest
                      [zaire_bundibugyo5, sudan_bundibugyo5, reston_bundibugyo5, taiForest_bundibugyo5, 0]  # bundibugyo
                     ])

np.savetxt("csv/VP301.csv", matrix5, delimiter=",")


matrix6 = np.asarray([[0, zaire_sudan6, zaire_restone6, zaire_taiforest6, zaire_bundibugyo6], # zaire
                      [zaire_sudan6, 0, sudan_restone6, sudan_taiforest6, sudan_bundibugyo6], # sudan
                      [zaire_restone6, sudan_restone6, 0, reston_taiforest6, reston_bundibugyo6], # reston
                      [zaire_taiforest6, sudan_taiforest6, reston_taiforest6, 0, taiForest_bundibugyo6], # tiforest
                      [zaire_bundibugyo6, sudan_bundibugyo6, reston_bundibugyo6, taiForest_bundibugyo6, 0]  # bundibugyo
                     ])

np.savetxt("csv/VP241.csv", matrix6, delimiter=",")


matrix7 = np.asarray([[0, zaire_sudan7, zaire_restone7, zaire_taiforest7, zaire_bundibugyo7], # zaire
                      [zaire_sudan7, 0, sudan_restone7, sudan_taiforest7, sudan_bundibugyo7], # sudan
                      [zaire_restone7, sudan_restone7, 0, reston_taiforest7, reston_bundibugyo7], # reston
                      [zaire_taiforest7, sudan_taiforest7, reston_taiforest7, 0, taiForest_bundibugyo7], # tiforest
                      [zaire_bundibugyo7, sudan_bundibugyo7, reston_bundibugyo7, taiForest_bundibugyo7, 0]  # bundibugyo
                     ])

np.savetxt("csv/L1.csv", matrix7, delimiter=",")


import csv

from Bio import Phylo
from Bio.Phylo.Consensus import adam_consensus, strict_consensus, majority_consensus
from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor

constructor = DistanceTreeConstructor()

#reading the csv files of Genes

NP = open('csv/NP1.csv', "rt", encoding="utf8")
VP35 = open('csv/VP351.csv', "rt", encoding="utf8")
VP40 = open('csv/VP401.csv', "rt", encoding="utf8")
GP = open('csv/GP1.csv', "rt", encoding="utf8")
VP30 = open('csv/VP301.csv', "rt", encoding="utf8")
VP24 = open('csv/VP241.csv', "rt", encoding="utf8")
L = open('csv/L1.csv', "rt", encoding="utf8")

#names to be displayed in on the trees

names = ['Zaire', 'Sudan', 'Reston', 'TaiForest', 'Bundibugyo']

read = csv.reader(NP)
n = 0
matrixNP = []
#filling the matrix of deistances of NP gene
for row in read:
     sub = []
     for i in range(n):
          sub.append(float(row[i]))
     sub.append(0.0)
     matrixNP.append(sub)
     n += 1     #the matrix has to be lower triangular
dm = DistanceMatrix(names=names, matrix=matrixNP)

NJtreeNP = constructor.nj(dm)
UPGMAtreeNP = constructor.upgma(dm)

NJtreeNP.ladderize()   # Flip branches so deeper clades are displayed at top
Phylo.draw(NJtreeNP)

UPGMAtreeNP.ladderize()   # Flip branches so deeper clades are displayed at top
Phylo.draw(UPGMAtreeNP)



read = csv.reader(VP35)
n = 0
matrixVP35 = []
#filling the matrix of deistances of VP35 gene
for row in read:
     sub = []
     for i in range(n):
          sub.append(float(row[i]))
     sub.append(0.0)
     matrixVP35.append(sub)
     n += 1
dmVP35 = DistanceMatrix(names=names, matrix=matrixVP35)

NJtreeVP35 = constructor.nj(dmVP35)
UPGMAtreeVP35 = constructor.upgma(dmVP35)

NJtreeVP35.ladderize()   # Flip branches so deeper clades are displayed at top
Phylo.draw(NJtreeVP35)

UPGMAtreeVP35.ladderize()   # Flip branches so deeper clades are displayed at top
Phylo.draw(UPGMAtreeVP35)



read = csv.reader(VP40)
n = 0
matrixVP40 = []
#filling the matrix of deistances of VP40 gene
for row in read:
     sub = []
     for i in range(n):
          sub.append(float(row[i]))
     sub.append(0.0)
     matrixVP40.append(sub)
     n += 1
dmVP40 = DistanceMatrix(names=names, matrix=matrixVP40)

NJtreeVP40 = constructor.nj(dmVP40)
UPGMAtreeVP40 = constructor.upgma(dmVP40)

NJtreeVP40.ladderize()   # Flip branches so deeper clades are displayed at top
Phylo.draw(NJtreeVP40)

UPGMAtreeVP40.ladderize()   # Flip branches so deeper clades are displayed at top
Phylo.draw(UPGMAtreeVP40)



read = csv.reader(GP)
n = 0
matrixGP = []
#filling the matrix of deistances of GP gene
for row in read:
     sub = []
     for i in range(n):
          sub.append(float(row[i]))
     sub.append(0.0)
     matrixGP.append(sub)
     n += 1
dmGP = DistanceMatrix(names=names, matrix=matrixGP)

NJtreeGP = constructor.nj(dmGP)
UPGMAtreeGP = constructor.upgma(dmGP)

NJtreeGP.ladderize()   # Flip branches so deeper clades are displayed at top
Phylo.draw(NJtreeGP)

UPGMAtreeGP.ladderize()   # Flip branches so deeper clades are displayed at top
Phylo.draw(UPGMAtreeGP)



read = csv.reader(VP30)
n = 0
matrixVP30 = []
#filling the matrix of deistances of VP30 gene
for row in read:
     sub = []
     for i in range(n):
          sub.append(float(row[i]))
     sub.append(0.0)
     matrixVP30.append(sub)
     n += 1
dmVP30 = DistanceMatrix(names=names, matrix=matrixVP30)

NJtreeVP30 = constructor.nj(dmVP30)
UPGMAtreeVP30 = constructor.upgma(dmVP30)

NJtreeVP30.ladderize()   # Flip branches so deeper clades are displayed at top
Phylo.draw(NJtreeVP30)

UPGMAtreeVP30.ladderize()   # Flip branches so deeper clades are displayed at top
Phylo.draw(UPGMAtreeVP30)



read = csv.reader(VP24)
n = 0
matrixVP24 = []
#filling the matrix of deistances of VP24 gene
for row in read:
     sub = []
     for i in range(n):
          sub.append(float(row[i]))
     sub.append(0.0)
     matrixVP24.append(sub)
     n += 1
dmVP24 = DistanceMatrix(names=names, matrix=matrixVP24)

NJtreeVP24 = constructor.nj(dmVP24)
UPGMAtreeVP24 = constructor.upgma(dmVP24)

NJtreeVP24.ladderize()   # Flip branches so deeper clades are displayed at top
Phylo.draw(NJtreeVP24)

UPGMAtreeVP24.ladderize()   # Flip branches so deeper clades are displayed at top
Phylo.draw(UPGMAtreeVP24)



read = csv.reader(L)
n = 0
matrixL = []
#filling the matrix of deistances of L gene
for row in read:
     sub = []
     for i in range(n):
          sub.append(float(row[i]))
     sub.append(0.0)
     matrixL.append(sub)
     n += 1
dmL = DistanceMatrix(names=names, matrix=matrixL)

NJtreeL = constructor.nj(dmL)
UPGMAtreeL = constructor.upgma(dmL)

NJtreeL.ladderize()   # Flip branches so deeper clades are displayed at top
Phylo.draw(NJtreeL)

UPGMAtreeL.ladderize()   # Flip branches so deeper clades are displayed at top
Phylo.draw(UPGMAtreeL)

trees = [NJtreeVP35, NJtreeGP, NJtreeL, NJtreeNP, NJtreeVP24, NJtreeVP30, NJtreeVP40]
treesUPGMA = [UPGMAtreeVP35, UPGMAtreeGP, UPGMAtreeL, UPGMAtreeNP, UPGMAtreeVP24, UPGMAtreeVP30, UPGMAtreeVP40]
consensus = majority_consensus(treesUPGMA, cutoff=0)
consensusNJ = majority_consensus(trees, cutoff=0)
Phylo.draw(consensus)
Phylo.draw(consensusNJ)