import numpy as np
from Bio import SeqIO, Phylo

import edlib

#       reading genomes from fasta files
from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor

for record in SeqIO.parse("fastaData/sequence.fasta", "fasta"):
    Marburg_genome = str(record.seq)

for record in SeqIO.parse("fastaData/Zaire_genome.fasta", "fasta"):
    Zaire_genome = str(record.seq)

for record in SeqIO.parse("fastaData/Sudan_genome.fasta", "fasta"):
    Sudan_genome = str(record.seq)

for record in SeqIO.parse("fastaData/Reston_genome.fasta", "fasta"):
    Reston_genome = str(record.seq)

for record in SeqIO.parse("fastaData/TaiForest_genome.fasta", "fasta"):
    TaiForest_genome = str(record.seq)

for record in SeqIO.parse("fastaData/Bundibugyo_genome.fasta", "fasta"):
    Bundibugyo_genome = str(record.seq)

#       finding the P for all pairs of viruses using global alignment
#       edlib.align: doing the alignment by Myers bit vector algorithm

a = edlib.align(Zaire_genome, Sudan_genome)
zaire_sudan = a["editDistance"] / a["locations"][0][1]
a = edlib.align(Zaire_genome, Reston_genome)
zaire_restone = a["editDistance"] / a["locations"][0][1]
a = edlib.align(Zaire_genome, TaiForest_genome)
zaire_taiforest = a["editDistance"] / a["locations"][0][1]
a = edlib.align(Zaire_genome, Bundibugyo_genome)
zaire_bundibugyo = a["editDistance"] / a["locations"][0][1]
a = edlib.align(Zaire_genome, Marburg_genome)
zaire_marburg = a["editDistance"] / a["locations"][0][1]

a = edlib.align(Sudan_genome, Reston_genome)
sudan_restone = a["editDistance"] / a["locations"][0][1]
a = edlib.align(Sudan_genome, TaiForest_genome)
sudan_taiforest = a["editDistance"] / a["locations"][0][1]
a = edlib.align(Sudan_genome, Bundibugyo_genome)
sudan_bundibugyo = a["editDistance"] / a["locations"][0][1]
a = edlib.align(Sudan_genome, Marburg_genome)
sudan_marburg = a["editDistance"] / a["locations"][0][1]

a = edlib.align(Reston_genome, TaiForest_genome)
reston_taiforest = a["editDistance"] / a["locations"][0][1]
a = edlib.align(Reston_genome, Bundibugyo_genome)
reston_bundibugyo = a["editDistance"] / a["locations"][0][1]
a = edlib.align(Reston_genome, Marburg_genome)
reston_marburg = a["editDistance"] / a["locations"][0][1]

a = edlib.align(TaiForest_genome, Bundibugyo_genome)
taiForest_bundibugyo = a["editDistance"] / a["locations"][0][1]
a = edlib.align(TaiForest_genome, Marburg_genome)
taiForest_marburg = a["editDistance"] / a["locations"][0][1]

a = edlib.align(Bundibugyo_genome, Marburg_genome)
bundibugyo_marburg = a["editDistance"] / a["locations"][0][1]

#       filling a matrix with calculated Ps

matrix = np.asarray([[0, zaire_sudan, zaire_restone, zaire_taiforest, zaire_bundibugyo, zaire_marburg], # zaire
                     [zaire_sudan, 0, sudan_restone, sudan_taiforest, sudan_bundibugyo, sudan_marburg], # sudan
                     [zaire_restone, sudan_restone, 0, reston_taiforest, reston_bundibugyo, reston_marburg], # reston
                     [zaire_taiforest, sudan_taiforest, reston_taiforest, 0, taiForest_bundibugyo, taiForest_marburg], # tiforest
                     [zaire_bundibugyo, sudan_bundibugyo, reston_bundibugyo, taiForest_bundibugyo, 0, bundibugyo_marburg],  # bundibugyo
                     [zaire_marburg, sudan_marburg, reston_marburg, taiForest_marburg, bundibugyo_marburg, 0]
                     ])

#       applying the Jukes Contour Model's result equation to the values of the matrix
for i in range(6):
    for j in range(6):
        if i != j:
            matrix[i][j] = ((-3.0 / 4.0) * np.log(1 - (4.0 / 3.0) * matrix[i][j]) * 1000) / 1.9

#       printing the results
print("Time Distance between Marburg        and      Zair :       ", matrix[0][5])
print("Time Distance between Marburg        and      Sudan :      ", matrix[1][5])
print("Time Distance between Marburg        and      Reston :     ", matrix[2][5])
print("Time Distance between Marburg        and      TaiForest :  ", matrix[3][5])
print("Time Distance between Marburg        and      Bundibugyo : ", matrix[4][5])
print("Time Distance between Zair           and      Sudan :      ", matrix[0][1])
print("Time Distance between Zair           and      Reston :     ", matrix[0][2])
print("Time Distance between Zair           and      TaiForest :  ", matrix[0][3])
print("Time Distance between Zair           and      Bundibugyo : ", matrix[0][4])
print("Time Distance between Reston         and      Sudan :      ", matrix[1][2])
print("Time Distance between TaiForest      and      Sudan :      ", matrix[1][3])
print("Time Distance between Bundibugyo     and      Sudan :      ", matrix[1][4])
print("Time Distance between Reston         and      Taiforest :  ", matrix[2][3])
print("Time Distance between Reston         and      Bundibugyo : ", matrix[2][4])
print("Time Distance between Bundibugyo     and      Taiforest :  ", matrix[3][4])
print("\n**All times are in Years**")

names = ['Zaire', 'Sudan', 'Reston', 'TaiForest', 'Bundibugyo', 'Marburg']

n = 0
matrixUPGMA = []
#    filling the matrix of distances of NP gene
for row in matrix:
     sub = []
     for i in range(n):
          sub.append(float(row[i]))
     sub.append(0.0)
     matrixUPGMA.append(sub)
     n += 1     # the matrix has to be lower triangular
dm = DistanceMatrix(names=names, matrix=matrixUPGMA)
constructor = DistanceTreeConstructor()

tree = constructor.upgma(dm)
Phylo.draw(tree)

print(tree.distance("Inner5", "Marburg"), " years")

