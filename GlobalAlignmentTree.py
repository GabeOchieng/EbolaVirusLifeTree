import numpy as np
from Bio import SeqIO
import edlib

#       reading genomes from fasta files

for record in SeqIO.parse("fastaData/Zaire_genome.fasta", "fasta"):
    Zaire_genome = str(record.seq)
    print(str(record.seq).__len__())

for record in SeqIO.parse("fastaData/Sudan_genome.fasta", "fasta"):
    Sudan_genome = str(record.seq)
    print(str(record.seq).__len__())

for record in SeqIO.parse("fastaData/Reston_genome.fasta", "fasta"):
    Reston_genome = str(record.seq)
    print(str(record.seq).__len__())

for record in SeqIO.parse("fastaData/TaiForest_genome.fasta", "fasta"):
    TaiForest_genome = str(record.seq)
    print(str(record.seq).__len__())

for record in SeqIO.parse("fastaData/Bundibugyo_genome.fasta", "fasta"):
    Bundibugyo_genome = str(record.seq)
    print(str(record.seq).__len__())

#       finding the edit distance between all pairs of viruses using global alignment
#       edlib.align: doing the alignment by Myers bit vector algorithm

a = edlib.align(Zaire_genome, Sudan_genome)
zaire_sudan = a["editDistance"]
a = edlib.align(Zaire_genome, Reston_genome)
zaire_restone = a["editDistance"]
a = edlib.align(Zaire_genome, TaiForest_genome)
zaire_taiforest = a["editDistance"]
a = edlib.align(Zaire_genome, Bundibugyo_genome)
zaire_bundibugyo = a["editDistance"]

a = edlib.align(Sudan_genome, Reston_genome)
sudan_restone = a["editDistance"]
a = edlib.align(Sudan_genome, TaiForest_genome)
sudan_taiforest = a["editDistance"]
a = edlib.align(Sudan_genome, Bundibugyo_genome)
sudan_bundibugyo = a["editDistance"]

a = edlib.align(Reston_genome, TaiForest_genome)
reston_taiforest = a["editDistance"]
a = edlib.align(Reston_genome, Bundibugyo_genome)
reston_bundibugyo = a["editDistance"]

a = edlib.align(TaiForest_genome, Bundibugyo_genome)
taiForest_bundibugyo = a["editDistance"]

#   fill the distance matrix
matrix = np.asarray([[0, zaire_sudan, zaire_restone, zaire_taiforest, zaire_bundibugyo], # zaire
                      [zaire_sudan, 0, sudan_restone, sudan_taiforest, sudan_bundibugyo], # sudan
                      [zaire_restone, sudan_restone, 0, reston_taiforest, reston_bundibugyo], # reston
                      [zaire_taiforest, sudan_taiforest, reston_taiforest, 0, taiForest_bundibugyo], # tiforest
                      [zaire_bundibugyo, sudan_bundibugyo, reston_bundibugyo, taiForest_bundibugyo, 0]  # bundibugyo
                     ])
#   write the matrix to a csv file
np.savetxt("csv/Global.csv", matrix, delimiter=",")