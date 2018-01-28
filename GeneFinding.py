from Bio import pairwise2, SeqIO

Marburg_Genes = []
i = 0

#       reading Marburg genes into a list : Marburg_genes
for record in SeqIO.parse("fastaData/Marburg_Genes.fasta", "fasta"):
    Marburg_Genes.append(str(record.seq))

#       reading Zaire genome data from the fasta files
for record in SeqIO.parse("fastaData/Zaire_genome.fasta", "fasta"):
    Zaire_genome = str(record.seq)

#       reading Sudan genome data from the fasta files
for record in SeqIO.parse("fastaData/Sudan_genome.fasta", "fasta"):
    Sudan_genome = str(record.seq)

#       reading Reston genome data from the fasta files
for record in SeqIO.parse("fastaData/Reston_genome.fasta", "fasta"):
    Reston_genome = str(record.seq)

#       reading TaiForest genome data from the fasta files
for record in SeqIO.parse("fastaData/TaiForest_genome.fasta", "fasta"):
    TaiForest_genome = str(record.seq)

#       reading Bundibugyo genome data from the fasta files
for record in SeqIO.parse("fastaData/Bundibugyo_genome.fasta", "fasta"):
    Bundibugyo_genome = str(record.seq)


#           finding the first gene in the genomes : NP

#       local alignment to find the gene, interval: 0-4000 (about 2 times of gene's length)
alignments = pairwise2.align.localms(Zaire_genome[0:4000], Marburg_Genes[0], 1, -1, -1, -1)
start = alignments[0][3]   # the start point of genes location
end = 4000 - len(alignments[0][0][alignments[0][4]:]) # the end point
zaire1 = Zaire_genome[start:end]  # substring of genome that the gene is located in it
with open("genes/zairNP.txt", "w") as text_file:
    text_file.write(zaire1)    # writing the gene to a file
# printing the start and end points to be used for next interval specifications
print(str(alignments[0][3]) + "  " + str(4000 - len(alignments[0][0][alignments[0][4]:])))

#       same for Sudan genome
alignments = pairwise2.align.localms(Sudan_genome[0:4000], Marburg_Genes[0], 1, -1, -1, -1)
start = alignments[0][3]
end = 4000 - len(alignments[0][0][alignments[0][4]:])
sudan1 = Sudan_genome[start:end]
with open("genes/sudanNP.txt", "w") as text_file:
    text_file.write(sudan1)
print(str(alignments[0][3]) + "  " + str(4000 - len(alignments[0][0][alignments[0][4]:])))

#       same for Reston genome
alignments = pairwise2.align.localms(Reston_genome[0:4000], Marburg_Genes[0], 1, -1, -1, -1)
start = alignments[0][3]
end = 4000 - len(alignments[0][0][alignments[0][4]:])
reston1 = Reston_genome[start:end]
with open("genes/restonNP.txt", "w") as text_file:
    text_file.write(reston1)
print(str(alignments[0][3]) + "  " + str(4000 - len(alignments[0][0][alignments[0][4]:])))

#       same for TaiForest genome
alignments = pairwise2.align.localms(TaiForest_genome[0:4000], Marburg_Genes[0], 1, -1, -1, -1)
start = alignments[0][3]
end = 4000 - len(alignments[0][0][alignments[0][4]:])
taiForest1 = TaiForest_genome[start:end]
with open("genes/taiForestNP.txt", "w") as text_file:
    text_file.write(taiForest1)
print(str(alignments[0][3]) + "  " + str(4000 - len(alignments[0][0][alignments[0][4]:])))

#       same for Bundibugyo genome
alignments = pairwise2.align.localms(Bundibugyo_genome[0:4000], Marburg_Genes[0], 1, -1, -1, -1)
start = alignments[0][3]
end = 4000 - len(alignments[0][0][alignments[0][4]:])
bundibugyo1 = Bundibugyo_genome[start:end]
with open("genes/bundibugyoNP.txt", "w") as text_file:
    text_file.write(bundibugyo1)
print(str(alignments[0][3]) + "  " + str(4000 - len(alignments[0][0][alignments[0][4]:])))


#           finding the second gene in the genomes : VP35

#       local alignment to find the gene, interval: 3000-6000 (3000: a little bit less than end
#                                                              points found in previuse gene)
alignments = pairwise2.align.localms(Zaire_genome[3000:6000], Marburg_Genes[1], 1, -1, -1, -1)
start = alignments[0][3] + 3000  # the start point of genes location
end = 6000 - len(alignments[0][0][alignments[0][4]:])  # the end point
zaire2 = Zaire_genome[start:end]
with open("genes/zairVP35.txt", "w") as text_file:
    text_file.write(zaire2)    # writing the gene to a file
# printing the start and end points to be used for next interval specifications
print(str(alignments[0][3] + 3000) + "  " + str(6000 - len(alignments[0][0][alignments[0][4]:])))

#       same for Sudan genome
alignments = pairwise2.align.localms(Sudan_genome[3000:6000], Marburg_Genes[1], 1, -1, -1, -1)
start = alignments[0][3] + 3000
end = 6000 - len(alignments[0][0][alignments[0][4]:])
sudan2 = Sudan_genome[start:end]
with open("genes/sudanVP35.txt", "w") as text_file:
    text_file.write(sudan2)
print(str(alignments[0][3] + 3000) + "  " + str(6000 - len(alignments[0][0][alignments[0][4]:])))

#       same for Reston genome
alignments = pairwise2.align.localms(Reston_genome[3000:6000], Marburg_Genes[1], 1, -1, -1, -1)
start = alignments[0][3] + 3000
end = 6000 - len(alignments[0][0][alignments[0][4]:])
reston2 = Reston_genome[start:end]
with open("genes/restonVP35.txt", "w") as text_file:
    text_file.write(reston2)
print(str(alignments[0][3] + 3000) + "  " + str(6000 - len(alignments[0][0][alignments[0][4]:])))

#       same for TaiForest genome
alignments = pairwise2.align.localms(TaiForest_genome[3000:6000], Marburg_Genes[1], 1, -1, -1, -1)
start = alignments[0][3] + 3000
end = 6000 - len(alignments[0][0][alignments[0][4]:])
taiForest2 = TaiForest_genome[start:end]
with open("genes/taiForestVP35.txt", "w") as text_file:
    text_file.write(taiForest2)
print(str(alignments[0][3] + 3000) + "  " + str(6000 - len(alignments[0][0][alignments[0][4]:])))

#       same for Bundibugyo genome
alignments = pairwise2.align.localms(Bundibugyo_genome[3000:6000], Marburg_Genes[1], 1, -1, -1, -1)
start = alignments[0][3] + 3000
end = 6000 - len(alignments[0][0][alignments[0][4]:])
bundibugyo2 = Bundibugyo_genome[start:end]
with open("genes/bundibugyoVP35.txt", "w") as text_file:
    text_file.write(bundibugyo2)
print(str(alignments[0][3] + 3000) + "  " + str(6000 - len(alignments[0][0][alignments[0][4]:])))


#           finding the third gene in the genomes : VP40

#       local alignment to find the gene, interval: 4100-8000 (4100: a little bit less than end
#                                                              points found in previuse gene)

alignments = pairwise2.align.localms(Zaire_genome[4100:8000], Marburg_Genes[2], 1, -1, -1, -1)
start = alignments[0][3] + 4100
end = 8000 - len(alignments[0][0][alignments[0][4]:])
zaire3 = Zaire_genome[start:end]
with open("genes/zairVP40.txt", "w") as text_file:
    text_file.write(zaire3)
print(str(alignments[0][3] + 4100) + "  " + str(8000 - len(alignments[0][0][alignments[0][4]:])))

#       same for Sudan genome
alignments = pairwise2.align.localms(Sudan_genome[4100:8000], Marburg_Genes[2], 1, -1, -1, -1)
start = alignments[0][3] + 4100
end = 8000 - len(alignments[0][0][alignments[0][4]:])
sudan3 = Sudan_genome[start:end]
with open("genes/sudanVP40.txt", "w") as text_file:
    text_file.write(sudan3)
print(str(alignments[0][3] + 4100) + "  " + str(8000 - len(alignments[0][0][alignments[0][4]:])))

#       same for Reston genome
alignments = pairwise2.align.localms(Reston_genome[4100:8000], Marburg_Genes[2], 1, -1, -1, -1)
start = alignments[0][3] + 4100
end = 8000 - len(alignments[0][0][alignments[0][4]:])
reston3 = Reston_genome[start:end]
with open("genes/restonVP40.txt", "w") as text_file:
    text_file.write(reston3)
print(str(alignments[0][3] + 4100) + "  " + str(8000 - len(alignments[0][0][alignments[0][4]:])))

#       same for TaiForest genome
alignments = pairwise2.align.localms(TaiForest_genome[4100:8000], Marburg_Genes[2], 1, -1, -1, -1)
start = alignments[0][3] + 4100
end = 8000 - len(alignments[0][0][alignments[0][4]:])
taiForest3 = TaiForest_genome[start:end]
with open("genes/taiForestVP40.txt", "w") as text_file:
    text_file.write(taiForest3)
print(str(alignments[0][3] + 4100) + "  " + str(8000 - len(alignments[0][0][alignments[0][4]:])))

#       same for Bundibugyo genome
alignments = pairwise2.align.localms(Bundibugyo_genome[4100:8000], Marburg_Genes[2], 1, -1, -1, -1)
start = alignments[0][3] + 4100
end = 8000 - len(alignments[0][0][alignments[0][4]:])
bundibugyo3 = Bundibugyo_genome[start:end]
with open("genes/bundibugyoVP40.txt", "w") as text_file:
    text_file.write(bundibugyo3)
print(str(alignments[0][3] + 4100) + "  " + str(8000 - len(alignments[0][0][alignments[0][4]:])))


#           finding the fourth gene in the genomes : GP

#       local alignment to find the gene, interval: 5500-9000 (5500: a little bit less than end
#                                                              points found in previuse gene)

alignments = pairwise2.align.localms(Zaire_genome[5500:9000], Marburg_Genes[3], 1, -1, -1, -1)
start = alignments[0][3] + 5500
end = 9000 - len(alignments[0][0][alignments[0][4]:])
zaire4 = Zaire_genome[start:end]
with open("genes/zairGP.txt", "w") as text_file:
    text_file.write(zaire4)
print(str(alignments[0][3] + 5500) + "  " + str(9000 - len(alignments[0][0][alignments[0][4]:])))

#       same for Sudan genome
alignments = pairwise2.align.localms(Sudan_genome[5500:9000], Marburg_Genes[3], 1, -1, -1, -1)
start = alignments[0][3] + 5500
end = 9000 - len(alignments[0][0][alignments[0][4]:])
sudan4 = Sudan_genome[start:end]
with open("genes/sudanGP.txt", "w") as text_file:
    text_file.write(sudan4)
print(str(alignments[0][3] + 5500) + "  " + str(9000 - len(alignments[0][0][alignments[0][4]:])))

#       same for Reston genome
alignments = pairwise2.align.localms(Reston_genome[5500:9000], Marburg_Genes[3], 1, -1, -1, -1)
start = alignments[0][3] + 5500
end = 9000 - len(alignments[0][0][alignments[0][4]:])
reston4 = Reston_genome[start:end]
with open("genes/restonGP.txt", "w") as text_file:
    text_file.write(reston4)
print(str(alignments[0][3] + 5500) + "  " + str(9000 - len(alignments[0][0][alignments[0][4]:])))

#       same for TaiForest genome                    5300: last end point was different from others
alignments = pairwise2.align.localms(TaiForest_genome[5300:9000], Marburg_Genes[3], 1, -1, -1, -1)
start = alignments[0][3] + 5300
end = 9000 - len(alignments[0][0][alignments[0][4]:])
taiForest4 = TaiForest_genome[start:end]
with open("genes/taiForestGP.txt", "w") as text_file:
    text_file.write(taiForest4)
print(str(alignments[0][3] + 5300) + "  " + str(9000 - len(alignments[0][0][alignments[0][4]:])))

#       same for Bundibugyo genome                    5900: last end point was different from others
alignments = pairwise2.align.localms(Bundibugyo_genome[5900:10000], Marburg_Genes[3], 1, -1, -1, -1)
start = alignments[0][3] + 5900
end = 10000 - len(alignments[0][0][alignments[0][4]:])
bundibugyo4 = Bundibugyo_genome[start:end]
with open("genes/bundibugyoGP.txt", "w") as text_file:
    text_file.write(bundibugyo4)
print(str(alignments[0][3] + 5900) + "  " + str(10000 - len(alignments[0][0][alignments[0][4]:])))


#           finding the fifth gene in the genomes : VP30

#       local alignment to find the gene, interval: 8500-11000 (8500: a little bit less than end
#                                                              points found in previuse gene)

alignments = pairwise2.align.localms(Zaire_genome[8500:11000], Marburg_Genes[4], 1, -1, -1, -1)
start = alignments[0][3] + 8500
end = 11000 - len(alignments[0][0][alignments[0][4]:])
zaire5 = Zaire_genome[start:end]
with open("genes/zairVP30.txt", "w") as text_file:
    text_file.write(zaire5)
print(str(alignments[0][3] + 8500) + "  " + str(11000 - len(alignments[0][0][alignments[0][4]:])))

#       same for Sudan genome
alignments = pairwise2.align.localms(Sudan_genome[8500:11000], Marburg_Genes[4], 1, -1, -1, -1)
start = alignments[0][3] + 8500
end = 11000 - len(alignments[0][0][alignments[0][4]:])
sudan5 = Sudan_genome[start:end]
with open("genes/sudanVP30.txt", "w") as text_file:
    text_file.write(sudan5)
print(str(alignments[0][3] + 8500) + "  " + str(11000 - len(alignments[0][0][alignments[0][4]:])))

#       same for Reston genome
alignments = pairwise2.align.localms(Reston_genome[8500:11000], Marburg_Genes[4], 1, -1, -1, -1)
start = alignments[0][3] + 8500
end = 11000 - len(alignments[0][0][alignments[0][4]:])
reston5 = Reston_genome[start:end]
with open("genes/restonVP30.txt", "w") as text_file:
    text_file.write(reston5)
print(str(alignments[0][3] + 8500) + "  " + str(11000 - len(alignments[0][0][alignments[0][4]:])))

#       same for TaiForest genome
alignments = pairwise2.align.localms(TaiForest_genome[8500:11000], Marburg_Genes[4], 1, -1, -1, -1)
start = alignments[0][3] + 8500
end = 11000 - len(alignments[0][0][alignments[0][4]:])
taiForest5 = TaiForest_genome[start:end]
with open("genes/taiForestVP30.txt", "w") as text_file:
    text_file.write(taiForest5)
print(str(alignments[0][3] + 8500) + "  " + str(11000 - len(alignments[0][0][alignments[0][4]:])))

#       same for Bundibugyo genome
alignments = pairwise2.align.localms(Bundibugyo_genome[8500:11000], Marburg_Genes[4], 1, -1, -1, -1)
start = alignments[0][3] + 8500
end = 11000 - len(alignments[0][0][alignments[0][4]:])
bundibugyo5 = Bundibugyo_genome[start:end]
with open("genes/bundibugyoVP30.txt", "w") as text_file:
    text_file.write(bundibugyo5)
print(str(alignments[0][3] + 8500) + "  " + str(11000 - len(alignments[0][0][alignments[0][4]:])))


#           finding the sixth gene in the genomes : VP24

#       local alignment to find the gene, interval: 9500-12000 (9500: a little bit less than end
#                                                              points found in previuse gene)

alignments = pairwise2.align.localms(Zaire_genome[9500:12000], Marburg_Genes[5], 1, -1, -1, -1)
start = alignments[0][3] + 9500
end = 12000 - len(alignments[0][0][alignments[0][4]:])
zaire6 = Zaire_genome[start:end]
with open("genes/zairVP24.txt", "w") as text_file:
    text_file.write(zaire6)
print(str(alignments[0][3] + 9500) + "  " + str(12000 - len(alignments[0][0][alignments[0][4]:])))

#       same for Sudan genome
alignments = pairwise2.align.localms(Sudan_genome[9500:12000], Marburg_Genes[5], 1, -1, -1, -1)
start = alignments[0][3] + 9500
end = 12000 - len(alignments[0][0][alignments[0][4]:])
sudan6 = Sudan_genome[start:end]
with open("genes/sudanVP24.txt", "w") as text_file:
    text_file.write(sudan6)
print(str(alignments[0][3] + 9500) + "  " + str(12000 - len(alignments[0][0][alignments[0][4]:])))

#       same for Reston genome
alignments = pairwise2.align.localms(Reston_genome[9500:12000], Marburg_Genes[5], 1, -1, -1, -1)
start = alignments[0][3] + 9500
end = 12000 - len(alignments[0][0][alignments[0][4]:])
reston6 = Reston_genome[start:end]
with open("genes/restonVP24.txt", "w") as text_file:
    text_file.write(reston6)
print(str(alignments[0][3] + 9500) + "  " + str(12000 - len(alignments[0][0][alignments[0][4]:])))

#       same for TaiForest genome
alignments = pairwise2.align.localms(TaiForest_genome[9500:12000], Marburg_Genes[5], 1, -1, -1, -1)
start = alignments[0][3] + 9500
end = 12000 - len(alignments[0][0][alignments[0][4]:])
taiForest6 = TaiForest_genome[start:end]
with open("genes/taiForestVP24.txt", "w") as text_file:
    text_file.write(taiForest6)
print(str(alignments[0][3] + 9500) + "  " + str(12000 - len(alignments[0][0][alignments[0][4]:])))

#       same for Bundibugyo genome
alignments = pairwise2.align.localms(Bundibugyo_genome[9500:12000], Marburg_Genes[5], 1, -1, -1, -1)
start = alignments[0][3] + 9500
end = 12000 - len(alignments[0][0][alignments[0][4]:])
bundibugyo6 = Bundibugyo_genome[start:end]
with open("genes/bundibugyoVP24.txt", "w") as text_file:
    text_file.write(bundibugyo6)
print(str(alignments[0][3] + 9500) + "  " + str(12000 - len(alignments[0][0][alignments[0][4]:])))


#           finding the seventh gene in the genomes : L

#       local alignment to find the gene, interval: 11000-15000 (11000: a little bit less than end
#                                                              points found in previuse gene)

alignments = pairwise2.align.localms(Zaire_genome[11000:15000], Marburg_Genes[6][0:2500], 1, -1, -1, -1)
alignments1 = pairwise2.align.localms(Zaire_genome[15000:], Marburg_Genes[6][5000:], 1, -1, -1, -1)
start = alignments[0][3] + 11000
end = len(Zaire_genome) - len(alignments1[0][0][alignments1[0][4]:])
zaire7 = Zaire_genome[start:end]
with open("genes/zairL.txt", "w") as text_file:
    text_file.write(zaire7)
print(str(alignments[0][3] + 11000) + "  " + str(len(Zaire_genome) - len(alignments1[0][0][alignments1[0][4]:])))

#       same for Sudan genome
alignments = pairwise2.align.localms(Sudan_genome[11000:15000], Marburg_Genes[6][0:2500], 1, -1, -1, -1)
alignments1 = pairwise2.align.localms(Sudan_genome[15000:], Marburg_Genes[6][5000:], 1, -1, -1, -1)
start = alignments[0][3] + 11000
end = len(Sudan_genome) - len(alignments1[0][0][alignments1[0][4]:])
sudan7 = Sudan_genome[start:end]
with open("genes/sudanL.txt", "w") as text_file:
    text_file.write(sudan7)
print(str(alignments[0][3] + 11000) + "  " + str(len(Sudan_genome) - len(alignments1[0][0][alignments1[0][4]:])))

#       same for Reston genome
alignments = pairwise2.align.localms(Reston_genome[11000:15000], Marburg_Genes[6][0:2500], 1, -1, -1, -1)
alignments1 = pairwise2.align.localms(Reston_genome[15000:], Marburg_Genes[6][5000:], 1, -1, -1, -1)
start = alignments[0][3] + 11000
end = len(Reston_genome) - len(alignments1[0][0][alignments1[0][4]:])
reston7 = Reston_genome[start:end]
with open("genes/restonL.txt", "w") as text_file:
    text_file.write(reston7)
print(str(alignments[0][3] + 11000) + "  " + str(len(Reston_genome) - len(alignments1[0][0][alignments1[0][4]:])))

#       same for TaiForest genome
alignments = pairwise2.align.localms(TaiForest_genome[11000:15000], Marburg_Genes[6][0:2500], 1, -1, -1, -1)
alignments1 = pairwise2.align.localms(TaiForest_genome[15000:], Marburg_Genes[6][5000:], 1, -1, -1, -1)
start = alignments[0][3] + 11000
end = len(TaiForest_genome) - len(alignments1[0][0][alignments1[0][4]:])
taiForest7 = TaiForest_genome[start:end]
with open("genes/taiForestL.txt", "w") as text_file:
    text_file.write(taiForest7)
print(str(alignments[0][3] + 11000) + "  " + str(len(TaiForest_genome) - len(alignments1[0][0][alignments1[0][4]:])))

#       same for Bundibugyo genome
alignments = pairwise2.align.localms(Bundibugyo_genome[11000:15000], Marburg_Genes[6][0:2500], 1, -1, -1, -1)
alignments1 = pairwise2.align.localms(Bundibugyo_genome[15000:], Marburg_Genes[6][5000:], 1, -1, -1, -1)
start = alignments[0][3] + 11000
end = len(Bundibugyo_genome) - len(alignments1[0][0][alignments1[0][4]:])
bundibugyo7 = Bundibugyo_genome[start:end]
with open("genes/bundibugyoL.txt", "w") as text_file:
    text_file.write(bundibugyo7)
print(str(alignments[0][3] + 11000) + "  " + str(len(Bundibugyo_genome) - len(alignments1[0][0][alignments1[0][4]:])))



