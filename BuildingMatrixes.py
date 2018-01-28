import numpy as np
from Bio import pairwise2

#       reading NP genes of viruses from the file written in GenesFinding

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

#           global alignment between NP genes of all 5 viruses, match scre = zero to
#           find the edit distance (the score is multiplied by -1 at the end)

a = pairwise2.align.globalms(zaire1, sudan1, 0, -1, -1, -1)
zaire_sudan1 = a[0][2] * -1
a = pairwise2.align.globalms(zaire1, reston1, 0, -1, -1, -1)
zaire_restone1 = a[0][2] * -1
a = pairwise2.align.globalms(zaire1, taiForest1, 0, -1, -1, -1)
zaire_taiforest1 = a[0][2] * -1
a = pairwise2.align.globalms(zaire1, bundibugyo1, 0, -1, -1, -1)
zaire_bundibugyo1 = a[0][2] * -1

a = pairwise2.align.globalms(sudan1, reston1, 0, -1, -1, -1)
sudan_restone1 = a[0][2] * -1
a = pairwise2.align.globalms(sudan1, taiForest1, 0, -1, -1, -1)
sudan_taiforest1 = a[0][2] * -1
a = pairwise2.align.globalms(sudan1, bundibugyo1, 0, -1, -1, -1)
sudan_bundibugyo1 = a[0][2] * -1

a = pairwise2.align.globalms(reston1, taiForest1, 0, -1, -1, -1)
reston_taiforest1 = a[0][2] * -1
a = pairwise2.align.globalms(reston1, bundibugyo1, 0, -1, -1, -1)
reston_bundibugyo1 = a[0][2] * -1

a = pairwise2.align.globalms(taiForest1, bundibugyo1, 0, -1, -1, -1)
taiForest_bundibugyo1 = a[0][2] * -1


#       reading VP35 genes of viruses from the file written in GenesFinding

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

#           global alignment between VP35 genes of all 5 viruses, match scre = zero to
#           find the edit distance (the score is multiplied by -1 at the end)

a = pairwise2.align.globalms(zaire2, sudan2, 0, -1, -1, -1)
zaire_sudan2 = a[0][2] * -1
a = pairwise2.align.globalms(zaire2, reston2, 0, -1, -1, -1)
zaire_restone2 = a[0][2] * -1
a = pairwise2.align.globalms(zaire2, taiForest2, 0, -1, -1, -1)
zaire_taiforest2 = a[0][2] * -1
a = pairwise2.align.globalms(zaire2, bundibugyo2, 0, -1, -1, -1)
zaire_bundibugyo2 = a[0][2] * -1

a = pairwise2.align.globalms(sudan2, reston2, 0, -1, -1, -1)
sudan_restone2 = a[0][2] * -1
a = pairwise2.align.globalms(sudan2, taiForest2, 0, -1, -1, -1)
sudan_taiforest2 = a[0][2] * -1
a = pairwise2.align.globalms(sudan2, bundibugyo2, 0, -1, -1, -1)
sudan_bundibugyo2 = a[0][2] * -1

a = pairwise2.align.globalms(reston2, taiForest2, 0, -1, -1, -1)
reston_taiforest2 = a[0][2] * -1
a = pairwise2.align.globalms(reston2, bundibugyo2, 0, -1, -1, -1)
reston_bundibugyo2 = a[0][2] * -1

a = pairwise2.align.globalms(taiForest2, bundibugyo2, 0, -1, -1, -1)
taiForest_bundibugyo2 = a[0][2] * -1


#       reading VP40 genes of viruses from the file written in GenesFinding

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

#           global alignment between VP40 genes of all 5 viruses, match scre = zero to
#           find the edit distance (the score is multiplied by -1 at the end)

a = pairwise2.align.globalms(zaire3, sudan3, 0, -1, -1, -1)
zaire_sudan3 = a[0][2] * -1
a = pairwise2.align.globalms(zaire3, reston3, 0, -1, -1, -1)
zaire_restone3 = a[0][2] * -1
a = pairwise2.align.globalms(zaire3, taiForest3, 0, -1, -1, -1)
zaire_taiforest3 = a[0][2] * -1
a = pairwise2.align.globalms(zaire3, bundibugyo3, 0, -1, -1, -1)
zaire_bundibugyo3 = a[0][2] * -1

a = pairwise2.align.globalms(sudan3, reston3, 0, -1, -1, -1)
sudan_restone3 = a[0][2] * -1
a = pairwise2.align.globalms(sudan3, taiForest3, 0, -1, -1, -1)
sudan_taiforest3 = a[0][2] * -1
a = pairwise2.align.globalms(sudan3, bundibugyo3, 0, -1, -1, -1)
sudan_bundibugyo3 = a[0][2] * -1

a = pairwise2.align.globalms(reston3, taiForest3, 0, -1, -1, -1)
reston_taiforest3 = a[0][2] * -1
a = pairwise2.align.globalms(reston3, bundibugyo3, 0, -1, -1, -1)
reston_bundibugyo3 = a[0][2] * -1

a = pairwise2.align.globalms(taiForest3, bundibugyo3, 0, -1, -1, -1)
taiForest_bundibugyo3 = a[0][2] * -1


#       reading GP genes of viruses from the file written in GenesFinding

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

#           global alignment between GP genes of all 5 viruses, match scre = zero to
#           find the edit distance (the score is multiplied by -1 at the end)

a = pairwise2.align.globalms(zaire4, sudan4, 0, -1, -1, -1)
zaire_sudan4 = a[0][2] * -1
a = pairwise2.align.globalms(zaire4, reston4, 0, -1, -1, -1)
zaire_restone4 = a[0][2] * -1
a = pairwise2.align.globalms(zaire4, taiForest4, 0, -1, -1, -1)
zaire_taiforest4 = a[0][2] * -1
a = pairwise2.align.globalms(zaire4, bundibugyo4, 0, -1, -1, -1)
zaire_bundibugyo4 = a[0][2] * -1

a = pairwise2.align.globalms(sudan4, reston4, 0, -1, -1, -1)
sudan_restone4 = a[0][2] * -1
a = pairwise2.align.globalms(sudan4, taiForest4, 0, -1, -1, -1)
sudan_taiforest4 = a[0][2] * -1
a = pairwise2.align.globalms(sudan4, bundibugyo4, 0, -1, -1, -1)
sudan_bundibugyo4 = a[0][2] * -1

a = pairwise2.align.globalms(reston4, taiForest4, 0, -1, -1, -1)
reston_taiforest4 = a[0][2] * -1
a = pairwise2.align.globalms(reston4, bundibugyo4, 0, -1, -1, -1)
reston_bundibugyo4 = a[0][2] * -1

a = pairwise2.align.globalms(taiForest4, bundibugyo4, 0, -1, -1, -1)
taiForest_bundibugyo4 = a[0][2] * -1


#       reading VP30 genes of viruses from the file written in GenesFinding

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

#           global alignment between VP30 genes of all 5 viruses, match scre = zero to
#           find the edit distance (the score is multiplied by -1 at the end)

a = pairwise2.align.globalms(zaire5, sudan5, 0, -1, -1, -1)
zaire_sudan5 = a[0][2] * -1
a = pairwise2.align.globalms(zaire5, reston5, 0, -1, -1, -1)
zaire_restone5 = a[0][2] * -1
a = pairwise2.align.globalms(zaire5, taiForest5, 0, -1, -1, -1)
zaire_taiforest5 = a[0][2] * -1
a = pairwise2.align.globalms(zaire5, bundibugyo5, 0, -1, -1, -1)
zaire_bundibugyo5 = a[0][2] * -1

a = pairwise2.align.globalms(sudan5, reston5, 0, -1, -1, -1)
sudan_restone5 = a[0][2] * -1
a = pairwise2.align.globalms(sudan5, taiForest5, 0, -1, -1, -1)
sudan_taiforest5 = a[0][2] * -1
a = pairwise2.align.globalms(sudan5, bundibugyo5, 0, -1, -1, -1)
sudan_bundibugyo5 = a[0][2] * -1

a = pairwise2.align.globalms(reston5, taiForest5, 0, -1, -1, -1)
reston_taiforest5 = a[0][2] * -1
a = pairwise2.align.globalms(reston5, bundibugyo5, 0, -1, -1, -1)
reston_bundibugyo5 = a[0][2] * -1

a = pairwise2.align.globalms(taiForest5, bundibugyo5, 0, -1, -1, -1)
taiForest_bundibugyo5 = a[0][2] * -1


#       reading VP24 genes of viruses from the file written in GenesFinding

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

#           global alignment between VP24 genes of all 5 viruses, match scre = zero to
#           find the edit distance (the score is multiplied by -1 at the end)

a = pairwise2.align.globalms(zaire6, sudan6, 0, -1, -1, -1)
zaire_sudan6 = a[0][2] * -1
a = pairwise2.align.globalms(zaire6, reston6, 0, -1, -1, -1)
zaire_restone6 = a[0][2] * -1
a = pairwise2.align.globalms(zaire6, taiForest6, 0, -1, -1, -1)
zaire_taiforest6 = a[0][2] * -1
a = pairwise2.align.globalms(zaire6, bundibugyo6, 0, -1, -1, -1)
zaire_bundibugyo6 = a[0][2] * -1

a = pairwise2.align.globalms(sudan6, reston6, 0, -1, -1, -1)
sudan_restone6 = a[0][2] * -1
a = pairwise2.align.globalms(sudan6, taiForest6, 0, -1, -1, -1)
sudan_taiforest6 = a[0][2] * -1
a = pairwise2.align.globalms(sudan6, bundibugyo6, 0, -1, -1, -1)
sudan_bundibugyo6 = a[0][2] * -1

a = pairwise2.align.globalms(reston6, taiForest6, 0, -1, -1, -1)
reston_taiforest6 = a[0][2] * -1
a = pairwise2.align.globalms(reston6, bundibugyo6, 0, -1, -1, -1)
reston_bundibugyo6 = a[0][2] * -1

a = pairwise2.align.globalms(taiForest6, bundibugyo6, 0, -1, -1, -1)
taiForest_bundibugyo6 = a[0][2] * -1


#       reading L genes of viruses from the file written in GenesFinding

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

#           global alignment between L genes of all 5 viruses, match scre = zero to
#           find the edit distance (the score is multiplied by -1 at the end)

a = pairwise2.align.globalms(zaire7, sudan7, 0, -1, -1, -1)
zaire_sudan7 = a[0][2] * -1
a = pairwise2.align.globalms(zaire7, reston7, 0, -1, -1, -1)
zaire_restone7 = a[0][2] * -1
a = pairwise2.align.globalms(zaire7, taiForest7, 0, -1, -1, -1)
zaire_taiforest7 = a[0][2] * -1
a = pairwise2.align.globalms(zaire7, bundibugyo7, 0, -1, -1, -1)
zaire_bundibugyo7 = a[0][2] * -1

a = pairwise2.align.globalms(sudan7, reston7, 0, -1, -1, -1)
sudan_restone7 = a[0][2] * -1
a = pairwise2.align.globalms(sudan7, taiForest7, 0, -1, -1, -1)
sudan_taiforest7 = a[0][2] * -1
a = pairwise2.align.globalms(sudan7, bundibugyo7, 0, -1, -1, -1)
sudan_bundibugyo7 = a[0][2] * -1

a = pairwise2.align.globalms(reston7, taiForest7, 0, -1, -1, -1)
reston_taiforest7 = a[0][2] * -1
a = pairwise2.align.globalms(reston7, bundibugyo7, 0, -1, -1, -1)
reston_bundibugyo7 = a[0][2] * -1

a = pairwise2.align.globalms(taiForest7, bundibugyo7, 0, -1, -1, -1)
taiForest_bundibugyo7 = a[0][2] * -1

#           building the NP distance matrix from the calculated distances

matrix1 = np.asarray([[0, zaire_sudan1, zaire_restone1, zaire_taiforest1, zaire_bundibugyo1], # zaire
                      [zaire_sudan1, 0, sudan_restone1, sudan_taiforest1, sudan_bundibugyo1], # sudan
                      [zaire_restone1, sudan_restone1, 0, reston_taiforest1, reston_bundibugyo1], # reston
                      [zaire_taiforest1, sudan_taiforest1, reston_taiforest1, 0, taiForest_bundibugyo1], # tiforest
                      [zaire_bundibugyo1, sudan_bundibugyo1, reston_bundibugyo1, taiForest_bundibugyo1, 0]  # bundibugyo
                     ])
#           writing the matrix to a csv file
np.savetxt("csv/NP.csv", matrix1, delimiter=",")


#           building the VP35 distance matrix from the calculated distances

matrix2 = np.asarray([[0, zaire_sudan2, zaire_restone2, zaire_taiforest2, zaire_bundibugyo2], # zaire
                      [zaire_sudan2, 0, sudan_restone2, sudan_taiforest2, sudan_bundibugyo2], # sudan
                      [zaire_restone2, sudan_restone2, 0, reston_taiforest2, reston_bundibugyo2], # reston
                      [zaire_taiforest2, sudan_taiforest2, reston_taiforest2, 0, taiForest_bundibugyo2], # tiforest
                      [zaire_bundibugyo2, sudan_bundibugyo2, reston_bundibugyo2, taiForest_bundibugyo2, 0]  # bundibugyo
                     ])
#           writing the matrix to a csv file
np.savetxt("csv/VP35.csv", matrix2, delimiter=",")


#           building the VP40 distance matrix from the calculated distances

matrix3 = np.asarray([[0, zaire_sudan3, zaire_restone3, zaire_taiforest3, zaire_bundibugyo3], # zaire
                      [zaire_sudan3, 0, sudan_restone3, sudan_taiforest3, sudan_bundibugyo3], # sudan
                      [zaire_restone3, sudan_restone3, 0, reston_taiforest3, reston_bundibugyo3], # reston
                      [zaire_taiforest3, sudan_taiforest3, reston_taiforest3, 0, taiForest_bundibugyo3], # tiforest
                      [zaire_bundibugyo3, sudan_bundibugyo3, reston_bundibugyo3, taiForest_bundibugyo3, 0]  # bundibugyo
                     ])
#           writing the matrix to a csv file
np.savetxt("csv/VP40.csv", matrix3, delimiter=",")


#           building the GP distance matrix from the calculated distances

matrix4 = np.asarray([[0, zaire_sudan4, zaire_restone4, zaire_taiforest4, zaire_bundibugyo4], # zaire
                      [zaire_sudan4, 0, sudan_restone4, sudan_taiforest4, sudan_bundibugyo4], # sudan
                      [zaire_restone4, sudan_restone4, 0, reston_taiforest4, reston_bundibugyo4], # reston
                      [zaire_taiforest4, sudan_taiforest4, reston_taiforest4, 0, taiForest_bundibugyo4], # tiforest
                      [zaire_bundibugyo4, sudan_bundibugyo4, reston_bundibugyo4, taiForest_bundibugyo4, 0]  # bundibugyo
                     ])
#           writing the matrix to a csv file
np.savetxt("csv/GP.csv", matrix4, delimiter=",")


#           building the VP30 distance matrix from the calculated distances

matrix5 = np.asarray([[0, zaire_sudan5, zaire_restone5, zaire_taiforest5, zaire_bundibugyo5], # zaire
                      [zaire_sudan5, 0, sudan_restone5, sudan_taiforest5, sudan_bundibugyo5], # sudan
                      [zaire_restone5, sudan_restone5, 0, reston_taiforest5, reston_bundibugyo5], # reston
                      [zaire_taiforest5, sudan_taiforest5, reston_taiforest5, 0, taiForest_bundibugyo5], # tiforest
                      [zaire_bundibugyo5, sudan_bundibugyo5, reston_bundibugyo5, taiForest_bundibugyo5, 0]  # bundibugyo
                     ])
#           writing the matrix to a csv file
np.savetxt("csv/VP30.csv", matrix5, delimiter=",")


#           building the VP24 distance matrix from the calculated distances

matrix6 = np.asarray([[0, zaire_sudan6, zaire_restone6, zaire_taiforest6, zaire_bundibugyo6], # zaire
                      [zaire_sudan6, 0, sudan_restone6, sudan_taiforest6, sudan_bundibugyo6], # sudan
                      [zaire_restone6, sudan_restone6, 0, reston_taiforest6, reston_bundibugyo6], # reston
                      [zaire_taiforest6, sudan_taiforest6, reston_taiforest6, 0, taiForest_bundibugyo6], # tiforest
                      [zaire_bundibugyo6, sudan_bundibugyo6, reston_bundibugyo6, taiForest_bundibugyo6, 0]  # bundibugyo
                     ])
#           writing the matrix to a csv file
np.savetxt("csv/VP24.csv", matrix6, delimiter=",")


#           building the L distance matrix from the calculated distances

matrix7 = np.asarray([[0, zaire_sudan7, zaire_restone7, zaire_taiforest7, zaire_bundibugyo7], # zaire
                      [zaire_sudan7, 0, sudan_restone7, sudan_taiforest7, sudan_bundibugyo7], # sudan
                      [zaire_restone7, sudan_restone7, 0, reston_taiforest7, reston_bundibugyo7], # reston
                      [zaire_taiforest7, sudan_taiforest7, reston_taiforest7, 0, taiForest_bundibugyo7], # tiforest
                      [zaire_bundibugyo7, sudan_bundibugyo7, reston_bundibugyo7, taiForest_bundibugyo7, 0]  # bundibugyo
                     ])
#           writing the matrix to a csv file
np.savetxt("csv/L.csv", matrix7, delimiter=",")