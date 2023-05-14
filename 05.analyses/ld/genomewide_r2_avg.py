# python script to average across all chromosomes worth of average r2 values per distance between base pairs
# only meant for this one instance of use in ld calculations

import pandas as pd

chr_list = ["Chr01", "Chr02", "Chr03", "Chr04", "Chr05", "Chr06", "Chr07", "Chr08", "Chr09", "Chr10", "Chr11"]

chr01 = open("Chr01_r2.csv", "r")
chr02 = open("Chr02_r2.csv", "r")
chr03 = open("Chr03_r2.csv", "r")
chr04 = open("Chr04_r2.csv", "r")
chr05 = open("Chr05_r2.csv", "r")
chr06 = open("Chr06_r2.csv", "r")
chr07 = open("Chr07_r2.csv", "r")
chr08 = open("Chr08_r2.csv", "r")
chr09 = open("Chr09_r2.csv", "r")
chr10 = open("Chr10_r2.csv", "r")
chr11 = open("Chr11_r2.csv", "r")

# skip header
chr01.readline()
chr02.readline()
chr03.readline()
chr04.readline()
chr05.readline()
chr06.readline()
chr07.readline()
chr08.readline()
chr09.readline()
chr10.readline()
chr11.readline()

genomewide_r2 = pd.DataFrame(data = None, columns = ["dist", "r2"])

for chr01_line in chr01:
    chr01_r2 = float(chr01_line.strip().split(",")[1])
    chr02_r2 = float(chr02.readline().strip().split(",")[1])
    chr03_r2 = float(chr03.readline().strip().split(",")[1])
    chr04_r2 = float(chr04.readline().strip().split(",")[1])
    chr05_r2 = float(chr05.readline().strip().split(",")[1])
    chr06_r2 = float(chr06.readline().strip().split(",")[1])
    chr07_r2 = float(chr07.readline().strip().split(",")[1])
    chr08_r2 = float(chr08.readline().strip().split(",")[1])
    chr09_r2 = float(chr09.readline().strip().split(",")[1])
    chr10_r2 = float(chr10.readline().strip().split(",")[1])
    chr11_r2 = float(chr11.readline().strip().split(",")[1])

    dist = chr01_line.strip().split(",")[0]
    avg_r2 = sum([chr01_r2, chr02_r2, chr03_r2, chr04_r2, chr05_r2, chr06_r2, chr07_r2, chr08_r2, chr09_r2, chr10_r2, chr11_r2])/11

    genomewide_r2.loc[len(genomewide_r2.index)] = [dist, avg_r2]

chr01.close()
chr02.close()
chr03.close()
chr04.close()
chr05.close()
chr06.close()
chr07.close()
chr08.close()
chr09.close()
chr10.close()
chr11.close()

genomewide_r2.to_csv("genomewide_r2.csv", index = False)