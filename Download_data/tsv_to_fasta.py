## download dataset from iReceptor
## extract only the sequence ID and the sequence column of interest
## for the input file
import argparse
import pandas as pd

parser = argparse.ArgumentParser()
# Add an argument
parser.add_argument('-i', '--input_tsv', type=str, required=True)
parser.add_argument('-o', '--output', type=str, required=True)
# Parse the argument
args = parser.parse_args()

x = open(args.input_tsv, "r")
z = x.readlines()
f = open(args.output, "w")
columns = z[0].split("\t")
id = columns.index("sequence_id")
sequence = columns.index("sequence")
junction = columns.index("junction")

l = len(z)
for i in range(1, l, 1):
    line = z[i].split("\t")
    id_i = line[id]
    seq = line[sequence]
    if seq == "" and line[junction] != "":
        seq = line[junction]
    f.write(">" + id_i + "\n" + seq + "\n")
f.close()
