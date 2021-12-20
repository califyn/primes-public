import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

DATA_DIR = '../../tcga/data'
FNAME = 'sample_matrix.txt'

mutations = []
names = []
genes = []
for subdir, dirs, files in os.walk(DATA_DIR):
    if dirs != []:
        continue
    fname = os.path.join(subdir, FNAME)
    dataframe = pd.read_csv(fname, sep='\t')
    if len(genes) == 0:
        genes = list(dataframe.columns[2:])
    else:
        assert(genes == list(dataframe.columns[2:])) # Make sure we have the same genes across samples.

    # Skip first two columns of data frame, convert to
    # binary numpy matrix, and transpose.
    mat = dataframe.iloc[:, 2:].to_numpy(dtype='bool').T
    ngenes, nsamples = mat.shape
    # print(f'{subdir}: {ngenes} genes x {nsamples} samples')

    names.append(subdir)
    mutations.append(mat)

names = [x[len(FNAME) - 1:] for x in names]
names = ["COAD" if x == "COADREAD" else x for x in names]





# Basic statistics
input("Press a key to continue.")

names.insert(0, "ALL")
mutations.insert(0, np.concatenate(tuple(mutations), axis=1))

input(np.shape(mutations[0]))

for i in range(0, len(mutations)):
    print("-------------------------------------------------------------------")
    print("Summary of cancer " + str(names[i]))

    print("Starting 20 x 20 subgrid:")
    print(mutations[i][:20, :20])

    print("Type:")
    print(mutations[i].dtype)

    print("Mean:")
    print(np.mean(mutations[i]))

    print("Std:")
    print(np.std(mutations[i]))

    print("10 most common mutated genes:")
    top_10 = np.arange(0, np.size(mutations[i], axis=0))[np.argsort(np.sum(mutations[i], axis=1))[::-1]][:10].astype('int64')
    print(top_10)
    gene_names = []
    for t in top_10:
        gene_names.append(genes[t])
    print(gene_names)

    print("Distribution of how often genes are mutated:")
    hist, bin_edges = np.histogram(np.sum(mutations[i], axis=1), 10)
    for j in range(0, len(hist)):
        print("[" + str(int(bin_edges[j])) + ", " + str(int(bin_edges[j + 1])) + "]: " + str(hist[j]))

    if i != len(mutations) - 1:
        input("Continue...?")
