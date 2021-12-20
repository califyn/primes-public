import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.linear_model import LinearRegression
#from sklearn.manifold import TSNE
from openTSNE import TSNE

DATA_DIR = '../../tcga/data'
FNAME = 'sample_matrix.txt'

mutations = []
names = []
genes = []
for subdir, dirs, files in os.walk(DATA_DIR):
    if dirs != []:
        continue

    fname = os.path.join(subdir, 'mutations_transposed.txt')
    dataframe = pd.read_csv(fname, sep='\t')
    mut = dataframe.iloc[1:, 1:].to_numpy(dtype='str')
    mut = np.where(mut == 'WT', 0, 1)
    print(mut)

    fname = os.path.join(subdir, 'cna_transposed.txt')
    dataframe = pd.read_csv(fname, sep='\t')
    cna = dataframe.iloc[1:, 1:].to_numpy(dtype='int64')
    cna_pos = np.where(cna > 1, 1, 0)
    cna_neg = np.where(cna < -1, 1, 0)

    full = np.concatenate((mut, cna_pos, cna_neg))

    ngenes, nsamples = full.shape
    print(f'{subdir}: {ngenes} genes x {nsamples} samples')

    names.append(subdir)
    mutations.append(full)


names = [x[len(FNAME) - 1:] for x in names]
names = ["COAD" if x == "COADREAD" else x for x in names]

patients = np.concatenate(tuple(mutations), axis=1)

print(np.shape(patients))
muts = np.sum(patients[:763, :], axis=0)
amps = np.sum(patients[763:1526, :], axis=0)

print(np.quantile(muts, 0.9))
muts = np.minimum(muts, np.quantile(muts, 0.98))
amps = np.minimum(amps, np.quantile(amps, 0.9))

print(np.corrcoef(muts, amps))

#model = LinearRegression(muts, amps)
#model.fit(muts, amps)
#x_new = np.linspace(0, 100, 100)
#y_new = model.predict(x_new[:, np.newaxis])

fig, axs = plt.subplots(figsize=(5,5))
axs.scatter(muts, amps, s=1, alpha=0.1)
#axs.plot(x_new, y_new)
plt.show()
