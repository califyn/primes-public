import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
#from sklearn.manifold import TSNE
from openTSNE import TSNE

DATA_DIR = '../../tcga/data'
FNAME = 'sample_matrix.txt'

mutations = []
names = []
genes = []
"""for subdir, dirs, files in os.walk(DATA_DIR):
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
    # dataframe = dataframe[['ZNRF3', 'PIK3CA', 'PTEN', 'CTNNB1', 'PRKAR1A', 'NPM1', 'FLT3', 'RUNX1', 'DNMT3A', 'CDKN1A', 'KDM6A', 'SDHC', 'PPARG', 'RB1', 'GATA3', 'CDH1', 'KRAS', 'TP53', 'BRAF', 'YAP1', 'CDKN2B', 'CDKN2A', 'APC', 'SOX9', 'SETD1B', 'FGF3', 'ATRX', 'GATA6', 'FHIT', 'EGFR', 'CDK4', 'CHIC2', 'CCND1', 'CASP8', 'NSD1', 'NOTCH1', 'VHL', 'PBRM1', 'FANCD2', 'MET', 'CUL3', 'IDH1', 'IDH2', 'CIC', 'ALB', 'AXIN1', 'KEAP1', 'STK11', 'NKX2-1', 'DCUN1D1', 'LRP1B', 'NFE2L2', 'TERC', 'NDRG1', 'SMAD4', 'HRAS', 'CSDE1', 'NRAS', 'ACVR2A', 'ARID1A', 'KIT', 'GTF2I', 'PPP2R1A', 'PIK3R1', 'SPOP']]
    mat = dataframe.to_numpy(dtype='bool').T
    ngenes, nsamples = mat.shape

    print(f'{subdir}: {ngenes} genes x {nsamples} samples')

    names.append(subdir)
    mutations.append(mat.T)"""
for subdir, dirs, files in os.walk(DATA_DIR):
    if dirs != []:
        continue

    #if len(genes) == 0:
    #    genes = list(dataframe.columns[2:])
    #else:
    #    assert(genes == list(dataframe.columns[2:])) # Make sure we have the same genes across samples.

    # Skip first two columns of data frame, convert to
    # binary numpy matrix, and transpose.
    fname = os.path.join(subdir, 'mutations_transposed.txt')
    dataframe = pd.read_csv(fname, sep='\t')
    mut = dataframe.iloc[1:, 1:].to_numpy(dtype='str')
    mut = np.where(mut == 'WT', 0, 1)
    print(mut)

    fname = os.path.join(subdir, 'cna_transposed.txt')
    dataframe = pd.read_csv(fname, sep='\t')
    cna = dataframe.iloc[1:, 1:].to_numpy(dtype='int64')
    cna_pos = np.where(cna > 0, 1, 0)
    cna_neg = np.where(cna < 0, -1, 0)

    full = np.concatenate((mut, cna_pos, cna_neg))

    ngenes, nsamples = full.shape
    print(f'{subdir}: {ngenes} genes x {nsamples} samples')

    names.append(subdir)
    mutations.append(full)

names = [x[len(FNAME) - 1:] for x in names]
names = ["COAD" if x == "COADREAD" else x for x in names]

#for i in range(10, len(mutations)):
#    del names[10]
#    del mutations[10]

# Signature creation
low = []
high = []
signatures = []
input(np.shape(mutations[0]))
for i in range(0, len(mutations)):
    #low.append(np.quantile(mutations[i], 0.1, axis=1))
    #high.append(np.quantile(mutations[i], 0.9, axis=1))
    signatures.append(np.mean(mutations[i], axis=1))
signatures = np.array(signatures)
input(np.shape(signatures))
#signatures = np.concatenate((low, mutations, high), axis=1)

# PCA
# pca = PCA(n_components=2)
# shrunk_signatures = pca.fit(signatures).transform(signatures).T
# print(np.shape(shrunk_signatures))

# t-SNE
pca = PCA(n_components=24)

#transpose_mut = []
#for mutation in mutations:
#    transpose_mut.append(mutation.T)
#v = np.concatenate(tuple(transpose_mut))


shrunk_signatures = pca.fit_transform(signatures).T

#shrunk_muts = []
#for mutation in mutations:
#    shrunk_muts.append(pca.transform(mutation))
#manifold = TSNE(n_components=2, init='random', perplexity=4)
tsne = TSNE(
    perplexity=3,
    metric="euclidean",
    n_jobs=8,
    random_state=42,
)

shrunk_signatures = tsne.fit(shrunk_signatures.T).T
print(np.shape(shrunk_signatures))
"""print(np.shape(shrunk_signatures))
input(type(shrunk_signatures))
for i in range(0, len(mutations)):
    shrunk_muts[i] = shrunk_signatures.transform(shrunk_muts[i])
shrunk_signatures = shrunk_signatures.T"""

# Display using matplotlib
fig, axs = plt.subplots(figsize=(5,5))
# axs.scatter(shrunk_signatures[0], shrunk_signatures[1])
for i in range(0, len(names)):
    #axs.scatter(np.concatenate((np.array([shrunk_signatures[0, i]]), shrunk_muts[i][:, 0])), np.concatenate((np.array([shrunk_signatures[1, i]]), shrunk_muts[i][:, 1])), c=[plt.get_cmap('tab10')(i % 20)] * (len(shrunk_muts[i]) + 1), s=[1] + [0.05] * len(shrunk_muts[i]))
    axs.scatter(shrunk_signatures[0, i], shrunk_signatures[1, i])
    axs.text(shrunk_signatures[0, i] * (1 + 0.001), shrunk_signatures[1, i] * (1 + 0.001), names[i])
plt.show()
