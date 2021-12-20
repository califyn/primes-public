import os
import matplotlib as mpl
mpl.use('tkagg')
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
#from sklearn.manifold import TSNE
from openTSNE import TSNE



DATA_DIR = '../../tcga/data'
FNAME = 'sample_matrix.txt'
CANCERS = ['BLCA', 'HNSC', 'LIHC']

def groupsep(i):
    """if i in ["LGG", "GBM"]:
        return 0
    elif i in ["OV", "UCEC", "CESC"]:
        return 1
    elif i in ["HNSC", "LUSC", "ESCA", "STAD"]:
        return 2
    elif i in ["BLCA", "BRCA", "LUAD", "LIHC", "SKCM", "SARC", "ACC", "COAD", "PAAD"]:
        return 3
    elif i in ["AML", "THCA", "KIRP", "KIRC", "PCPG", "THYM", "TEST"]:
        return 4"""

    if i in ["AML", "THCA", "KIRP", "KIRC", "PCPG", "THYM", "TEST"]:
        return 0
    else:
        return 1

mutations = []
names = []
genes = []
for subdir, dirs, files in os.walk(DATA_DIR):
    if dirs != []:
        continue

    # Skip first two columns of data frame, convert to
    # binary numpy matrix, and transpose.
    fname = os.path.join(subdir, 'mutations_transposed.txt')
    dataframe = pd.read_csv(fname, sep='\t')
    if len(genes) == 0:
        genes = list(dataframe.iloc[1:, 0])
    mut = dataframe.iloc[1:, 1:].to_numpy(dtype='str')
    mut = np.where(mut == 'WT', 0, 1)
    print(mut)

    #if len(genes) == 0:
    #    genes = list(dataframe.columns[2:])
    #else:
    #    assert(genes == list(dataframe.columns[2:])) # Make sure we have the same genes across samples.
    #print(genes)

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
patients = np.concatenate(tuple(mutations), axis=1).T

shrunk_patients = np.load("patients_tsne-interim-mut.npy")

num_mutated = np.sum(np.abs(patients[:, 763:]), axis=1)
dist = np.histogram(num_mutated, bins=np.max(num_mutated))
print(dist)
print(np.sum(dist[0]))

markers = [0]
#num_markers = 5
#markers.append(np.quantile(num_mutated, 0.9).item()**(1.0 / (num_markers - 1))) # slightly incorrect b/c of the +1 thing
#for i in range(2, num_markers):
#    markers.append(markers[i - 1] * markers[1])
#input(np.quantile(num_mutated, 0.9))
max = np.quantile(num_mutated, 0.9)
num_mutated = np.log(num_mutated + 1)
#input(np.quantile(num_mutated, 0.9))
num_mutated = num_mutated / np.quantile(num_mutated, 0.9)

#input(np.sum(patients[:, 92]))

#num_mutated = num_mutated / 10

fig, axs = plt.subplots(figsize=(8,5))
plt.subplots_adjust(right=0.75)

color = []
label = []
#cancer_list = ["LGG", "SKCM", "COAD", "LUSC", "LUAD", "BLCA", "BRCA", "LIHC", "SARC", "PAAD"]
gene_num = 71
#for i in range(0, np.size(num_mutated, axis=0)):
#    color.append(plt.get_cmap('viridis')(int(min([254, (patients[i, gene_num] + patients[i, gene_num + 763*1] + patients[i, gene_num + 763*2]) * 255]))))
for i in range(0, np.size(num_mutated, axis=0)):
    color.append(plt.get_cmap('viridis')(int(min([254, num_mutated[i] * 255]))))
"""for i in range(0, len(mutations)):
    print(np.size(mutations[i], axis=1))
    for j in range(0, np.size(mutations[i], axis=1)):
        if names[i] in cancer_list:
            color.append(plt.get_cmap('tab10')(cancer_list.index(names[i])))
            label.append(names[i])
        else:
            color.append([0.9, 0.9, 0.9, 0])
            label.append("Other")"""
    #axs.scatter(0, i, c=color[len(color) - 2:len(color) - 1], s=1, label=names[i])
        #color.append(plt.get_cmap('tab10'))

#for i in range(0, len(mutations)):
#    for j in range(0, np.size(mutations[i], axis=1)):
#        color.append(plt.get_cmap('tab20')(i % 20))
#    axs.scatter(0, 0, c=color[len(color) - 2:len(color) - 1], s=2, label=names[i])

#axs.scatter([0, 0], [0, 1])
axs.scatter(shrunk_patients[:, 0], shrunk_patients[:, 1], s=1, c=color)
#for i, txt in enumerate(names):
#    ax.annotate(txt, (z[i], y[i]))
#axs.legend(bbox_to_anchor=(1.05, 1),ncol=2)


cmap = mpl.cm.viridis
norm = mpl.colors.LogNorm(vmin=1, vmax=max)
#norm = mpl.colors.BoundaryNorm([0, 0.5, 1], cmap.N)

fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap),
             ax=axs, label='CNA')

plt.show()
