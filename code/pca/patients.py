import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
#from sklearn.manifold import TSNE
from openTSNE import TSNE

DATA_DIR = '../../tcga/data'
FNAME = 'sample_matrix.txt'
CANCERS = ['ESCA', 'STAD', 'HNSC']

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
    cna_pos = np.where(cna > 1, 1, 0)
    cna_neg = np.where(cna < -1, 1, 0)

    full = np.concatenate((mut, cna_pos, cna_neg))
    #full = mut

    ngenes, nsamples = full.shape
    print(f'{subdir}: {ngenes} genes x {nsamples} samples')

    names.append(subdir)
    mutations.append(full)
"""for subdir, dirs, files in os.walk(DATA_DIR):
    if dirs != []:
        continue
    fname = os.path.join(subdir, FNAME)
    dataframe = pd.read_csv(fname, sep='\t')
    # Skip first two columns of data frame, convert to
    # binary numpy matrix, and transpose.
    mat = dataframe.iloc[:, 2:].to_numpy(dtype='bool').T
    ngenes, nsamples = mat.shape
    print(f'{subdir}: {ngenes} genes x {nsamples} samples')

    names.append(subdir)
    mutations.append(mat)"""


names = [x[len(FNAME) - 1:] for x in names]
names = ["COAD" if x == "COADREAD" else x for x in names]

"""sorted_names = ["LGG", "GBM", "OV", "UCEC", "CESC", "HNSC", "LUSC", "ESCA", "STAD", "BLCA", "BRCA", "LUAD", "LIHC", "SKCM", "SARC", "ACC", "COAD", "PAAD", "AML", "THCA", "KIRP", "KIRC", "PCPG", "THYM", "TEST"]
for i in range(0, len(names)):
    sorted_names[i] = names.index(sorted_names[i])

new_mutations = []
new_names = []
for i in range(0, len(names)):
    new_mutations.append(mutations[sorted_names[i]])
    new_names.append(names[sorted_names[i]])
mutations = new_mutations
names = new_names"""

"""new_mutations = []
for i in range(0, len(mutations)):
    if names[len(mutations) - 1 - i] in CANCERS:
        new_mutations.insert(0, mutations[len(mutations) - 1 - i])
    else:
        del names[len(mutations) - 1 - i]
mutations = new_mutations"""

"""for i in range(10, len(mutations)):
    del names[10]
    del mutations[10]"""

cumulative_size = np.zeros(len(mutations)+1).astype('int64')
for i in range(1, len(mutations) + 1):
    cumulative_size[i] = cumulative_size[i - 1] + np.size(mutations[i - 1], axis=1)
print(cumulative_size)

patients = np.concatenate(tuple(mutations), axis=1).T
# patients = patients - np.mean(patients, axis=0)


"""c = 0.0625 * 0.125
importance = np.load("../toe/importance.npy")
limit = np.quantile(importance, 1 - c)
print(limit)
print(np.min(importance))
print(np.shape(patients))
to_delete = []
for i in range(0, np.size(patients, axis=1)):
    if importance[i] < limit:
        to_delete.append(i)
patients = np.delete(patients, to_delete, 1)
input(np.shape(patients))"""

## Option 1: Pure PCA
# pca = PCA(n_components=2)
#shrunk_patients = pca.fit(patients).transform(patients)

## Option 2: PCA+tSNE
shrunk_patients = patients
if np.size(patients, axis=1) > 100:
    pca = PCA(n_components=100)
    # shrunk_patients = patients
    shrunk_patients = pca.fit_transform(patients)
    print(shrunk_patients)

# np.random.shuffle(shrunk_patients)
#manifold = TSNE(n_components=2,perplexity=90,random_state=42)
tsne = TSNE(
    perplexity=30,
    metric="euclidean",
    n_jobs=8,
    random_state=42,
    n_iter=250
)

shrunk_patients = tsne.fit(shrunk_patients)
#manifold = MDS()
#shrunk_patients = manifold.fit_transform(shrunk_patients)

print(np.shape(shrunk_patients))

np.save("patients_tsne-interim-cna.npy", shrunk_patients)


"""num_mutated = np.sum(np.abs(patients[:, :763]), axis=1)
num_mutated = np.log(num_mutated + 1)
num_mutated = num_mutated / np.quantile(num_mutated, 0.9)
#num_mutated = num_mutated / 10

color = []
for i in range(0, np.size(num_mutated, axis=0)):
    color.append(plt.get_cmap('viridis')(int(min([254, num_mutated[i] * 255]))))"""

fig, axs = plt.subplots(figsize=(5,5))
print(cumulative_size)
print(names)
for i in range(0, len(mutations)):
    plt.get_cmap('tab10')(groupsep(names[i]))
    begin = cumulative_size[i]
    end = cumulative_size[i + 1]
    axs.scatter(shrunk_patients[begin:end, 0], shrunk_patients[begin:end, 1], s=1, label=names[i], c=plt.get_cmap("tab10")(i % 10))
    #for j in range(cumulative_size[i], cumulative_size[i + 1]):
    #    axs.scatter(shrunk_patients[j, 0], shrunk_patients[j, 1], s=1, label=names[i], c=plt.get_cmap('viridis')(min([254, num_mutated[j] * 10])))
    #for j in range(cumulative_size[i],cumulative_size[i + 1]):
    #    axs.text(shrunk_patients[j, 0] * (1 + 0.001), shrunk_patients[j, 1] * (1 + 0.001) , names[i], fontsize=4)
#axs.scatter(shrunk_patients[:, 0], shrunk_patients[:, 1], s=1, c=color)
#, c=plt.get_cmap('viridis')(np.min(254, num_mutated * 255))

"""for i in range(0, np.size(shrunk_patients, axis=0)):
    if -33 < shrunk_patients[i, 0] and shrunk_patients[i, 0] <= -28:
        if 6 < shrunk_patients[i, 1] and shrunk_patients[i, 1] <= 21:
            strr = str("%04d" % i) + ":"
            for j in range(0, np.size(patients, axis=1)):
                if patients[i, j] == 0:
                    strr = strr + " "
                else:
                    #strr = strr + str(int(patients[i, j]))
                    strr = strr + "1"
                #if j % 150 == 0:
                #    strr = strr + "\n"
                #    strr = strr + "     "
            strr = strr + "| "
            for j in range(0, 2):
                strr = strr + str(np.around(shrunk_patients[i, j], 2)) + " "
            axs.scatter(shrunk_patients[i, 0], shrunk_patients[i, 1], s=1, c='k')
            print(strr)"""

#for i in range(0, len(patients)):
#    if num_mutated[i] == 1:
#        axs.scatter(shrunk_patients[i, 0], shrunk_patients[i, 1], s=1, c='k')
axs.legend()
plt.show()
