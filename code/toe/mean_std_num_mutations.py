import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

DATA_DIR = '../../tcga/data'
FNAME = 'sample_matrix.txt'

mutations = []
names = []
for subdir, dirs, files in os.walk(DATA_DIR):
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
    mutations.append(mat)

names = [x[len(FNAME) - 1:] for x in names]
names = ["COAD" if x == "COADREAD" else x for x in names]

fig, axs = plt.subplots(2, 1, figsize=(25,6))
fig.tight_layout()

# Sort cancer types by average # of mutations.
average_mutations = [np.mean(x) for x in mutations]
names_sorted = [x for _, x in sorted(zip(average_mutations, names), key=lambda pair: pair[0])] # thanks stack overflow!
average_mutations = sorted(average_mutations)

axs[0].bar(names_sorted, average_mutations)
axs[0].title.set_text("Cancers sorted by p(mutated)")

# Sort cancer types by variability in # of mutations

variance_mutations = [np.sum(x, axis=0) for x in mutations]
variance_mutations = [np.std(x) / np.mean(x) for x in variance_mutations]
names_sorted = [x for _, x in sorted(zip(variance_mutations, names), key=lambda pair: pair[0])]
variance_mutations = sorted(variance_mutations)

axs[1].bar(names_sorted, variance_mutations)
axs[1].title.set_text("Cancers sorted by variability in # of mutations")

plt.show()
