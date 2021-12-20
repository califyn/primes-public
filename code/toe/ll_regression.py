import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.linear_model import LogisticRegressionCV, LinearRegression, Lasso
from sklearn.multiclass import OneVsRestClassifier
from sklearn.metrics import confusion_matrix#, plot_confusion_matrix
from sklearn.model_selection import train_test_split

DATA_DIR = '../../tcga/data'


# make matrices as diagonal as possible
def blockify(mat, niter=10, verbose=True):
    n = mat.shape[0]
    glperm = np.asarray(range(n))
    perm = np.random.permutation(range(n))
    mat = mat[np.ix_(perm, perm)]
    glperm = glperm[perm]

    dmat = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            dmat[i, j] = 1.0 * abs(i - j)
    obj = lambda m: np.sum(m * (m >= 0) * dmat)
    for i in range(niter):
        if verbose:
            print(f'Block diagonalizing: {i+1}/{niter}')
        for k in np.random.permutation(range(n)):
            perms = []
            for s in range(n):
                perm = list(range(n))
                perm[k] = s
                perm[s] = k
                perms.append(perm)
                perm = list(range(n))
                xdel = perm[k]
                del perm[k]
                perm.insert(s, xdel)
                perms.append(perm)
            curobj = obj(mat)
            imps = []
            for perm in perms:
                mat = mat[np.ix_(perm, perm)]
                newobj = obj(mat)
                imps.append(curobj - newobj)
                invperm = np.argsort(perm)
                mat = mat[np.ix_(invperm, invperm)]
            imps = np.asarray(imps)
            if max(imps) > 0:
                s = np.argmax(imps)
                perm = perms[np.argmax(imps)]
                mat = mat[np.ix_(perm, perm)]
                glperm = glperm[perm]
    if verbose:
        print()

    return glperm

names = []
genes = []

cancers = ['']
mutations = []
FNAME = 'sample_matrix.txt'
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

# make 200 samples for each cancer
num_sample = 200
for i in range(0, len(mutations)):
    #print(np.shape(mutations[i]))
    while np.size(mutations[i], axis=1) < num_sample:
        mutations[i] = np.concatenate((mutations[i], mutations[i]), axis=1)
    rand_index = np.random.choice(np.size(mutations[i], axis=1), size=num_sample, replace=False)
    mutations[i] = mutations[i][:, rand_index]
    print(np.shape(mutations[i]))

X = np.concatenate(tuple(mutations), axis=1).T
y = []
for i in range(0, len(mutations)):
    for j in range(0, np.size(mutations[i], axis=1)):
        #y.append([1 if x == i else 0 for x in range(0, len(names))])
        y.append(i)
y = np.array(y)
print(y)
print(X.shape)
print(y.shape)

X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.33, random_state=20)

"""#reg = Lasso(alpha=0.5).fit(X, y)
reg = OneVsRestClassifier(LogisticRegression()).fit(X, y)
#print(reg.predict(X))

linear_results = []
for i in range(0, len(mutations)):
    predicted = reg.predict(mutations[i].T)
    linear_results.append(np.where(predicted == i, 0, 1))
#estimators = reg.estimators_
#linear_coefs = []
#for estimator in estimators:
#    linear_coefs.append(estimator.coef_)
linear_coefs = reg.coef_"""

# reg = OneVsRestClassifier(LogisticRegression(penalty='l1', C=1, solver='saga')).fit(X_train, y_train)
# reg = OneVsRestClassifier(LogisticRegression()).fit(X, y)


reg = LogisticRegressionCV(penalty='l1', solver='liblinear').fit(X_train, y_train)
#print(reg.predict(X))

"""logistic_results = []
logistic_results_full = []
results_true = []
for i in range(0, len(mutations)):
    predicted = reg.predict(mutations[i].T)
    #logistic_results.append(np.mean(predicted[:, i]))
    logistic_results.append(np.mean(np.where(predicted == i, 1, 0)))
    #logistic_results_full.append(np.sum(np.multiply(np.arange(0, len(names)) + 1,reg.predict(mutations[i].T)), axis=1))
    logistic_results_full.append(predicted)
    #logistic_results_full.append(np.sum(predicted, axis=0))
    results_true.append(np.full((np.size(mutations[i], axis=1)),i))
    #input(predicted[:50, :20].astype('int64'))
    #print(np.shape(logistic_results_full[i]))
    print(logistic_results_full[i])
logistic_coefs = []
#estimators = reg.estimators_
#top = []
#for estimator in estimators:
#   logistic_coefs.append(estimator.coef_)"""

predict = reg.predict(X_test)
#print(reg.coef_)
#input(np.shape(reg.coef_)) # n class n feature
#np.save("ll_coef.npy", reg.coef_)
#importance = np.max(reg.coef_ * np.array([np.sum(np.abs(mutations[i]), axis=1) for i in range(0, len(mutations))]), axis=0).T
#print(importance)
#input(np.shape(importance)) # n class n feature
#np.save("importance.npy", importance)

    #good_genes = np.where(np.abs(estimator.coef_) in np.sort(np.abs(estimator.coef_))[-10:])
    #print(np.sort(np.abs(estimator.coef_))[-10:])
    #input(good_genes)
    #for gene in good_genes:
    #    print(gene)
    #    print(genes[int(gene)])

#for i in range(0, len(mutations)):
#    print(names[i])
#    print(top[i])

"""select_genes = []
for i in range(0, len(mutations)):
    print(names[i])
    print("-----------------")
    list = []
    for j in range(0,len(genes)):
        if logistic_coefs[i][0][j] != 0:
            list.append((genes[j], logistic_coefs[i][0][j]))
    sort = sorted(list, key=lambda x: x[1] * x[1])[-5:]
    sort.reverse()
    for pair in sort:
        print(str(pair[1]) + ": " + pair[0])
        if pair[0] not in select_genes:
            select_genes.append(pair[0])
print(select_genes)
input("continue...")"""


# Fig shows 1st classification acc
"""fig, ax = plt.subplots(figsize=(20,5))
x = np.arange(len(mutations))
width = 0.35

ax.bar(x - width/2, linear_results, width, label='Linear')
ax.bar(x + width/2, logistic_results, width, label='Logistic')
ax.set_xticks(x)
ax.set_xticklabels(names)
ax.set_ylabel("Proportion 1st guess correct")
ax.legend()

ax.hlines(y=1/len(mutations), xmin=-width, xmax=np.max(x)+width, color='k', linestyle='--')
plt.show()"""

# Fig shows 1st classification acc, one method only
"""fig, ax = plt.subplots(figsize=(20,5))
x = np.arange(len(mutations))
width = 0.5

ax.bar(x, logistic_results, width, label='Logistic')

ax.hlines(y=1/len(mutations), xmin=-width, xmax=np.max(x)+width, color='k', linestyle='--')
plt.show()"""

# Fig shows histogram of coefficients
"""fig, axs = plt.subplots(nrows=4, ncols=6)
num_bars = 10
coef_histogram = []
for coef in logistic_coefs:
    coef_histogram.append(np.histogram(coef)[0])
for i in range(0, 4):
    for j in range(0, 6):
        print(np.shape(coef_histogram[4 * i + 6]))
        axs[i][j].bar(np.arange(0, num_bars), coef_histogram[4 * i + j])
plt.show()"""

# Confusion matrix
fig, ax = plt.subplots(figsize=(10,10))
#y_true = np.concatenate(tuple(results_true))
#y_pred = np.concatenate(tuple(logistic_results_full))
y_true = y_test
y_pred = predict
print(y_true[:50])
print(y_pred[:50])
confusion = confusion_matrix(y_true, y_pred, normalize='true')
#input(confusion)

perm = blockify(confusion)
confusion = confusion[np.ix_(perm, perm)]
new_names = []
print(perm)
for i in range(0, len(names)):
    new_names.append(names[int(perm[i])])
names = new_names

#confusion = confusion / np.sum(confusion, axis=0)


norm = plt.Normalize(confusion.min()-1, confusion.max()+1)
colours = plt.cm.afmhot(confusion)

ax.axis('off')

plt.subplots_adjust(left=0.2, top=0.5)

tbl = ax.table(cellText = np.around(confusion, 2), cellColours = colours, rowLabels=names, colLabels=names, loc='top')
tbl.set_fontsize(20)

# plot_confusion_matrix(reg, X, y, labels=names, normalize='true')
print(reg.C_)
plt.show()
# Alternative plot: Makeup up of indices from 1 to 25 (don't think this works very well)
"""
separate_labels = []
for i in range(0, len(mutations)):
    to_add = []
    for j in range(0, np.size(mutations[i], axis=1)):
        to_add.append([1 if x == i else 0 for x in range(0, len(names))])
    separate_labels.append(to_add)

results = np.zeros((len(mutations), len(mutations)))
for i in range(0, len(mutations)):
    predict = reg.predict_proba(mutations[i].T)
    #print(predict[:10])
    results[i] = np.bincount(np.argsort(predict, axis=1)[:, i], minlength=len(mutations))[::-1] / np.size(predict, axis=0)
print(results[0])

fig, ax = plt.subplots(figsize=(20,5))
results = results.T
print(np.sum(results[0, :]))
partial_sum = np.zeros((len(mutations)))
for result in results:
    ax.bar(np.arange(0, len(mutations)), result,  bottom=partial_sum)
    partial_sum = partial_sum + result
plt.show()"""
