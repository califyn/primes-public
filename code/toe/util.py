import numpy as np
import sklearn.cluster as skc


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


if __name__ == '__main__':
    mat = np.array([
        [1, 0, 1, 0],
        [0, 0.5, 0, 0.7],
        [0.9, 0, 0.8, 0],
        [0, 0.6, 0, 1]
    ])

    # Get permutation (you may need to use this to also
    # permute the matrix labels)
    perm = blockify(mat)
    # Permute rows and columns of `mat`
    block = mat[np.ix_(perm, perm)]
    
    print('Original:')
    print(mat)
    print('Permuted:')
    print(block)
