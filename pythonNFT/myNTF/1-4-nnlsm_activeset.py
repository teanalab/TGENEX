def nnlsm_activeset(A, B, overwrite=0, isInputProd=0, init=None):
    """
    Nonnegativity Constrained Least Squares with Multiple Righthand Sides
         using Active Set method

    This function solves the following problem: given A and B, find X such that
               minimize || AX-B ||_F^2 where X>=0 elementwise.

    Reference:
         Charles L. Lawson and Richard J. Hanson,
               Solving Least Squares Problems,
               Society for Industrial and Applied Mathematics, 1995
         M. H. Van Benthem and M. R. Keenan,
               Fast Algorithm for the Solution of Large-scale
               Non-negativity-constrained Least Squares Problems,
               J. Chemometrics 2004; 18: 441-450

    Based on the Matlab version written by Jingu Kim (jingu.kim@gmail.com)
                  School of Computational Science and Engineering,
                  Georgia Institute of Technology

    Parameters
    ----------
    A : input matrix (m x n) (by default),
        or A'*A (n x n) if isInputProd==1

    B : input matrix (m x k) (by default),
        or A'*B (n x k) if isInputProd==1

    overwrite : (optional, default:0)
        if turned on, unconstrained least squares solution is computed
        in the beginning

    isInputProd : (optional, default:0)
        if turned on, use (A'*A,A'*B) as input instead of (A,B)

    init : (optional) initial value for X

    Returns
    -------
    X : the solution (n x k)

    Y : A'*A*X - A'*B where X is the solution (n x k)
    """

    if isInputProd:
        AtA = A
        AtB = B
    else:
        AtA = A.T.dot(A)
        AtB = A.T.dot(B)

    n, k = AtB.shape
    MAX_ITER = n * 5

    # set initial feasible solution
    if overwrite:
        X = normalEqComb(AtA, AtB)
        PassSet = (X > 0).copy()
        NotOptSet = any(X < 0)
    elif init is not None:
        X = init
        X[X < 0] = 0
        PassSet = (X > 0).copy()
        NotOptSet = ones((1, k), dtype=np.bool)
    else:
        X = zeros((n, k))
        PassSet = zeros((n, k), dtype=np.bool)
        NotOptSet = ones((1, k), dtype=np.bool)

    Y = zeros((n, k))
    if (~NotOptSet).any():
        Y[:, ~NotOptSet] = AtA.dot(X[:, ~NotOptSet]) - AtB[:, ~NotOptSet]
    NotOptCols = find(NotOptSet)

    bigIter = 0

    while NotOptCols.shape[0] > 0:
        bigIter = bigIter + 1
        # set max_iter for ill-conditioned (numerically unstable) case
        if ((MAX_ITER > 0) & (bigIter > MAX_ITER)):
            break

        Z = normalEqComb(AtA, AtB[:, NotOptCols], PassSet[:, NotOptCols])

        Z[abs(Z) < 1e-12] = 0  # for numerical stability.

        InfeaSubSet = Z < 0
        InfeaSubCols = find(any(InfeaSubSet, axis=0))
        FeaSubCols = find(all(~InfeaSubSet, axis=0))

        if InfeaSubCols.shape[0] > 0:               # for infeasible cols
            ZInfea = Z[:, InfeaSubCols]
            InfeaCols = NotOptCols[InfeaSubCols]

            Alpha = zeros((n, InfeaSubCols.shape[0]))
            Alpha[:] = np.inf

            ij = np.argwhere(InfeaSubSet[:, InfeaSubCols])
            i = ij[:, 0]
            j = ij[:, 1]

            InfeaSubIx = np.ravel_multi_index((i, j), Alpha.shape)
            if InfeaCols.shape[0] == 1:
                InfeaIx = np.ravel_multi_index((i,
                                                InfeaCols * ones((len(j), 1),
                                                                 dtype=int)),
                                               (n, k))
            else:
                InfeaIx = np.ravel_multi_index((i, InfeaCols[j]), (n, k))

            Alpha.ravel()[InfeaSubIx] = X.ravel()[InfeaIx] / \
                (X.ravel()[InfeaIx] - ZInfea.ravel()[InfeaSubIx])

            minVal, minIx = np.min(Alpha, axis=0), np.argmin(Alpha, axis=0)
            Alpha[:, :] = kron(ones((n, 1)), minVal)

            X[:, InfeaCols] = X[:, InfeaCols] + \
                              Alpha * (ZInfea - X[:, InfeaCols])

            IxToActive = np.ravel_multi_index((minIx, InfeaCols), (n, k))

            X.ravel()[IxToActive] = 0
            PassSet.ravel()[IxToActive] = False

        if FeaSubCols.shape[0] > 0:  # for feasible cols

            FeaCols = NotOptCols[FeaSubCols]
            X[:, FeaCols] = Z[:, FeaSubCols]
            Y[:, FeaCols] = AtA.dot(X[:, FeaCols]) - AtB[:, FeaCols]

            Y[abs(Y) < 1e-12] = 0   # for numerical stability.

            NotOptSubSet = (Y[:, FeaCols] < 0) & ~PassSet[:, FeaCols]

            NewOptCols = FeaCols[all(~NotOptSubSet, axis=0)]
            UpdateNotOptCols = FeaCols[any(NotOptSubSet, axis=0)]

            if UpdateNotOptCols.shape[0] > 0:
                minIx = np.argmin(Y[:, UpdateNotOptCols] * \
                                  ~PassSet[:, UpdateNotOptCols], axis=0)
                idx = np.ravel_multi_index((minIx, UpdateNotOptCols), (n, k))
                PassSet.ravel()[idx] = True

            NotOptSet.T[NewOptCols] = False
            NotOptCols = find(NotOptSet)

    return X, Y
