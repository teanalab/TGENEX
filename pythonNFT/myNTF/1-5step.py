def nnlsm_blockpivot(A, B, isInputProd=0, init=None):
    """
    Nonnegativity Constrained Least Squares with Multiple Righthand Sides
         using Block Principal Pivoting method

    This function solves the following problem: given A and B, find X such that
               minimize || AX-B ||_F^2 where X>=0 elementwise.

    Reference:
        Jingu Kim and Haesun Park. Fast Nonnegative Matrix Factorization:
            An Activeset-like Method and Comparisons.
            SIAM Journal on Scientific Computing, 33(6), pp. 3261-3281, 2011.

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
    MAX_BIG_ITER = n * 5

    # set initial feasible solution
    X = zeros((n, k))
    if init is None:
        Y = - AtB
        PassiveSet = zeros((n, k), dtype=np.bool)
    else:
        PassiveSet = (init > 0).copy()
        X = normalEqComb(AtA, AtB, PassiveSet)
        Y = AtA.dot(X) - AtB

    # parameters
    pbar = 3
    P = zeros((1, k))
    P[:] = pbar
    Ninf = zeros((1, k))
    Ninf[:] = n + 1

    NonOptSet = (Y < 0) & ~PassiveSet
    InfeaSet = (X < 0) & PassiveSet
    NotGood = (np.sum(NonOptSet, axis=0) + \
               np.sum(InfeaSet, axis=0))[np.newaxis, :]
    NotOptCols = NotGood > 0

    bigIter = 0

    while find(NotOptCols).shape[0] > 0:

        bigIter = bigIter + 1
        # set max_iter for ill-conditioned (numerically unstable) case
        if ((MAX_BIG_ITER > 0) & (bigIter > MAX_BIG_ITER)):
            break

        Cols1 = NotOptCols & (NotGood < Ninf)
        Cols2 = NotOptCols & (NotGood >= Ninf) & (P >= 1)
        Cols3Ix = find(NotOptCols & ~Cols1 & ~Cols2)

        if find(Cols1).shape[0] > 0:
            P[Cols1] = pbar
            NotGood[Cols1]
            Ninf[Cols1] = NotGood[Cols1]
            PassiveSet[NonOptSet & tile(Cols1, (n, 1))] = True
            PassiveSet[InfeaSet & tile(Cols1, (n, 1))] = False

        if find(Cols2).shape[0] > 0:
            P[Cols2] = P[Cols2] - 1
            PassiveSet[NonOptSet & tile(Cols2, (n, 1))] = True
            PassiveSet[InfeaSet & tile(Cols2, (n, 1))] = False

        if Cols3Ix.shape[0] > 0:
            for i in range(Cols3Ix.shape[0]):
                Ix = Cols3Ix[i]
                toChange = np.max(find(NonOptSet[:, Ix] | InfeaSet[:, Ix]))
                if PassiveSet[toChange, Ix]:
                    PassiveSet[toChange, Ix] = False
                else:
                    PassiveSet[toChange, Ix] = True

        Z = normalEqComb(AtA, AtB[:, NotOptCols.flatten()],
                         PassiveSet[:, NotOptCols.flatten()])
        X[:, NotOptCols.flatten()] = Z[:]
        X[abs(X) < 1e-12] = 0  # for numerical stability.
        Y[:, NotOptCols.flatten()] = AtA.dot(X[:, NotOptCols.flatten()]) - \
                                     AtB[:, NotOptCols.flatten()]
        Y[abs(Y) < 1e-12] = 0  # for numerical stability.

        # check optimality
        NotOptMask = tile(NotOptCols, (n, 1))
        NonOptSet = NotOptMask & (Y < 0) & ~PassiveSet
        InfeaSet = NotOptMask & (X < 0) & PassiveSet
        NotGood = (np.sum(NonOptSet, axis=0) +
                   np.sum(InfeaSet, axis=0))[np.newaxis, :]
        NotOptCols = NotGood > 0

    return X, Y


def getGradient(X, F, nWay, r):
    grad = []
    for k in range(nWay):
        ways = range(nWay)
        ways.remove(k)
        XF = X.uttkrp(F, k)
        # Compute the inner-product matrix
        FF = ones((r, r))
        for i in ways:
            FF = FF * (F[i].T.dot(F[i]))
        grad.append(F[k].dot(FF) - XF)
    return grad


def getProjGradient(X, F, nWay, r):
    pGrad = []
    for k in range(nWay):
        ways = range(nWay)
        ways.remove(k)
        XF = X.uttkrp(F, k)
        # Compute the inner-product matrix
        FF = ones((r, r))
        for i in ways:
            FF = FF * (F[i].T.dot(F[i]))
        grad = F[k].dot(FF) - XF
        grad[~((grad < 0) | (F[k] > 0))] = 0.
        pGrad.append(grad)
    return pGrad


class anls_asgroup(object):

    def initializer(self, X, F, nWay, orderWays):
        F[orderWays[0]] = zeros(F[orderWays[0]].shape)
        FF = []
        for k in range(nWay):
            FF.append((F[k].T.dot(F[k])))
        return F, FF

    def iterSolver(self, X, F, FF_init, nWay, r, orderWays):
        # solve NNLS problems for each factor
        for k in range(nWay):
            curWay = orderWays[k]
            ways = range(nWay)
            ways.remove(curWay)
            XF = X.uttkrp(F, curWay)
            # Compute the inner-product matrix
            FF = ones((r, r))
            for i in ways:
                FF = FF * FF_init[i]  # (F[i].T.dot(F[i]))
            ow = 0
            Fthis, temp = nnlsm_activeset(FF, XF.T, ow, 1, F[curWay].T)
            F[curWay] = Fthis.T
            FF_init[curWay] = (F[curWay].T.dot(F[curWay]))
        return F, FF_init
