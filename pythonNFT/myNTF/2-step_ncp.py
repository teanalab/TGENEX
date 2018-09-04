
class anls_bpp(object):

    def initializer(self, X, F, nWay, orderWays):
        F[orderWays[0]] = zeros(F[orderWays[0]].shape)
        FF = []
        for k in range(nWay):
            FF.append((F[k].T.dot(F[k])))
        return F, FF

    def iterSolver(self, X, F, FF_init, nWay, r, orderWays):
        for k in range(nWay):
            curWay = orderWays[k]
            ways = range(nWay)
            ways.remove(curWay)
            XF = X.uttkrp(F, curWay)
            # Compute the inner-product matrix
            FF = ones((r, r))
            for i in ways:
                FF = FF * FF_init[i]  # (F[i].T.dot(F[i]))
            Fthis, temp = nnlsm_blockpivot(FF, XF.T, 1, F[curWay].T)
            F[curWay] = Fthis.T
            FF_init[curWay] = (F[curWay].T.dot(F[curWay]))
        return F, FF_init


def getStopCriterion(pGrad, nWay, nr_grad_all):
    retVal = np.sum(np.linalg.norm(pGrad[i], 'fro') ** 2
                    for i in range(nWay))
    return np.sqrt(retVal) / nr_grad_all


def getRelError(X, F_kten, nWay, nr_X):
    error = nr_X ** 2 + F_kten.norm() ** 2 - 2 * F_kten.innerprod(X)
    return np.sqrt(max(error, 0)) / nr_X

