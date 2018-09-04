
def nonnegative_tensor_factorization(X, r, method='anls_bpp',
                                     tol=1e-4, stop_criterion=1,
                                     min_iter=20, max_iter=200, max_time=1e6,
                                     init=None, orderWays=None):
    """
    Nonnegative Tensor Factorization (Canonical Decomposition / PARAFAC)

    Based on the Matlab version written by Jingu Kim (jingu.kim@gmail.com)
               School of Computational Science and Engineering,
               Georgia Institute of Technology

    This software implements nonnegativity-constrained low-rank approximation
    of tensors in PARAFAC model. Assuming that a k-way tensor X and target rank
    r are given, this software seeks F1, ... , Fk by solving the following
    problem:

    minimize
        || X- sum_(j=1)^r (F1_j o F2_j o ... o Fk_j) ||_F^2 +
              G(F1, ... , Fk) + H(F1, ..., Fk)
    where
        G(F1, ... , Fk) = sum_(i=1)^k ( alpha_i * ||Fi||_F^2 ),
        H(F1, ... , Fk) = sum_(i=1)^k ( beta_i sum_(j=1)^n || Fi_j ||_1^2 ).
    such that
        Fi >= 0 for all i.

    To use this software, it is necessary to first install scikit_tensor.

    Reference:
         Fast Nonnegative Tensor Factorization with an Active-set-like Method.
         Jingu Kim and Haesun Park.
         In High-Performance Scientific Computing: Algorithms and Applications,
         Springer, 2012, pp. 311-326.

    Parameters
    ----------
    X : tensor' object of scikit_tensor
        Input data tensor.

    r : int
        Target low-rank.

    method : string, optional
        Algorithm for solving NMF. One of the following values:
         'anls_bpp' 'anls_asgroup' 'hals' 'mu'
         See above paper (and references therein) for the details
         of these algorithms.
         Default is 'anls_bpp'.

    tol : float, optional
        Stopping tolerance. Default is 1e-4.
        If you want to obtain a more accurate solution,
        decrease TOL and increase MAX_ITER at the same time.

    min_iter : int, optional
        Minimum number of iterations. Default is 20.

    max_iter : int, optional
        Maximum number of iterations. Default is 200.

    init : A cell array that contains initial values for factors Fi.
            See examples to learn how to set.

    Returns
    -------
        F : a 'ktensor' object that represent a factorized form of a tensor.

    Examples
    --------
        F = nonnegative_tensor_factorization(X, 5)
        F = nonnegative_tensor_factorization(X, 10, tol=1e-3)
        F = nonnegative_tensor_factorization(X, 7, init=Finit, tol=1e-5)
    """

    nWay = len(X.shape)

    if orderWays is None:
        orderWays = np.arange(nWay)

    # set initial values
    if init is not None:
        F_cell = init
    else:
        Finit = [np.random.rand(X.shape[i], r) for i in range(nWay)]
        F_cell = Finit

    grad = getGradient(X, F_cell, nWay, r)

    nr_X = X.norm()
    nr_grad_all = np.sqrt(np.sum(np.linalg.norm(grad[i], 'fro') ** 2
                                 for i in range(nWay)))

    if method == "anls_bpp":
        method = anls_bpp()
    elif method == "anls_asgroup":
        method = anls_asgroup()
    else:
        raise Exception("Unknown method")

    # Execute initializer
    F_cell, FF_init = method.initializer(X, F_cell, nWay, orderWays)

    tStart = time.time()

    if stop_criterion == 2:
        F_kten = ktensor(F_cell)
        rel_Error = getRelError(X, ktensor(F_cell), nWay, nr_X)

    if stop_criterion == 1:
        pGrad = getProjGradient(X, F_cell, nWay, r)
        SC_PGRAD = getStopCriterion(pGrad, nWay, nr_grad_all)

    # main iterations
    for iteration in range(max_iter):
        cntu = True

        F_cell, FF_init = method.iterSolver(X, F_cell, FF_init, nWay, r, orderWays)
        F_kten = ktensor(F_cell)

        if iteration >= min_iter:

            if time.time() - tStart > max_time:
                cntu = False

            else:

                if stop_criterion == 1:
                    pGrad = getProjGradient(X, F_cell, nWay, r)
                    SC_PGRAD = getStopCriterion(pGrad, nWay, nr_grad_all)
                    if SC_PGRAD < tol:
                        cntu = False

                elif stop_criterion == 2:
                    prev_rel_Error = rel_Error
                    rel_Error = getRelError(X, F_kten, nWay, nr_X)
                    SC_DIFF = np.abs(prev_rel_Error - rel_Error)
                    if SC_DIFF < tol:
                        cntu = False
                else:
                    rel_Error = getRelError(X, F_kten, nWay, nr_X)
                    if rel_Error < 1:
                        cntu = False

        if not cntu:
            break

    return F_kten



