from numpy.random import rand

# -----------------------------------------------
# Creating a synthetic 3rd-order tensor
# -----------------------------------------------
N1 = 20
N2 = 25
N3 = 30

R = 10

# Random initialization
np.random.seed(42)

A_org = np.random.rand(N1, R)
A_org[A_org < 0.4] = 0

B_org = rand(N2, R)
B_org[B_org < 0.4] = 0

C_org = rand(N3, R)
C_org[C_org < 0.4] = 0


X_ks = ktensor([A_org, B_org, C_org])
X = X_ks.totensor()


X_ks.lmbda



# -----------------------------------------------
# Tentative initial values
# -----------------------------------------------
A0 = np.random.rand(N1, R)
B0 = np.random.rand(N2, R)
C0 = np.random.rand(N3, R)
D0 = np.random.rand(N4, R)

Finit = [A0, B0, C0, D0]

# -----------------------------------------------
# Uncomment only one of the following
# -----------------------------------------------
X_approx_ks = nonnegative_tensor_factorization(X, R)

#     X_approx_ks = nonnegative_tensor_factorization(X, R,
#                                                    min_iter=5, max_iter=20)
#
#     X_approx_ks = nonnegative_tensor_factorization(X, R,
#                                                    method='anls_asgroup')
#
#     X_approx_ks = nonnegative_tensor_factorization(X, R,
#                                                    tol=1e-7, max_iter=300)
#
#     X_approx_ks = nonnegative_tensor_factorization(X, R,
#                                                    init=Finit)

# -----------------------------------------------
# Approximation Error
# -----------------------------------------------
X_approx = X_approx_ks.totensor()
X_err = (X - X_approx).norm() / X.norm()
print "Error:", X_err

