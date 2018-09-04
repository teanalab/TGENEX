from numpy import random
from nonnegfac.nmf import NMF
from nonnegfac.nmf import NMF_ANLS_BLOCKPIVOT
from nonnegfac.nmf import NMF_ANLS_AS_NUMPY
from nonnegfac.nmf import NMF_ANLS_AS_GROUP
from nonnegfac.nmf import NMF_HALS
from nonnegfac.nmf import NMF_MU

W_org = random.rand(300, 10)
H_org = random.rand(300, 10)
A = W_org.dot(H_org.T)
print '\nTesting NMF().run() ...\n'
W, H, info = NMF().run(A, 10)
