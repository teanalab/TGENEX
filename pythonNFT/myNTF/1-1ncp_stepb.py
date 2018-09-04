#  Copyright (C) 2013 Istituto per l'Interscambio Scientifico I.S.I.
#  You can contact us by email (isi@isi.it) or write to:
#  ISI Foundation, Via Alassio 11/c, 10126 Torino, Italy.
#
#  This work is licensed under a Creative Commons 4.0
#  Attribution-NonCommercial-ShareAlike License
#  You may obtain a copy of the License at
#  http://creativecommons.org/licenses/by-nc-sa/4.0/
#
#  This program was written by Andre Panisson <panisson@gmail.com> at
#  the Data Science Lab of the ISI Foundation,
#  and its development was partly supported by
#  the EU FET project MULTIPLEX (grant number 317532).
#
'''
Nonnegative Tensor Factorization, based on the Matlab source code
available at Jingu Kim's (jingu.kim@gmail.com) home page:

    https://github.com/kimjingu/nonnegfac-matlab
    https://github.com/kimjingu/nonnegfac-python

Requires the installation of Numpy and Scikit-Tensor
    (https://github.com/mnick/scikit-tensor).

For examples, see main() function.

This code comes with no guarantee or warranty of any kind.
Created on Nov 2013

@author: Andre Panisson
@contact: panisson@gmail.com
@organization: ISI Foundation, Torino, Italy
'''

import numpy as np
from numpy import zeros, ones, diff, kron, tile, any, all, linalg
import numpy.linalg as nla
import time
from sktensor import ktensor


def find(condition):
    "Return the indices where ravel(condition) is true"
    res, = np.nonzero(np.ravel(condition))
    return res
