#change working directory
import os
import sys
print os.getcwd()
os.chdir('/Users/diamac/Documents/CLINGEN/9source/pythonNFT/myNTF/')

#read clinical and mutation data
import csv
import numpy as np

def iter_loadtxt(filename, num_rows, num_cols, delimiter=';', skiprows=0, skipcols=0, dtype=float):
    def iter_func():
        with open(filename, 'r') as infile:
            for _ in range(skiprows):
                next(infile)
            for line in infile:
                line = line.rstrip().split(delimiter)[skipcols:]
                for item in line:
                    yield dtype(item)
    data = np.fromiter(iter_func(), dtype=dtype, count=num_rows*num_cols)
    data = data.reshape((num_rows, num_cols))
    return data
    
#read matrices
clinical = iter_loadtxt('../data/10-binaC.csv', 456, 38, delimiter=';', skiprows=1, skipcols=0, dtype=int)
mutation = iter_loadtxt('../data/10-binaM.csv', 2714, 456, delimiter=';', skiprows=1, skipcols=0, dtype=int)


#tensor construction
from sktensor import dtensor, cp_als, ktensor

T = np.zeros((2714, 456, 38)).astype(int)
for i in range(0, 38):
    for j in range(0, 456):
        patient_i_j =  clinical[j,i]
        if patient_i_j == 1:
            T[:,j,i] = mutation[:,j]

np.shape(T)
dT = dtensor(T)
np.shape(dT)
rankT = 10
X_approx_ks = nonnegative_tensor_factorization(dT, rankT)


#https://github.com/mnick/scikit-tensor/blob/master/sktensor/ktensor.py
X_approx_ks.lmbda
X_approx_ks.U[0].shape #genes
genes = X_approx_ks.U[0]
sum(genes[:,0])
X_approx_ks.U[1].shape #patients
patients = X_approx_ks.U[1]
sum(patients[:,0])
X_approx_ks.U[2].shape #clinical
clinicalFactor = X_approx_ks.U[2]
sum(clinicalFactor[:,0])

folderName = '/Users/diamac/Documents/CLINGEN/9source/pythonNFT/temp/'

import numpy as np
np.savetxt(folderName+'clinicalMay.csv', X_approx_ks.U[2], delimiter=',')
np.savetxt(folderName+'genesMay.csv', X_approx_ks.U[0], delimiter=',')
np.savetxt(folderName+'patientsMay.csv', X_approx_ks.U[1], delimiter=',')


#https://stackoverflow.com/questions/4529815/saving-an-object-data-persistence
import pickle


#def writingMainFiles ():
#write files
with open('may2018_Tensor_GPC.pkl', 'wb') as output:
    pickle.dump(dT, output, pickle.HIGHEST_PROTOCOL)
    
with open('may2018_X_approx_ks.pkl', 'wb') as output:
    pickle.dump(X_approx_ks, output, pickle.HIGHEST_PROTOCOL)
        
with open('may2018_IntGPC.pkl', 'wb') as output:
    pickle.dump(T, output, pickle.HIGHEST_PROTOCOL)


#store all the environment
import shelve

filename='shelve.out'
my_shelf = shelve.open(filename,'n') # 'n' for new

for key in dir():
    try:
        my_shelf[key] = globals()[key]
    except TypeError:
        #
        # __builtins__, my_shelf, and imported modules can not be shelved.
        #
        print('ERROR shelving: {0}'.format(key))
my_shelf.close()
#To restore:

my_shelf = shelve.open(filename)
for key in my_shelf:
    globals()[key]=my_shelf[key]
my_shelf.close()

#Another way to save variables
import dill                            #pip install dill --user
filename = 'globalsave.pkl'
dill.dump_session(filename)

# and to load the session again:
dill.load_session(filename)