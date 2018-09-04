#read clinical and mutation data
import csv
import numpy as np
reader = csv.reader(open('../data/10-binaC.csv', 'r'), delimiter=";", quotechar='"')
x = list(reader)
clinical = np.array(x)
clinical[0,0]


#https://stackoverflow.com/questions/4529815/saving-an-object-data-persistence
import pickle

folderName = '/Users/diamac/Documents/CLINGEN/9source/pythonNFT/temp/'

#read output files
with open('Documents/current/CLINGEN/9source/pythonNFT/temp/feb2018_X_approx_ks.pkl', 'rb') as input:
    X_approx_ks = pickle.load(input)


#https://github.com/mnick/scikit-tensor/blob/master/sktensor/ktensor.py
X_approx_ks.lmbda
X_approx_ks.U[0].shape #clinical
X_approx_ks.U[1].shape #genes
X_approx_ks.U[2].shape #patients


#save a csv file for each factor matrix
#From Fedor suggestion
#https://github.com/teanalab/viz-and-rank/blob/master/tensorboard-save.py#L13
#I copied it from https://gist.github.com/seumasmorrison/4563724
#it's explained here http://stackoverflow.com/a/8964779/1135883
#you can try https://docs.scipy.org/doc/numpy/reference/generated/numpy.genfromtxt.html also
#but it will use more memory
import numpy as np
np.savetxt(folderName+'clinical.csv', X_approx_ks.U[0], delimiter=',')
np.savetxt(folderName+'genes.csv', X_approx_ks.U[1], delimiter=',')
np.savetxt(folderName+'patients.csv', X_approx_ks.U[2], delimiter=',')



def writingMainFiles ():
    #write files
    with open('may2018_Tensor_GPC.pkl', 'wb') as output:
        pickle.dump(Tensor_GPC, output, pickle.HIGHEST_PROTOCOL)
    
    with open('feb2018_X_approx_ks.pkl', 'wb') as output:
        pickle.dump(X_approx_ks, output, pickle.HIGHEST_PROTOCOL)
        
    with open('feb2018_IntGPC.pkl', 'wb') as output:
        pickle.dump(IntGPC, output, pickle.HIGHEST_PROTOCOL)
    
    
