blas_lib='/Applications/MATLAB_R2018a.app/bin/maci64/libmwblas.dylib' 
mex('-DDEFINEUNIX','-largeArrayDims','mtimesx.c',blas_lib)
