require(NMF)

runNMF <- function(Amatrix)
{
  for (rankNMF in c(3:80) )
  {
    res <- nmf(Amatrix, rankNMF, "snmf/l", seed = 123456)
    
    # get matrix W
    w <- basis(res)
    dim(w)
    
    # get matrix H
    h <- coef(res)
    dim(h)
    
    patientNMF_PxM <- w
    mutationNMF_PxM <- h
    
    save( patientNMF_PxM, mutationNMF_PxM,
          file=paste0("temp/NMF_",rankNMF,"_PxM.RData") )
  }
}
