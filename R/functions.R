

tpm_approx <- function(data,length){

  RPK = data / (length/1000)
  TPM = RPK/(colSums(RPK)/1e6)
  return(TPM)
}

cpm_approx <- function(data){

  CPM = (data/colSums(data))* 1e6
  return(CPM)
}

depth_approx <- function(data){

  CPM = data / colSums(data)
  return(CPM)
}
