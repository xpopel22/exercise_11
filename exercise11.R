#setwd("V:/MPA-BTB/MPA-PRG/exercise_11")

#library("Biostrings")


viterbi <- function(AAsequence, hmm){
  length_seq <- length(AAsequence)
  length_states <- length(hmm[["N"]])
  length_M <- length(hmm[["M"]])
  states <- rep("A",length_seq)
  tmpmat <- rep(0,length_seq*length_states)
  mat <- matrix(tmpmat,nrow=length_states,ncol=length_seq)
  print(mat)
  
  # mat[1,1] <- eA(X) + max(prob(A,znak)+prob(odkud))
  index_B <- which(hmm[["M"]] == as.character(AAsequence[1]))
  for (n in 1:length_states){
    mat[n,1] <- hmm[["B"]][n,index_B] + hmm[["pi"]][n]
  }
  
  index_B <- which(hmm[["M"]] == as.character(AAsequence[2]))
  maximum_prob <- mat[1,1]+hmm[["A"]][1,1]
  for (m in 2:length_states){
    maximum_prob <- max(maximum_prob, mat[m,1]+hmm[["A"]][m,1])
  }
  mat[1,2] <- hmm[["B"]][1,index_B] + maximum_prob
  
  cnt <- 0
  for (i in 1:length_seq){
    for (j in 1:length_states){
      cnt <- cnt + 1
    }
  }
  
  return_list <- list(
    prob_matrix = mat,
    hidden_states = states
  )
  return(return_list)
}


load("./HMM1.Rdata")
print(HMM1["M"])
hmm <- HMM1

AAsequence <- AAString("TGA")
viterbi(AAsequence, hmm)
