#setwd("V:/MPA-BTB/MPA-PRG/exercise_11")

library("Biostrings")


viterbi_old <- function(AAsequence, hmm){
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


viterbi <- function(AAsequence, hmm){
  
  T <- length(AAsequence)            # délka sekvence
  K <- length(hmm$N)                 # počet stavů
  
  # převod sekvence na char (DŮLEŽITÉ)
  obs <- as.character(AAsequence)
  obs <- unlist(strsplit(obs, split = ""))

  # mapování znaků na indexy emisí
  sym_index <- match(obs, hmm$M)
  print(sym_index)
  # Viterbi matice (log-pravděpodobnosti)
  V <- matrix(-Inf, nrow = K, ncol = T)
  
  # Backpointer
  back <- matrix(0, nrow = K, ncol = T)
  
  ## === Inicializace (t = 1) ===
  for (j in 1:K){
    V[j,1] <- hmm$pi[j] + hmm$B[j, sym_index[1]]
    back[j,1] <- 0
  }
  
  ## === Rekurze (t = 2..T) ===
  for (t in 2:T){
    for (j in 1:K){
      
      scores <- numeric(K)
      for (i in 1:K){
        scores[i] <- V[i,t-1] + hmm$A[i,j]
      }
      
      best_i <- which.max(scores)
      V[j,t] <- scores[best_i] + hmm$B[j, sym_index[t]]
      back[j,t] <- best_i
    }
  }
  
  ## === Terminace ===
  last_state <- which.max(V[,T])
  
  ## === Zpětný průchod ===
  states <- integer(T)
  states[T] <- last_state
  
  for (t in (T-1):1){
    states[t] <- back[states[t+1], t+1]
  }
  
  state_names <- hmm$N[states]
  
  return(list(
    log_prob = max(V[,T]),
    prob_matrix = V,
    hidden_states = state_names
  ))
}



load("./HMM1.Rdata")
load("./HMM2.Rdata")
load("./HMM3.Rdata")
print(HMM1["M"])
hmm <- HMM3

AAsequence <- AAString("AGTCTGG")
viterbi(AAsequence, hmm)
res <- viterbi(AAsequence, hmm)

res$hidden_states
res$log_prob
