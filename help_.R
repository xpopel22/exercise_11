CountMutations <- function(vector) {
  # counting mutations 
  
  mutations <- 0
  
  # duplications
  num_dupl <- length(vector[duplicated(vector)])
  mutations <- mutations + num_dupl
  vector <- unique(vector)
  
  
  sorted_vec <- vector[order(vector)]
  indels <- c()
  inserts <- c()
  
  # inserts
  if (!(length(vector)==max(vector))) {
    for (i in length(vector):3) {
      if (sorted_vec[i] != sorted_vec[i-1]+1 & sorted_vec[i] != sorted_vec[i-2]+2) {
        inserts <- c(inserts,sorted_vec[i])
      }
    }
  }
  mutations <- mutations + length(inserts)
  
  # indels
  for (i in 1:length(vector)) {
    if (!(i %in% vector)){
      indels <- c(indels,i)
    }
  }
  mutations <- mutations + length(indels)
  
  # replacement
  where_is_insert <- which(inserts == vector)
  vector[where_is_insert] <- indels
  
  # reversals
  mutations <- mutations + breakpointSort(vector)
  return(mutations)
}


findSorted <- function(vector) {
  comparing_value = vector[1]
  for (i in (2:length(vector))) {
    if (vector[i] != comparing_value + 1) {
      return(i)
    }
    comparing_value = vector[i]
  }
  return(0)
}

indicateAscending <- function(vector){
  # ascending: 1
  # descending: 0
  ind_vector <- rep(0,length(vector))
  ind_vector[1] <- 1
  ind_vector[length(ind_vector)] <- 1
  for (i in (1:(length(vector)-1))){
    if (vector[i+1] == vector[i] + 1) {
      ind_vector[i] = 1 
      ind_vector[i+1] = 1
    }
  }
  return(ind_vector)
}

breakpointSort <- function(vector){
  # breakpoint method
  vector <- c(0,vector,max(vector)+1)
  i <- 0
  while (findSorted(vector) != 0) {
    parts <- indicateAscending(vector)
    inds_des <- which(parts %in% c(0))
    min_index <- inds_des[which.min(vector[inds_des])]
    first_breakpoint <- findSorted(vector)
    part_to_reverse <- c(first_breakpoint:min_index)
    
    if (length(part_to_reverse) == 1) {
      # finding last breakpoint
      for (j in (length(vector):1)) {
        if (vector[j] != vector[j-1] + 1) last_breakpoint <- j
        break
      }
      
      part_to_reverse <- c(first_breakpoint:last_breakpoint)
      vector[part_to_reverse] <- vector[last_breakpoint:first_breakpoint]
    }
    else vector[part_to_reverse] <- vector[min_index:first_breakpoint]
    i <- i +1
  }
  return(i)
} 

CountMutations(c(1,1,1,2,3,8))
