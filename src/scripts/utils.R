##########################################################################################
#
# gCIS pipeline utils
# Dec 3, 2022
# Anders E. & Patryk S.
# Taylor lab
# 
##########################################################################################



# Notes -------------------------------------------------------------------
# all credit to Patryk Skowron for writing the original versions of all these functions
# which are in some cases entirely intact here, e.g.
# cluster_insertions(), most of neg_binom_thresh(), etc etc.

# Palette -----------------------------------------------------------------
pal <- c(blue = "#345b94", purple = "#902A7A", red = "#d2403b",
         yellow = "#F1E902", green = "#1D9849", orange = "#F68B1F", grey = "#ADAFB2")
tableau10 <- c("#4E79A7", "#F28E2B", "#E15759", "#76B7B2", "#59A14F", "#EDC948", "#B07AA1", "#FF9DA7", "#9C755F", "#BAB0AC") # from {ggthemes}
tableau20 <- c('#4E79A7', '#A0CBE8', '#F28E2B', '#FFBE7D', '#59A14F', '#8CD17D', '#B6992D', '#F1CE63', '#499894', '#86BCB6', '#E15759', '#FF9D9A', '#79706E', '#BAB0AC', '#D37295', '#FABFD2', '#B07AA1', '#D4A6C8', '#9D7660', '#D7B5A6')
phase_pal <- c("#2B83BA", "#FDAE61", "#D7191C")

# Misc functions -----------------------------------------------------------------
`%notin%` <- Negate(`%in%`)

#-------------------------------------------
##    CLUSTER INSERTIONS 
#-------------------------------------------
# Assumes that:
# (1) the input is sorted by chr and location 
# (2) either the plus or minus strand (not both) 
# Goes through each insertion and creates a new cluster if the next insertion is more than 'window' 
# O(n) Linear Time

cluster_insertions <- function(tbl, window){
  
  tbl$cluster[1]=1
  if(nrow(tbl) == 1){return(tbl)} #Edge case where there is only one insertion found in the sample
  
  for (row in 2:nrow(tbl)){
    if(tbl$loc[row] - tbl$loc[row - 1] < window &  tbl$chr[row] == tbl$chr[row - 1])
    { tbl$cluster[row] = tbl$cluster[row -1]} else {tbl$cluster[row] = tbl$cluster[row -1] + 1}
  }
  return(tbl)
}

# for use within a dtplyr pipeline:
cluster_insertions_group <- function(tbl, groups){
  window <- 5
  
  tbl$cluster[1]=1
  if(nrow(tbl) == 1){return(tbl)} #Edge case where there is only one insertion found in the sample
  
  for (row in 2:nrow(tbl)){
    if(tbl$loc[row] - tbl$loc[row - 1] < window &  tbl$chr[row] == tbl$chr[row - 1])
    { tbl$cluster[row] = tbl$cluster[row -1]} else {tbl$cluster[row] = tbl$cluster[row -1] + 1}
  }
  return(tbl)
}

#-------------------------------------------
##    COLLAPSE INSERTIONS 
#-------------------------------------------
# Insertions with the same cluster are collapsed into a single insertion
# This new insertion will have a new count = sum of counts in cluster
# This new insertion will have a new ULP = sum of ULP in cluster

collapse_clust <-function(tbl, type="shear"){
  
  if(type=="shear"){
    #Sort by ULP then count - pick this location for the 'true' insertion
    out = tbl %>% arrange(desc(ULP), desc(count)) %>% head(. , n=1)
    
    #Sum up the ULP in cluster and set as new ULP
    out$ULP=sum(tbl$ULP)
    
    #Sum up the read counts in cluster and set as new count
    out$count=sum(tbl$count)
    
  } else if (type=="restriction"){
    
    out = tbl %>% arrange(desc(count)) %>% head(. , n=1)
    
    #Sum up the read counts in cluster and set as new count
    out$count=sum(tbl$count)
  }
  
  return(out)
}  

#-------------------------------------------
##    CALCULATE CLONALITY
#-------------------------------------------
# Within each library calculate the clonal by dividing each insertions ULP by the highest ULP in the library

calc_clonality <-function(tbl){
  
  highest_clonality = max(tbl$ULP)
  
  tbl$clonality = tbl$ULP/highest_clonality
  
  return(tbl)
}

#-------------------------------------------
##    MERGE IRL AND IRR LIBRARIES
#-------------------------------------------
#Sort the IRL and IRR insertions based on ULP and count (in cases where both orientations are found)
#Pick the highest ULP (if tied pick highest read count) orientation


merged_counts <- function(tbl, type="shear"){
  
  if (type == "shear"){
    
    out = tbl %>% arrange(desc(ULP), desc(count)) %>% head(. , n=1) %>% select(-orientation, -Sample_key, -cluster)
    
  } else if (type == "restriction"){
    
    out = tbl %>% arrange(desc(count)) %>% head(. , n=1) %>% select(-orientation, -Sample_key, -cluster)
    
  }
  
  return(out)
  
}

#-------------------------------------------
##    Top 1 percent read count filter
#-------------------------------------------
top_1per_max_insertion = function(tbl)
{
  threshold = max(tbl$count)*0.01
  
  return(tbl %>% mutate(top_1per_max_insertion=threshold))
  
}

#-------------------------------------------
##    Top 0.1 percent read sum
#-------------------------------------------
sum_reads_0.1_per = function(tbl){
  
  threshold = sum(tbl$count)*0.001
  
  return(tbl %>% mutate(sum_reads_0.1_per=threshold))
  
}


#-------------------------------------------
##    95% percentile of negative binomial
#-------------------------------------------
neg_binom_thresh <-function(tbl){
  
  #print(paste(tbl$patient, tbl$sample, sep="_"))
  
  #Calculate the 95th percentile of a negative binomial model fitted on the background insertions (counts 1-3)
  input = tbl %>% arrange(desc(count)) %>% group_by(count) %>% summarise(numer_of_sites=n()) %>% dplyr::select(numer_of_sites, count) %>% t(.) %>% as.matrix(.)
  
  count=sum(input[1,])
  reads=sum(input[2,])
  
  actual=input[1,1:3]
  size=.01
  prob=.001
  bestprob=prob
  bestsize=size
  points=input[2,1:3]-1
  bestsum=""
  while(size<=1){
    prob=.001
    while(prob<=1){
      sim=dnbinom(points, size, prob)
      sim=sim*count
      sum=sumsquared(actual, sim)
      if(bestsum==""){
        bestsum=sum
      }
      if(sum<bestsum){
        bestsum=sum
        bestsize=size
        bestprob=prob
      }
      prob=prob+.001
    }
    size=size+.01
  }
  
  threshold=qnbinom(0.99, bestsize, bestprob) + 1
  
  #print(threshold)
  
  #Use the threshold to flag the data
  return(tbl %>% mutate(neg_binom_thresh=threshold))
  
}

# for use in a dtplyr pipeline
neg_binom_thresh_groups <- function(tbl, groups){
  
  #print(paste(tbl$patient, tbl$sample, sep="_"))
  
  #Calculate the 95th percentile of a negative binomial model fitted on the background insertions (counts 1-3)
  input = tbl %>% as_tibble() %>% arrange(desc(count)) %>% group_by(count) %>% summarise(numer_of_sites=n()) %>% dplyr::select(numer_of_sites, count) %>% t(.) %>% as.matrix(.)
  
  count=sum(input[1,])
  reads=sum(input[2,])
  
  if(ncol(input) < 3){
    return(tbl %>% as_tibble() %>% mutate(neg_binom_thresh=0)) #need to clean this later
  }
  
  actual=input[1,1:3]  
  points=input[2,1:3]-1
  size=.01
  prob=.001
  bestprob=prob
  bestsize=size
  bestsum=""
  while(size<=1){
    prob=.001
    while(prob<=1){
      sim=dnbinom(points, size, prob)
      sim=sim*count
      sum=sumsquared(actual, sim)
      if(bestsum==""){
        bestsum=sum
      }
      if(sum<bestsum){
        bestsum=sum
        bestsize=size
        bestprob=prob
      }
      prob=prob+.001
    }
    size=size+.01
  }
  
  threshold=qnbinom(0.99, bestsize, bestprob) + 1
  
  #print(threshold)
  
  #Use the threshold to flag the data
  return(tbl %>% as_tibble() %>% mutate(neg_binom_thresh=threshold))
  
}
#-------------------------------------------
##    Error between model and data 
#-------------------------------------------

sumsquared <- function(real, test){
  sum=0
  for( i in 1:length(real)){
    sum=sum+ (real[i]-test[i])*(real[i]-test[i])
  }
  return(sum)
}


#-------------------------------------------
##    Jaccard matrix
#-------------------------------------------

jaccard_matrix <- function(in_dt, ref_dt, mice){
  require(foreach)
  
  # initialize matrix
  mtx <- as.data.frame(matrix(NA,length(mice),length(mice)))
  rownames(mtx) <- mice
  colnames(mtx) <- mice
  
  # nested loop
  foreach(i = mice) %do% {
    foreach(j = mice) %do% {
      mtx[i,j] <- length(intersect((in_dt %>% filter.(patient == i) %>% pull.(insertion_map)), 
                                   (ref_dt %>% filter.(patient == j) %>% pull.(insertion_map)))) / 
        length(unique(c((in_dt %>% filter.(patient == i) %>% pull.(insertion_map)), 
                        (ref_dt %>% filter.(patient == j) %>% pull.(insertion_map)))))
    }
  }
  
  return(mtx)
}
