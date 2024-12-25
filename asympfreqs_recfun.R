#' We are computing the scaled equilibrium frequencies 
#' y_0/y_0,y_1/y_1,... (omitting the * from the ms)
#' of the non-core copies, denoted as I in the manuscript
#' this implements the recursion Eq. 8 from the manuscript
#' if we iterate it using y0=1, we already get yj/y0 when computing yj
#' Arguments:
#' until: Compute recursively ratios up to this number (included)
#' d: duplication rate
#' s: Selection coefficient, selective disadvantage per copy of gene
#' model either "ccm" or "stm" (anything but "ccm" goes to "stm", 
#' as currently coded)
#' io: Number of common (core) copies across the population

dupfreq_asymratio <- function(until=50,d=0.001,s=0.0005,
                              model1="ccm",io=1){
#' Summands of Eq. 8
if (model1=="ccm"){
#S <- 1; T <- 0  
duptrans_prob <- function(k,l){
  if (l %in% (k:(2*k))){
  out1 <- choose(k,l-k)*d^(l-k)*(1-d)^(2*k-l)} else {out1 <- 0}
  return(out1)}} else {
#S <- 0; T <- 1    
duptrans_prob <- function(k,l){
  if (l %in% c(k,k+1)){
  out1 <- d^(l-k)*(1-d)^(k-l+1)} else {out1 <- 0}
  return(out1)}  }


y <- c(1,rep(0,until))  #Iterate limit y, start in y_0=1
for (j in 1:until){
  #duplication also depends on io, and vector in R starts at 1
  #we just have -1 for i though since we compute y[j+1] 
  #in the loop
  real_dupprob <- function(i){duptrans_prob(io+i-1,io+j)}
  summands <- 1:j  #Summation indices of Eq. 5
  y[j+1] <- sum(y[summands]*sapply(summands,real_dupprob)*(1-s)^(0:(j-1)))
  y[j+1] <- y[j+1]/(duptrans_prob(io,io)-duptrans_prob(io+j,io+j)*(1-s)^j)
}
return(y)
}
