#' ML grid search for CCM and STM based on the 
#' Poisson and negative binomial approximations
#' of the equilibria frequencies. These have essentially 2
#' parameters - i_o and d/s resp. d/(d+s). The latter parameters do 
#' only depend on the ratio of d to s. This means that we can just set 
#' s to 1 and d to how many times s it should be. 

run_human <- TRUE
run_fish <- FALSE
#' ML function
lml_f <- function(counts,d1,io,mod1){
                 out1 <- switch(mod1,
          "ccm"=sum(sapply(counts-io,dnbinom,
                            size=io,prob=1/(d1+1),log=TRUE)),
          "stm"=sum(sapply(counts-io,dpois,lambda=d1,log=TRUE))
          )
          return(out1)}

#' compute ML estimates for human gene families
if (run_human){
pops <- c("CEU","CHB","YRI")
for (pop1 in pops){
#adjust path in next line
count_fams <- read.csv(paste0("Copy_numbers-format-corrected-integers-transpose-",
                       pop1,".csv"),
                       sep="")

res_df <- data.frame("bestmodel"=rep("PO",ncol(count_fams)-1),
                     "lBFmodels"=rep(0,ncol(count_fams)-1),
                     "io"=rep(0,ncol(count_fams)-1),
                     "d"=rep(0,ncol(count_fams)-1),
                     "unique"=rep(FALSE,ncol(count_fams)-1))
rownames(res_df) <- colnames(count_fams)[-1]
for (g1 in 1:(ncol(count_fams)-1)){
  #' define parameter grid
  s <- 1
  d <- c(((400:4)/4)^(-1),seq(1.25,100,0.25))
  i_o <- 1:min(count_fams[,g1+1])
  mods <- c("ccm","stm")
  params1 <- expand.grid(d,i_o,mods)
  params1$Var3 <- as.character(params1$Var3)
  names(params1) <- c("d","io","model")
  
    cat(g1,"\n")
  ml_res <- sapply(1:nrow(params1),function(i){lml_f(counts = count_fams[,g1+1],
                                                    d1 = params1[i,1],
                                                    io = params1[i,2],
                                                    mod1 = params1[i,3])
  })
  lBF <- max(ml_res[params1$model=="ccm"]) - max(ml_res[params1$model=="stm"])
  max_pos <- which.max(ml_res)
  if (is.finite(lBF)){
    if (lBF > 1){res_df$bestmodel[g1] <- "ccm"} else {
      if (lBF < -1){res_df$bestmodel[g1] <- "stm"} else {res_df$bestmodel[g1] <- "both"}}
  } else {res_df$bestmodel[g1] <- ifelse(max(ml_res[params1$model=="ccm"])==0,
                                         "none","ccm")}
  res_df$lBFmodels[g1] <- lBF
  res_df$io[g1] <- params1$io[max_pos]
  res_df$d[g1] <- params1$d[max_pos]
  res_df$unique[g1] <- all(ml_res[-max_pos]<ml_res[max_pos])
}
write.table(res_df,file = paste0("copymodel_",pop1,"_fine.txt"))
}}
#' Zebrafish
if (run_fish){
  counts_dr <- list("FN_CGN"=list("genefam"="Fisnacht","pop"="CGN",
                                 counts=
                                   c(269,267,299,264,271,245,284,276)),
                   "B302_CGN"=list("genefam"="B302","pop"="CGN",
                                   counts=
                                     c(88,79,77,
                                       73,90,85,
                                       88,88)),
                   "FN_CHT"=list("genefam"="Fisnacht","pop"="CHT",
                                 counts=
                                   c(387,532,418,459,418,
                                     410,389,410,420,404,
                                     381,167,348,534,305,
                                     355,374,383)),
                   "B302_CHT"=list("genefam"="B302","pop"="CHT",
                                   counts=
                                     c(128,174,137,
                                       164,143,143,
                                       134,139,140,
                                       147,122,66,
                                       121,174,86,
                                       101,117,102)),
                   "FN_KG"=list("genefam"="Fisnacht","pop"="KG",
                                counts=
                                  c(485,444,518,382,402,
                                    408,477,351,431,429,
                                    448,445,464,454,496,
                                    282,231,223,280,365)),
                   "B302_KG"=list("genefam"="B302","pop"="KG",
                                  counts=
                                    c(159,161,180,
                                      137,143,137,
                                      162,118,141,
                                      147,164,154,
                                      161,160,180,
                                      75,60,61,
                                      73,129
                                    )),
                   "FN_SN"=list("genefam"="Fisnacht","pop"="SN",
                                counts=
                                  c(521,399,385,405,456,
                                    396,344,534,468,562,
                                    543,399,443,450,491,
                                    335,378,524,358)),
                   "B302_SN"=list("genefam"="B302","pop"="SN",
                                  counts=
                                    c(153,126,143,
                                      120,142,121,
                                      108,169,146,
                                      187,160,129,
                                      118,128,144,
                                      102,113,153,
                                      109
                                    )),
                   "FN_DP"=list("genefam"="Fisnacht","pop"="DP",
                                counts= c(61,278,336,
                                          285,286,294,
                                          258,250,298,
                                          435,279,282,
                                          261,252,328,
                                          289,260,309,
                                          235,201)),
                   "B302_DP"=list("genefam"="B302","pop"="DP",
                                  counts=c(12,57,75,
                                           67,60,60,
                                           71,51,65,
                                           134,61,61,
                                           55,63,79,
                                           61,48,76,
                                           55,46
                                  )),
                   "FN_TU"=list("genefam"="Fisnacht","pop"="TU",
                                counts=c(238,229,246,
                                         264,120,126,
                                         78,112)),
                   "B302_TU"=list("genefam"="B302","pop"="TU",
                                  counts=c(68,73,76,
                                           73,35,45,
                                           28,35
                                  ))
  )
counts_dr$"FN_wp" <- list("genefam"="Fisnacht","pop"="wp",
                          counts=c(counts_dr$FN_CHT$counts,
                                   counts_dr$FN_DP$counts,
                                   counts_dr$FN_KG$counts,
                                   counts_dr$FN_SN$counts))
counts_dr$"B302_wp" <- list("genefam"="B302","pop"="wp",
                          counts=c(counts_dr$B302_CHT$counts,
                                   counts_dr$B302_DP$counts,
                                   counts_dr$B302_KG$counts,
                                   counts_dr$B302_SN$counts))
 
res_df <- data.frame("bestmodel"=rep("PO",length(counts_dr)),
                     "lBFmodels"=rep(0,length(counts_dr)),
                     "io"=rep(0,length(counts_dr)),
                     "d"=rep(0,length(counts_dr)),
                     "unique"=rep(FALSE,length(counts_dr)))

rownames(res_df) <- names(counts_dr)

ml_list <- vector("list",length(counts_dr))
names(ml_list) <- names(counts_dr)

for (g1 in 1:length(counts_dr)){
  #Parameters
  s <- 1
  d <- c(((400:4)/4)^(-1),seq(1.25,100,0.25))
  i_o <- 1:min(counts_dr[[g1]][["counts"]])
  mods <- c("ccm","stm")
  params1 <- expand.grid(d,i_o,mods)
  params1$Var3 <- as.character(params1$Var3)
  names(params1) <- c("d","io","model")
  cat(g1,"\n")
  ml_res <- sapply(1:nrow(params1),function(i){lml_f(counts = counts_dr[[g1]][["counts"]],
                                                    d1 = params1[i,1],
                                                    io = params1[i,2],
                                                    mod1 = params1[i,3])
  })
  ml_list[[g1]] <- ml_res
  lBF <- max(ml_res[params1$model=="ccm"])-max(ml_res[params1$model=="stm"])
  max_pos <- which.max(ml_res)
  #To get CCM best estimate
  #max_ccm <- max(ml_res[params1$model=="ccm"])
  #params1[which(ml_res==max_ccm),]

  if (is.finite(lBF)){
    if (lBF > 1){res_df$bestmodel[g1] <- "ccm"} else {
      if (lBF < -1){res_df$bestmodel[g1] <- "stm"} else {res_df$bestmodel[g1] <- "both"}}
  } else {res_df$bestmodel[g1] <- ifelse(max(ml_res[params1$model=="ccm"])==0,
                                         "none","ccm")}
  res_df$lBFmodels[g1] <- lBF
  res_df$io[g1] <- params1$io[max_pos]
  res_df$d[g1] <- params1$d[max_pos]
  res_df$unique[g1] <- all(ml_res[-max_pos]<ml_res[max_pos])
}
write.table(res_df,file = "copymodel_dr_fine_recount.txt")
}