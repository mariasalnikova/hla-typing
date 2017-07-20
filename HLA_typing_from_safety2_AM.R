#Download all nessesary packages
library(dplyr)
library(data.table)
library(cowsay)
library(Biostrings)
library(stringr)
library(igraph)
library(tibble)

get_HLA_result<-function (donor) {  say("START!")

    print(format(Sys.time(), "%a %b %d %X %Y"))
  
  if (sum(sapply(donor$safety2, nrow)[c(1,2,6,7,12,13,17,18)])>0)
    #if (donor$safety2!="fail")
    {
    safety3<-get_safety3(lapply(donor$safety2, function(x){x[((x$parents==1&x$freq>0.001)|(x$parents==2&x$freq>0.05))&!grepl(pattern = "NULL",x$assembled), ]}))
    safety4<-get_safety4(safety3)
    safety5<-get_safety5(DT1=safety4$DT1, DQB=safety4$DQB, DRB=safety4$DRB)
    list(safety3=safety3, safety4=safety4, safety5=safety5)
  }
}

get_safety3 <- function(safety2) {
  HLA_allele_list<-HLA_base
  HLA_allele_Iclass<-select(HLA_allele_list, Allele, Sequence, HLA_class)%>%filter(HLA_class=="A"|HLA_class=="B"|HLA_class=="C")%>%select(Allele, Sequence)
  HLA_allele_DRB<-select(HLA_allele_list, Allele, Sequence, HLA_class)%>%filter(grepl(HLA_allele_list$HLA_class, pattern="DRB*", fixed=F))%>%select(Allele, Sequence)
  HLA_allele_DQB<-select(HLA_allele_list, Allele, Sequence, HLA_class)%>%filter(grepl(HLA_allele_list$HLA_class, pattern="DQB*", fixed=F))%>%select(Allele, Sequence)
  
  N150<-paste0(rep(x = "N", times=30), collapse="")
  HLA_allele_Iclass$Sequence<-paste0(N150, HLA_allele_Iclass$Sequence, N150)
  HLA_allele_DRB$Sequence<-paste0(N150, HLA_allele_DRB$Sequence, N150)
  HLA_allele_DQB$Sequence<-paste0(N150, HLA_allele_DQB$Sequence, N150)
  
  HLA_allele_Iclass$Iamp1<-0
  HLA_allele_Iclass$Iamp1_inv<-0
  HLA_allele_Iclass$Iamp2<-0
  HLA_allele_Iclass$Iamp2_inv<-0
  HLA_allele_Iclass$Iamp1alt<-0
  HLA_allele_Iclass$Iamp1alt_inv<-0
  HLA_allele_Iclass$Iamp2alt<-0
  HLA_allele_Iclass$Iamp2alt_inv<-0
   
  HLA_allele_DQB$IIamp1_DQB<-0
  HLA_allele_DQB$IIamp1_DQB_inv<-0
  HLA_allele_DQB$IIamp2<-0
  HLA_allele_DQB$IIamp2_inv<-0
  HLA_allele_DQB$IIamp1_DQBalt<-0
  HLA_allele_DQB$IIamp1_DQBalt_inv<-0
  HLA_allele_DQB$IIamp2_DQBalt<-0
  HLA_allele_DQB$IIamp2_DQBalt_inv<-0
  
  HLA_allele_DRB$IIamp1_Others<-0
  HLA_allele_DRB$IIamp1_Others_inv<-0
  HLA_allele_DRB$IIamp2<-0
  HLA_allele_DRB$IIamp2_inv<-0
  HLA_allele_DRB$IIamp1_DRBalt<-0
  HLA_allele_DRB$IIamp1_DRBalt_inv<-0
  HLA_allele_DRB$IIamp2_DRBalt<-0
  HLA_allele_DRB$IIamp2_DRBalt_inv<-0
  
  HLA_Iclass_amps<-select(HLA_allele_Iclass, Allele)
  HLA_DQB_amps<-select(HLA_allele_DQB, Allele)
  HLA_DRB_amps<-select(HLA_allele_DRB, Allele)
  
  HLA_Iclass_amps$Iamp1<-0
  HLA_Iclass_amps$Iamp1_inv<-0
  HLA_Iclass_amps$Iamp2<-0
  HLA_Iclass_amps$Iamp2_inv<-0
  HLA_Iclass_amps$Iamp1alt<-0
  HLA_Iclass_amps$Iamp1alt_inv<-0
  HLA_Iclass_amps$Iamp2alt<-0
  HLA_Iclass_amps$Iamp2alt_inv<-0
  
  HLA_DQB_amps$IIamp1_DQB<-0
  HLA_DQB_amps$IIamp1_DQB_inv<-0
  HLA_DQB_amps$IIamp2<-0
  HLA_DQB_amps$IIamp2_inv<-0
  HLA_DQB_amps$IIamp1_DQBalt<-0
  HLA_DQB_amps$IIamp1_DQBalt_inv<-0
  HLA_DQB_amps$IIamp2_DQBalt<-0
  HLA_DQB_amps$IIamp2_DQBalt_inv<-0
  
  HLA_DRB_amps$IIamp1_Others<-0
  HLA_DRB_amps$IIamp1_Others_inv<-0
  HLA_DRB_amps$IIamp2<-0
  HLA_DRB_amps$IIamp2_inv<-0
  HLA_DRB_amps$IIamp1_DRBalt<-0
  HLA_DRB_amps$IIamp1_DRBalt_inv<-0
  HLA_DRB_amps$IIamp2_DRBalt<-0
  HLA_DRB_amps$IIamp2_DRBalt_inv<-0

  HLA1_dnastr<-DNAStringSet(HLA_allele_Iclass$Sequence)
  HLADQB_dnastr<- DNAStringSet(HLA_allele_DQB$Sequence)
  HLADRB_dnastr<- DNAStringSet(HLA_allele_DRB$Sequence)
  
    allele_freq_in_ampsI <- function(what,where, indices) {
    
    #what$assembled<-gsub(what$assembled, pattern = "N", replacement = ".")
    
    
    for (i in 1:nrow(what)) {
      # where[grepl(x=HLA_allele_Iclass$Sequence,
      #             pattern = what$assembled[i], 
      #             fixed=F)]<-where[grepl(x=HLA_allele_Iclass$Sequence, 
      #                                    pattern = what$assembled[i])]+what$freq[i]
      # 
      # indices[grepl(x=HLA_allele_Iclass$Sequence,
      #               pattern = what$assembled[i], 
      #               fixed=F)] <- paste0(indices[grepl(x=HLA_allele_Iclass$Sequence,
      #                                                 pattern = what$assembled[i], 
      #                                                 fixed=F)],",", i)
      VC<-vcountPattern(pattern = what$assembled[i], subject = HLA1_dnastr,
                        max.mismatch = 0, with.indels = F, fixed=F)
      
      where[VC!=0]<-where[VC!=0]+what$freq[i]

      indices[VC!=0] <- paste0(indices[VC!=0],",", i)
      
    }
    list(where,indices)
  }
  
  allele_freq_in_ampsDQB <- function(what,where, indices) {
    
    #what$assembled<-gsub(what$assembled, pattern = "N", replacement = ".")
    
    for (i in 1:nrow(what)) {
      # where[grepl(x=HLA_allele_DQB$Sequence,
      #             pattern = what$assembled[i], 
      #             fixed=F)]<-where[grepl(x=HLA_allele_DQB$Sequence, 
      #                                    pattern = what$assembled[i])]+what$freq[i]
      # 
      # indices[grepl(x=HLA_allele_DQB$Sequence,
      #               pattern = what$assembled[i], 
      #               fixed=F)] <- paste0(indices[grepl(x=HLA_allele_DQB$Sequence,
      #                                                 pattern = what$assembled[i], 
      #                                                 fixed=F)],",", i)
      
      VC<-vcountPattern(pattern = what$assembled[i], subject = HLADQB_dnastr,
                        max.mismatch = 0, with.indels = F, fixed=F)
      
      where[VC!=0]<-where[VC!=0]+what$freq[i]
      
      indices[VC!=0] <- paste0(indices[VC!=0],",", i)
    }
    list(where,indices)
  }
  
  allele_freq_in_ampsDRB <- function(what,where, indices) {
    
    #what$assembled<-gsub(what$assembled, pattern = "N", replacement = ".")
    
    for (i in 1:nrow(what)) {
      # where[grepl(x=HLA_allele_DRB$Sequence,
      #             pattern = what$assembled[i], 
      #             fixed=F)]<-where[grepl(x=HLA_allele_DRB$Sequence, 
      #                                    pattern = what$assembled[i])]+what$freq[i]
      # 
      # indices[grepl(x=HLA_allele_DRB$Sequence,
      #               pattern = what$assembled[i], 
      #               fixed=F)] <- paste0(indices[grepl(x=HLA_allele_DRB$Sequence,
      #                                                 pattern = what$assembled[i], 
      #                                                 fixed=F)],",", i)
      VC<-vcountPattern(pattern = what$assembled[i], subject = HLADRB_dnastr,
                        max.mismatch = 0, with.indels = F, fixed=F)
      
      where[VC!=0]<-where[VC!=0]+what$freq[i]
      
      indices[VC!=0] <- paste0(indices[VC!=0],",", i)
      
      }
    list(where,indices)
  }
  
  
  n_readsI<-rep(0,times=8)
  n_readsDQB<-rep(0,times=8)
  n_readsDRB<- rep(0,times=8)
  names(n_readsI)<-colnames(HLA_allele_Iclass)[3:10]
  names(n_readsDQB)<-colnames(HLA_allele_DQB)[3:10]
  names(n_readsDRB)<-colnames(HLA_allele_DRB)[3:10]
  
  for (colname in intersect(colnames(HLA_allele_Iclass)[3:10], names(safety2))) {
    if (nrow(safety2[[colname]])>0) 
      {
      tempres<-allele_freq_in_ampsI(safety2[[colname]], HLA_allele_Iclass[[colname]], HLA_Iclass_amps[[colname]])
      HLA_allele_Iclass[[colname]]<-tempres[[1]]
      HLA_Iclass_amps[[colname]]<- tempres[[2]]
      n_readsI[colname]<-sum(safety2[[colname]]$readnumber)
      }
  }
  for (colname in intersect(colnames(HLA_allele_DQB)[3:10], names(safety2))) {
    if (nrow(safety2[[colname]])>0) 
    {
      tempres<-allele_freq_in_ampsDQB(safety2[[colname]], HLA_allele_DQB[[colname]], HLA_DQB_amps[[colname]])
      HLA_allele_DQB[[colname]]<-tempres[[1]]
      HLA_DQB_amps[[colname]]<- tempres[[2]]
      n_readsDQB[colname]<-sum(safety2[[colname]]$readnumber)
    }
  }
  for (colname in intersect(colnames(HLA_allele_DRB)[3:10], names(safety2))) {
    if (nrow(safety2[[colname]])>0) 
      {
      tempres<-allele_freq_in_ampsDRB(safety2[[colname]], HLA_allele_DRB[[colname]], HLA_DRB_amps[[colname]])
      HLA_allele_DRB[[colname]]<-tempres[[1]]
      HLA_DRB_amps[[colname]]<- tempres[[2]]
      n_readsDRB[colname]<-sum(safety2[[colname]]$readnumber)
    }
  }
  
  safety3<-list(HLA_allele_Iclass=HLA_allele_Iclass, HLA_Iclass_amps=HLA_Iclass_amps, 
                HLA_allele_DQB=HLA_allele_DQB, HLA_DQB_amps=HLA_DQB_amps,
                HLA_allele_DRB=HLA_allele_DRB, HLA_DRB_amps=HLA_DRB_amps,
                n_readsI=n_readsI, n_readsDQB=n_readsDQB, n_readsDRB=n_readsDRB)
  safety3
}


get_safety4<-function(safety3) {
  safety3$HLA_allele_Iclass$amps<-apply(safety3$HLA_Iclass_amps[2:9], MARGIN = 1, paste0, collapse = "_")
  safety3$HLA_allele_DQB$amps<-apply(safety3$HLA_DQB_amps[2:9], MARGIN = 1, paste0, collapse = "_")
  safety3$HLA_allele_DRB$amps<-apply(safety3$HLA_DRB_amps[2:9], MARGIN = 1, paste0, collapse = "_")
  
  for (i in 1:nrow(safety3$HLA_allele_Iclass)) {
    safety3$HLA_allele_Iclass$no_amps[i]<-paste0(colnames(safety3$HLA_Iclass_amps[which(safety3$HLA_Iclass_amps[i, ]==0)]), collapse = "   ")
  }
  for (i in 1:nrow(safety3$HLA_allele_DQB)) {
    safety3$HLA_allele_DQB$no_amps[i]<-paste0(colnames(safety3$HLA_DQB_amps[which(safety3$HLA_DQB_amps[i, ]==0)]), collapse = "   ")
  }
  for (i in 1:nrow(safety3$HLA_allele_DRB)) {
    safety3$HLA_allele_DRB$no_amps[i]<-paste0(colnames(safety3$HLA_DRB_amps[which(safety3$HLA_DRB_amps[i, ]==0)]), collapse = "   ")
  }
  
  #safety3$HLA_allele_Iclass[3:10]<-apply(safety3$HLA_allele_Iclass[3:10], MARGIN = 2, prop.table)
  DT1<-as.data.table(safety3$HLA_allele_Iclass)
  DQB<-as.data.table(safety3$HLA_allele_DQB)
  DRB<-as.data.table(safety3$HLA_allele_DRB)
  DT1<-DT1[amps!="0_0_0_0_0_0_0_0", .(Allele=paste0(Allele, collapse=" "),
                                    Iamp1=unique(Iamp1),
                                    Iamp1_inv=unique(Iamp1_inv),
                                    Iamp2=unique(Iamp2),
                                    Iamp2_inv=unique(Iamp2_inv),
                                    Iamp1alt=unique(Iamp1alt),
                                    Iamp1alt_inv=unique(Iamp1alt_inv),
                                    Iamp2alt=unique(Iamp2alt),
                                    Iamp2alt_inv=unique(Iamp2alt_inv), 
                                    no_amps=unique(no_amps)),by=amps]
  
  DQB<-DQB[amps!="0_0_0_0_0_0_0_0", .(Allele=paste0(Allele, collapse=" "),
                                      IIamp1_DQB=unique(IIamp1_DQB),
                                      IIamp1_DQB_inv=unique(IIamp1_DQB_inv),
                                      IIamp2=unique(IIamp2),
                                      IIamp2_inv=unique(IIamp2_inv),
                                      IIamp1_DQBalt=unique(IIamp1_DQBalt),
                                      IIamp1_DQBalt_inv=unique(IIamp1_DQBalt_inv),
                                      IIamp2_DQBalt=unique(IIamp2_DQBalt),
                                      IIamp2_DQBalt_inv=unique(IIamp2_DQBalt_inv),
                                      no_amps=unique(no_amps)), by=amps]
  
  DRB<-DRB[amps!="0_0_0_0_0_0_0_0", .(Allele=paste0(Allele, collapse=" "),
                                      IIamp1_Others=unique(IIamp1_Others),
                                      IIamp1_Others_inv=unique(IIamp1_Others_inv),
                                      IIamp2=unique(IIamp2),
                                      IIamp2_inv=unique(IIamp2_inv),
                                      IIamp1_DRBalt=unique(IIamp1_DRBalt),
                                      IIamp1_DRBalt_inv=unique(IIamp1_DRBalt_inv),
                                      IIamp2_DRBalt=unique(IIamp2_DRBalt),
                                      IIamp2_DRBalt_inv=unique(IIamp2_DRBalt_inv),
                                      no_amps=unique(no_amps)), by=amps]
  DT1<-as.data.frame(DT1)
  DQB<-as.data.frame(DQB)
  DRB<-as.data.frame(DRB)
  
  DT1$sumreads<-apply(DT1[3:10], 1, function(x) {sum(x*safety3$n_readsI)})
  DQB$sumreads<-apply(DQB[3:10], 1, function(x) {sum(x*safety3$n_readsI)})
  DRB$sumreads<-apply(DRB[3:10], 1, function(x) {sum(x*safety3$n_readsI)})
  
  DT1$meanreads<-apply(DT1[3:10], 1, function(x) {mean(x*safety3$n_readsI)})
  DQB$meanreads<-apply(DQB[3:10], 1, function(x) {mean(x*safety3$n_readsI)})
  DRB$meanreads<-apply(DRB[3:10], 1, function(x) {mean(x*safety3$n_readsI)})
  
  DT1$medianreads<-apply(DT1[3:10], 1, function(x) {median(x*safety3$n_readsI)})
  DQB$medianreads<-apply(DQB[3:10], 1, function(x) {median(x*safety3$n_readsI)})
  DRB$medianreads<-apply(DRB[3:10], 1, function(x) {median(x*safety3$n_readsI)})
  
  Ps<-rep(0.05, times=8)
  for (i in 3:10) {
    DT1[,i][DT1[,i]>0]<- (1-Ps[i-2])#(1-exp(-safety3$Ns[i-2]*Ps[i-2]))
    DT1[,i][DT1[,i]==0]<- Ps[i-2]#exp(-safety3$Ns[i-2]*Ps[i-2])
    DQB[,i][DQB[,i]>0]<- (1-Ps[i-2])
    DQB[,i][DQB[,i]==0]<- Ps[i-2]
    DRB[,i][DRB[,i]>0]<- (1-Ps[i-2])
    DRB[,i][DRB[,i]==0]<- Ps[i-2]
  }
  
  DT1$Score<-apply(DT1[3:10], 1, prod)
  DQB$Score<-apply(DQB[3:10], 1, prod)
  DRB$Score<-apply(DRB[3:10], 1, prod)

  amps_list1<-str_split(DT1$amps, pattern = "_")
  #str_split(amps_list1[[1]], pattern = ",")
  for (i in 1:length(amps_list1)) {
    amps_list1[[i]]<-str_split(amps_list1[[i]], pattern = "," )
  }
  amps_listDQB<-str_split(DQB$amps, pattern = "_")
  #str_split(amps_listDQB[[1]], pattern = ",")
  for (i in 1:length(amps_listDQB)) {
    amps_listDQB[[i]]<-str_split(amps_listDQB[[i]], pattern = "," )
  }
  amps_listDRB<-str_split(DRB$amps, pattern = "_")
  #str_split(amps_listDRB[[1]], pattern = ",")
  for (i in 1:length(amps_listDRB)) {
    amps_listDRB[[i]]<-str_split(amps_listDRB[[i]], pattern = "," )
  }
  isdaddy <- function(list1, list2) {
    a<-length(setdiff(list2[[1]], list1[[1]]))==0
    b<-length(setdiff(list2[[2]], list1[[2]]))==0
    c<-length(setdiff(list2[[3]], list1[[3]]))==0
    d<-length(setdiff(list2[[4]], list1[[4]]))==0
    e<-length(setdiff(list2[[5]], list1[[5]]))==0
    f<-length(setdiff(list2[[6]], list1[[6]]))==0
    g<-length(setdiff(list2[[7]], list1[[7]]))==0
    h<-length(setdiff(list2[[8]], list1[[8]]))==0
    daddy<-sum(c(a,b,c,d,e,f,g,h))==8
    daddy
  }
  res1<-matrix(0, ncol=nrow(DT1), nrow=nrow(DT1))
  resDQB<-matrix(0, ncol=nrow(DQB), nrow=nrow(DQB))
  resDRB<-matrix(0, ncol=nrow(DRB), nrow=nrow(DRB))
 
   for (i in 1:length(amps_list1))
    for (j in 1:length(amps_list1)) {
      res1[i,j]<-isdaddy(amps_list1[[i]], amps_list1[[j]])
    } 
  for (i in 1:length(amps_listDQB))
    for (j in 1:length(amps_listDQB)) {
      resDQB[i,j]<-isdaddy(amps_listDQB[[i]], amps_listDQB[[j]])
    } 
  for (i in 1:length(amps_listDRB))
    for (j in 1:length(amps_listDRB)) {
      resDRB[i,j]<-isdaddy(amps_listDRB[[i]], amps_listDRB[[j]])
    } 
  
  gr1<-igraph::simplify(graph_from_adjacency_matrix(res1))
  grDQB<-igraph::simplify(graph_from_adjacency_matrix(resDQB))
  grDRB<-igraph::simplify(graph_from_adjacency_matrix(resDRB))
  
  DT1$degree<-degree(gr1, mode="in")
  DQB$degree<-degree(grDQB, mode="in")
  DRB$degree<-degree(grDRB, mode="in")
  
  list(DT1=DT1, gr1=gr1, DQB=DQB, grDQB=grDQB, DRB=DRB, grDRB=grDRB)
}

get_safety5<-function(DT1, DQB, DRB) {HLA_typing_resultsIclass<-filter(DT1, DT1$degree==0)%>%select(Allele, Score, no_amps, sumreads, meanreads, medianreads)
HLA_typing_resultsDQB<-filter(DQB, DQB$degree==0)%>%select(Allele, Score, no_amps, sumreads, meanreads, medianreads)
HLA_typing_resultsDRB<-filter(DRB, DRB$degree==0)%>%select(Allele, Score, no_amps, sumreads, meanreads, medianreads)
HLA_all<-list(HLA_typing_resultsIclass, HLA_typing_resultsDQB, HLA_typing_resultsDRB)
HLA_typing_results<-do.call(rbind, HLA_all)

HLA_typing_results[grepl(HLA_typing_results$Allele, pattern = "A*", fixed = T), ]$Score<-prop.table(HLA_typing_results[grepl(HLA_typing_results$Allele, pattern = "A*", fixed = T), ]$Score)
HLA_typing_results[grepl(HLA_typing_results$Allele, pattern = "B*",fixed = T), ]$Score<-prop.table(HLA_typing_results[grepl(HLA_typing_results$Allele, pattern = "B*", fixed = T), ]$Score)
HLA_typing_results[grepl(HLA_typing_results$Allele, pattern = "C*", fixed = T), ]$Score<-prop.table(HLA_typing_results[grepl(HLA_typing_results$Allele, pattern = "C*", fixed = T), ]$Score)
HLA_typing_results[grepl(HLA_typing_results$Allele, pattern = "DQB", fixed = T), ]$Score<-prop.table(HLA_typing_results[grepl(HLA_typing_results$Allele, pattern = "DQB", fixed = T), ]$Score)
HLA_typing_results[grepl(HLA_typing_results$Allele, pattern = "DRB1", fixed = T), ]$Score<-prop.table(HLA_typing_results[grepl(HLA_typing_results$Allele, pattern = "DRB1", fixed = T), ]$Score)
HLA_typing_results[grepl(HLA_typing_results$Allele, pattern = "DRB[2-9]", fixed = F), ]$Score<-prop.table(HLA_typing_results[grepl(HLA_typing_results$Allele, pattern = "DRB[2-9]", fixed = F), ]$Score)


HLA_typing_results<-select(HLA_typing_results, Allele, Score, no_amps, sumreads, meanreads, medianreads)
}

#TypingResults2<-lapply(amplicones_miseq_new_lists_all, get_HLA_result)

# Mega_table_hlares<-do.call(rbind, lapply(TypingResults2, function(x){x$safety5}))
# Mega_table_hlares<-rownames_to_column(df = Mega_table_hlares)
# Mega_table_hlares$donor<-sapply(str_split(Mega_table_hlares$rowname, fixed(".")), function(x) {x[[1]]})
# Mega_table_hlares<-select(Mega_table_hlares, donor, Allele, Score, no_amps, sumreads, meanreads, medianreads)

signs3<-function(string) {
  paste(str_split(string, pattern = fixed(":"), simplify = T )[1:3], sep=":", collapse = ":") 
  
}

gettidyHLA<-function (longstr) {
  paste(unique(sapply(str_split(longstr, pattern=fixed(" "), simplify = T), signs3)), collapse = " ")
}

# Mega_table_hlares$tidyAllele<-sapply(Mega_table_hlares$Allele, gettidyHLA)
# Mega_table_hlares$tidyAllele<-gsub(Mega_table_hlares$tidyAllele, pattern = ":NA", replacement = "", fixed = T)
# Mega_table_hlares<-Mega_table_hlares[Mega_table_hlares$Score>1e-5, ]
# Mega_table_hlares<-select(Mega_table_hlares, donor, tidyAllele, Allele, Score, no_amps, sumreads, meanreads, medianreads)
# Luci<-get_HLA_result(amplicones_miseq_new_lists_all$Luci)
# SL<-get_HLA_result(amplicones_miseq_new_lists_all$SL)


