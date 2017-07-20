load("HLA_base.rda")
library(igraph)
library(stringdist)
library(Biostrings)
library(data.table)
library(dplyr)

HLA_set<-DNAStringSet(HLA_base$Sequence)

#Function for reverse complement
revcomp<-function (.seq) 
{
  rc.table <- c(A = "T", T = "A", C = "G", G = "C",N="N")
  sapply(strsplit(.seq, "", T, F, T), function(l) paste0(rc.table[l][length(l):1], 
                                                         collapse = ""), USE.NAMES = F)
}

#Get table with read1 and read2
gzreader<-function(read1,read2){
  #supernaive reader:
  #gzread1 <- gzfile(read1, open = "r")
  #gzread2 <- gzfile(read2, open = "r")
  Fastq1<-readLines(read1);
  Fastq2<-readLines(read2);
  res<-data.frame(read1=Fastq1[(1:length(Fastq1))%%4==2],read2=Fastq2[(1:length(Fastq2))%%4==2],stringsAsFactors = F)
  rm(Fastq1)
  rm(Fastq2)
  #close(gzread1)
  #close(gzread2)
  res
}

get_groups_complicated<-function(readed,threshold=100,read_length=250,cons_length=150,primer_length=20){
  readed$read1<-substr(readed$read1,primer_length,read_length)#trim them
  readed$read2<-substr(readed$read2,primer_length,read_length)#trim them
  readed$cons1<-substr(readed$read1,primer_length,cons_length)
  readed<-readed[duplicated(readed$cons1),]
  readed<-split(readed,f = readed$cons1) 
  #merge them
  if (length(readed)>0)
    readed<-readed[(rank(-sapply(readed,nrow),ties.method = "first")<=threshold)&(sapply(readed,nrow)>2)]
  readed
}
merge_cons_complicated<-function(readed){
  consensus1<-character(length(readed))
  consensus2<-character(length(readed))
  if(length(readed)>0)
    for (i in 1:length(readed))
    {#print(i)
      consensus1[i]<-docons_corr(readed[[i]]$read1)
      consensus2[i]<-docons_corr(readed[[i]]$read2)
      
    }
  #consensus
  data.frame(read1=consensus1,read2cons=consensus2,readnumber=sapply(readed,nrow),stringsAsFactors = F)
}

#merge them, sort them, filter unique in crosstab fashion. trim them. 
#Trim reads, get nonunique read1, get champions for read 1
get_groups<-function(readed,threshold=100,read_length=250,primer_length=20){
  readed$read1<-substr(readed$read1,primer_length,read_length)#trim them
  readed$read2<-substr(readed$read2,primer_length,read_length)#trim them
  readed<-readed[duplicated(readed$read1),]
  readed<-split(readed,f = readed$read1) 
  #merge them
  if (length(readed)>0)
    readed<-readed[(rank(-sapply(readed,nrow),ties.method = "first")<=threshold)&(sapply(readed,nrow)>2)]
  readed
}


alphabet<-c(A=1,G=2,C=3,T=4,N=5)
#Function for consensus
docons<-function(vec,thres=0.7){
  str_mat<-do.call(rbind,(strsplit(vec,"")))
  paste0(sapply(apply(str_mat,MARGIN = 2,table, simplify=F),
                function(x)
                  {if(max(prop.table(x))>thres)
                    {names(x)[x==max(x)][1]}
                  else{"N"}}),collapse = "")
}

docons_corr<-function(vec,thres=0.7){
  str_mat<-do.call(rbind,(strsplit(vec,"")))
  paste0(sapply(lapply(1:ncol(str_mat), function(y) {table(str_mat[, y])}),
                function(x)
                {if(max(prop.table(x))>thres)
                {names(x)[x==max(x)][1]}
                  else{"N"}}),collapse = "")
}


#Get consensus for read2
merge_cons<-function(readed){
  consensus<-character(length(readed))
  if(length(readed)>0)
    for (i in 1:length(readed))
    {#print(i)
      consensus[i]<-docons_corr(readed[[i]]$read2)
    }
  consensus
  data.frame(read1=names(readed),read2cons=consensus,readnumber=sapply(readed,nrow),stringsAsFactors = F)
}

N_content<-function(str){
  table(c("N",unlist(strsplit(str,""))))["N"]-1
}

#Search for overlap
LCS_comp2<-function(a,b){#one side overlap a b
  a<-alphabet[unlist(strsplit(a,""))]
  b<-alphabet[unlist(strsplit(b,""))]
  res<-a
  for (i in 1:length(a))
  {
    res[i]<-suppressWarnings(max(rle(a[i:length(a)]-b[1:(length(a)-i+1)])$lengths[rle(a[i:length(a)]-b[1:(length(a)-i+1)])$values==0],na.rm=T))
  }
  c(max(res),which(res==max(res)))
}


#Merge overlap reads 1 and 2
get_overlap_merged<-function(merged){
  temp<-list();
  if(nrow(merged)==0)return(merged)
  
  for (i in 1:nrow(merged))temp[[i]]<-LCS_comp2(merged$read1[i],revcomp(merged$read2cons[i]))
  temp<-do.call(rbind,temp)
  merged$Overlap_max_aligned<-temp[,1]
  merged$Overlap_start<-temp[,2]
  merged$assembled<-paste0(merged$read1,substr(revcomp(merged$read2cons),nchar(merged$read1)-merged$Overlap_start+2,nchar(merged$read2cons)))
  merged
}

get_overlap_merged_fix<-function(merged,shift=1, rlength=231){
  temp<-list();
  if(nrow(merged)==0)return(merged)
  #for (i in 1:nrow(merged))temp[[i]]<-LCS_comp2(merged$read1[i],revcomp(merged$read2cons[i]))
  #temp<-do.call(rbind,temp)
  merged$Overlap_max_aligned<-shift
  merged$Overlap_start<-rlength-shift+1
  merged$assembled<-paste0(merged$read1,substr(revcomp(merged$read2cons),nchar(merged$read1)-merged$Overlap_start+2,nchar(merged$read2cons)))
  merged
}

#Merge nonoverlap read 1 and 2 
get_nonoverlap_merged<-function(merged,insert=1){
  if(nrow(merged)==0)return(merged)
  merged$Overlap_max_aligned<-insert
  merged$Overlap_start<-nchar(merged$read1[1])+insert
  #merged$assembled<-paste0(merged$read1,substr(revcomp(merged$read2cons),nchar(merged$read1)-merged$Overlap_start+2,nchar(merged$read2cons)))
  N_vec<-paste0(rep("N", times=insert), collapse="")
  merged$assembled<-paste0(merged$read1,N_vec,revcomp(merged$read2cons))
  merged
}

#Get assembled exact from base
get_exact<-function(merged){
  if (nrow(merged)==0)return(merged)
  resM<-character(nrow(merged))
  for (i in 1:length(resM))
    resM[i]<-paste0(HLA_base[grepl(gsub("N","[ATGC]",merged$assembled[i],fixed = T),x = HLA_base$Sequence,fixed=F),]$Allele,collapse = " ")
  merged$Exact<-resM
  merged
}

#Get reversed assembled exact from base
get_exact_rev<-function(merged){
  if (nrow(merged)==0)return(merged)
  
  resM<-character(nrow(merged))
  for (i in 1:length(resM))
    resM[i]<-paste0(HLA_base[grepl(gsub("N","[ATGC]",revcomp(merged$assembled[i]),fixed = T),x = HLA_base$Sequence,fixed=F),]$Allele,collapse = " ")
  merged$Exact_rev<-resM
  merged
}

#Make seq graph
make.sequence.graph <- function (.data, .name = '',max_errs=1) {
  G <- graph.empty(n = length(.data), directed=F)
  tmp<-stringdistmatrix(.data,.data,method="hamming")
  G <- add.edges(G, t(which(tmp==max_errs,arr.ind=T)))
  G <- igraph::simplify(G)
  G <- set.vertex.attribute(G, 'label', V(G), .data)
  #print(G)
  G
}


#Search for NGS and PCR mistakes
fathers_and_children<-function(readed){#graph and cluster analysis
  if (nrow(readed)==0)return(readed)
  Gr<-make.sequence.graph(readed$read1)
  readed$neighbours_read1<-degree(Gr)
  cl<-clusters(Gr)
  readed$clusters<-cl$membership
  readed$freq<-prop.table(readed$readnumber)
  parents<-numeric(nrow(readed))
  for (clust in unique(cl$membership))
    parents[readed$clusters==clust]<-rank(-readed[readed$clusters==clust,]$neighbours_read1,ties.method = "min")
  readed$parents<-parents
  readed
}

#Ungay the sequences
straight_inverse<-function(readed){
  if (nrow(readed)==0)return(readed)
  readed$Exact<-readed$Exact_rev
  readed$assembled<-revcomp(readed$assembled)
  readed
}

#Find the butt of amps
make_overlaps_one_side_grep<-function(seq1,seq2,length=50){
  Overlapstart<-regexpr(gsub("N","[AGCT]",substr(seq2,1,length) ,fixed=T),seq1)
  Overlaplength<-nchar(seq1)-Overlapstart+1
  over1<-c(Overlaplength,Overlapstart)
  if(over1[2]!=-1)
  {
    data.frame(seq=paste0(substr(seq1,1,over1[2]-1),seq2),
               Overlap_length=over1[1],Overlap_start=over1[2], 
               Overlap_mismatches=stringdist(substr(seq1,over1[2],nchar(seq1)),substr(seq2,1,nchar(seq1)-over1[2]+1),method = "hamming")-max(N_content(substr(seq1,over1[2],nchar(seq1))),N_content(substr(seq2,1,nchar(seq1)-over1[2]+1))),stringsAsFactors = F)
  }
}

#Merge amp1 and amp2
intersect_amplicones3<-function(first_amp_list,second_amp_list){
  #res<-list()
  res<-data.frame()
  if ((nrow(first_amp_list)==0)|(nrow(second_amp_list)==0)) return (res)
  for (i in 1:nrow(first_amp_list)){#res[[i]]<-list();#print(i)
    for (j in 1:nrow(second_amp_list))
      # if (i<j)
      #if (sum(duplicated(c(unlist(strsplit(as.character(exact_list[i,]$Exact)," ")),unlist(strsplit(as.character(exact_list[j,]$Exact)," ")))))!=0)
    {
      over<-make_overlaps_one_side_grep(first_amp_list[i,]$assembled,second_amp_list[j,]$assembled)
      if(length(over)!=0){
        over$first_reads=first_amp_list[i,]$readnumber;
        over$second_reads=second_amp_list[j,]$readnumber;
        over$first_reads_freq=first_amp_list[i,]$freq;
        over$second_reads_freq=second_amp_list[j,]$freq;
        over$mean_geom_freq=sqrt(first_amp_list[i,]$freq*second_amp_list[j,]$freq)
        over$geom_mean_reads=sqrt(first_amp_list[i,]$readnumber*second_amp_list[j,]$readnumber)
        over$geom_min_coverage=min(first_amp_list[i,]$readnumber,second_amp_list[j,]$readnumber)
        over$first_assembled=first_amp_list[i,]$assembled
        over$second_assembled=second_amp_list[j,]$assembled
        over$Ns=sapply(over$seq,N_content)
        over$Length=nchar(over$seq)
        #over$Exact_match<-paste0(HLA_base[grepl(gsub("N","[AGCT]",over$seq,fixed=T),x = HLA_base$Sequence),]$Allele,collapse = " ")
        over$Exact_match<-paste0(HLA_base[vcountPattern(pattern = over$seq,subject = HLA_set,algorithm = "auto",fixed = F,max.mismatch = 0,min.mismatch = 0,with.indels = F)!=0,]$Allele,collapse = " ")
        over$One_mismatch<-paste0(HLA_base[vcountPattern(pattern = over$seq,subject = HLA_set,algorithm = "auto",fixed = F,max.mismatch = 1,min.mismatch = 1,with.indels = F)!=0,]$Allele,collapse = " ")
        over$Two_mismatch<-paste0(HLA_base[vcountPattern(pattern = over$seq,subject = HLA_set,algorithm = "auto",fixed = F,max.mismatch = 2,min.mismatch = 2,with.indels = F)!=0,]$Allele,collapse = " ")
        res<-rbind(res,over)
        #res[[i]][[j]]<-over
      }
    }
  }
  res
  #do.call(rbind,do.call(rbind,res))
}

#Aggregate_data
Alleles_bayes_container<-function(safety2){
  tstluk2_all<-Alleles_bayes(do.call(rbind,safety2[grepl("amp1",names(safety2))]),do.call(rbind,safety2[grepl("amp2",names(safety2))]))
  bayesian_alleles_report(tstluk2_all)
}

mismatch_vector<-function(seq,err_prob=0.001,max_mismatch=3){ #returns a vector???
  mism<-rep(0.0000000001,nrow(HLA_base)) #set very low prob 
  mism[vcountPattern(pattern = seq,subject = HLA_set,algorithm = "auto",fixed = F,max.mismatch = 0,min.mismatch = 0,with.indels = F)!=0]<-1 #if something is exactly matched it gets weight 1
  for (i in 1:max_mismatch)
    mism[vcountPattern(pattern = seq,subject = HLA_set,algorithm = "auto",fixed = F,max.mismatch = i,min.mismatch = i,with.indels = F)!=0]<-err_prob**i #if something is not exactly matched it gets weight err_prob**i
  #mism<-mism/sum(mism,na.rm = T)
  prop.table(mism)
}

mismatch_matrix<-function(amplist,err_prob=0.001,max_mismatch=3){
  do.call(cbind,lapply(amplist$assembled,mismatch_vector,err_prob=0.001,max_mismatch=3))
}

Alleles_bayes<-function(amp1list, amp2list){
  #get mismat.
  amp1list<-amp1list[!grepl("NULL",amp1list$assembled),]
  amp2list<-amp2list[!grepl("NULL",amp2list$assembled),]
  
  #fi get fi
  amp1mism<-mismatch_matrix(amp1list)
  amp2mism<-mismatch_matrix(amp2list)
  amp1maps<-(apply(amp1mism,MARGIN = 1,function(x){paste0(which(x>0.001),collapse=" ")}))
  amp2maps<-(apply(amp2mism,MARGIN = 1,function(x){paste0(which(x>0.001),collapse=" ")}))
  fi<-amp1list$freq%*%t(amp1mism) #why is that? it could be replaced with some prob, that given seq is error from bigger seq (matrix n_seq by n_seq), ratios? and then we simply multiply each by rowsums?== integrate over n2
  bi<-amp2list$freq%*%t(amp2mism) #PROB THAT IT IS TRUE SEQUENCE. 
  
  #list(as.vector(fi),as.vector(bi),as.vector(fi*bi))
  as.data.table(data.frame(fi=as.vector(fi),bi=as.vector(bi),fibi=as.vector(fi)*as.vector(bi),allele=HLA_base$Allele,class=HLA_base$HLA_class,a1=amp1maps,a2=amp2maps,a1_a2=paste(amp1maps,amp2maps,sep="_")))
  #get mi.
  #tstluk[class=="C",,][prop.table(fibi)>0.01,,][,list(amps=paste0(allele,collapse=" "),fibi=sum(fibi)),a1_a2]
}
bayesian_alleles_report<-function(allelelist){
  a<-allelelist[class=="A",,][,,][,list(amps=paste0(allele,collapse=" "),fibi=sum(fibi)),a1_a2][,rel_fibi:=prop.table(fibi),][order(-rel_fibi)]
  b<-allelelist[class=="B",,][,,][,list(amps=paste0(allele,collapse=" "),fibi=sum(fibi)),a1_a2][,rel_fibi:=prop.table(fibi),][order(-rel_fibi)]
  c<-allelelist[class=="C",,][,,][,list(amps=paste0(allele,collapse=" "),fibi=sum(fibi)),a1_a2][,rel_fibi:=prop.table(fibi),][order(-rel_fibi)]
  DRB<-allelelist[grepl("DRB",class),,][,,][,list(amps=paste0(allele,collapse=" "),fibi=sum(fibi)),a1_a2][,rel_fibi:=prop.table(fibi),][order(-rel_fibi)]
  DQB<-allelelist[class=="DQB1",,][,,][,list(amps=paste0(allele,collapse=" "),fibi=sum(fibi)),a1_a2][,rel_fibi:=prop.table(fibi),][order(-rel_fibi)]
  list(A=a,B=b,C=c,DRB=DRB,DQB=DQB)
}

HLA_amplicones_full<-function(read1,read2,threshold=100,read_length=250){
  print("File read start")
  print(format(Sys.time(), "%a %b %d %X %Y"))
  readed<-gzreader(read1,read2)
  readed<-readed[nchar(readed$read1)>200&nchar(readed$read2)>200, ]
  nrow(readed)
  readed<-list(Iamp1=readed[grepl("CCCTGACC[GC]AGACCTG",substr(readed$read1,1,20)),],
               Iamp2=readed[grepl("CGACGGCAA[AG]GATTAC",substr(readed$read1,1,20)),],
               IIamp1_DQB=readed[grepl("AG[GT]CTTTGCGGATCCC",substr(readed$read1,1,20)),],
               IIamp1_Others=readed[grepl("CTGAGCTCCC[GC]ACTGG",substr(readed$read1,1,20)),],
               IIamp2=readed[grepl("GGAACAGCCAGAAGGA",substr(readed$read1,1,20)),],
               Iamp1alt=readed[grepl("TC[CT]CACTCCATGAGGTATTTC|TCCCACTCCATGAAGTATTTC",substr(readed$read1,1,22)),],
               Iamp2alt=readed[grepl("GGCAA[AG]GATTACATCGCC|GGCAAGGATTACATCGCT",substr(readed$read1,1,22)),],
               IIamp1_DQBalt=readed[grepl("TGAGGGCAGAGAC[CT]CTCC",substr(readed$read1,1,22)),],
               IIamp1_DRBalt=readed[grepl("TGACAGTGACACTGATGG|TGACAGTGACATTGACGG",substr(readed$read1,1,22)),],
               IIamp2_DRBalt=readed[grepl("GAGAGCTTCAC[AG]GTGCAG",substr(readed$read1,1,22)),],#&grepl("TG[CT]TCTGGGCAGATTCAG",substr(readed$read2,1,22))
               IIamp2_DQBalt=readed[grepl("ACCATCTCCCCATCCAG",substr(readed$read1,1,22)),],#&grepl("TG[CT]TCTGGGCAGATTCAG",substr(readed$read2,1,22))
               Iamp1_inv=readed[grepl("GGGCCGCCTCC[AC]ACTTG|GGGCCGTCTCCCACTTG|GGACCGCCTCCCACTTG",substr(readed$read1,1,20)),],
               Iamp2_inv=readed[grepl("[CT]GGTGG[AG]CTGGGAAGA",substr(readed$read1,1,20)),],
               IIamp1_DQB_inv=readed[grepl("[CT]CAGCAGGTTGTGGTG|CCAG[GC]AGGTT[AG]TGGTG",substr(readed$read1,1,20)),],#&grepl("AG[GT]CTTTGCGGATCCC",substr(readed$read2,1,20)),],
               IIamp1_Others_inv=readed[grepl("[CT]CAGCAGGTTGTGGTG|CCAG[GC]AGGTT[AG]TGGTG",substr(readed$read1,1,20)),],#&grepl("CTGAGCTCCC[GC]ACTGG",substr(readed$read2,1,20)),],
               IIamp2_inv=readed[grepl("CCAC[GT]TGGCAGGTGTA|CCACTTGGCAAGTGTA",substr(readed$read1,1,20)),],
               Iamp1alt_inv=readed[grepl("GAGC[GC]ACTCCACGCAC|GAGCCCGTCCACGCAC",substr(readed$read1,1,22)),],
               Iamp2alt_inv=readed[grepl("TCAGGGTGAGGGGCT|TCAGGGTGCAGGGCT",substr(readed$read1,1,22)),],
               IIamp1_DQBalt_inv=readed[grepl("GTCCAGTCACC[AG]TTCCTA",substr(readed$read1,1,22)),],
               IIamp1_DRBalt_inv=readed[grepl("CAG[CT]CTTCTCTTCCTGGC",substr(readed$read1,1,22)),],
               IIamp2_DRBalt_inv=readed[grepl("TGCTCTGTGCAGATTCAG",substr(readed$read1,1,22)),],#&grepl("TG[CT]TCTGGGCAGATTCAG",substr(readed$read1,1,22)),],
               IIamp2_DQBalt_inv=readed[grepl("TG[CT]TCTGGGCAGATTCAG",substr(readed$read1,1,22)),]#&grepl("TG[CT]TCTGGGCAGATTCAG",substr(readed$read1,1,22)),]
  )
  ampreads<-sapply(readed,nrow)
  print(ampreads)
  allreads<-sum(sapply(readed,nrow))
  print(allreads)
  print("File read end. Looking for champions!")
  print(format(Sys.time(), "%a %b %d %X %Y"))
  readed<-lapply(readed,get_groups,threshold = threshold,read_length = read_length, primer_length=20)
    readed<-lapply(readed, merge_cons)
  print("We are the champions!")
  print(format(Sys.time(), "%a %b %d %X %Y"))
  nchamp<-sapply(readed,nrow)
  champ_reads<-sapply(readed, function(x) {sum(x$readnumber)})
  print(nchamp)
  #readed<-lapply(readed,get_overlap_merged)
  
  readed$Iamp1<-get_overlap_merged_fix(readed$Iamp1, shift=29)
  readed$Iamp1_inv<-get_overlap_merged_fix(readed$Iamp1_inv, shift=29)
  readed$Iamp2<-get_overlap_merged_fix(readed$Iamp2, shift=9)
  readed$Iamp2_inv<-get_overlap_merged_fix(readed$Iamp2_inv, shift=9)
  readed$IIamp1_DQB<-get_overlap_merged_fix(readed$IIamp1_DQB, shift=71)
  readed$IIamp1_DQB_inv<-get_overlap_merged_fix(readed$IIamp1_DQB_inv, shift = 71)
  readed$IIamp1_Others<-get_overlap_merged_fix(readed$IIamp1_Others, shift = 130)
  readed$IIamp1_Others_inv<-get_overlap_merged_fix(readed$IIamp1_Others_inv, shift=130)
  readed$IIamp2<-get_overlap_merged_fix(readed$IIamp2, shift = 155)
  readed$IIamp2_inv<-get_overlap_merged_fix(readed$IIamp2_inv, shift = 155)
  readed$Iamp2alt<-rbind(get_overlap_merged_fix(readed$Iamp2alt, shift=9), get_overlap_merged_fix(readed$Iamp2alt, shift=40))
  readed$Iamp2alt_inv<-rbind(get_overlap_merged_fix(readed$Iamp2alt_inv, shift = 9), get_overlap_merged_fix(readed$Iamp2alt_inv, shift = 40))
  readed$IIamp1_DQBalt<-get_overlap_merged_fix(readed$IIamp1_DQBalt, shift = 32)
  readed$IIamp1_DQBalt_inv<-get_overlap_merged_fix(readed$IIamp1_DQBalt_inv, shift = 32)
  readed$IIamp1_DRBalt<-get_overlap_merged_fix(readed$IIamp1_DRBalt, shift=35)
  readed$IIamp1_DRBalt_inv<-get_overlap_merged_fix(readed$IIamp1_DRBalt_inv, shift = 35)
  readed$IIamp2_DQBalt<-get_overlap_merged_fix(readed$IIamp2_DQBalt, shift = 205)
  readed$IIamp2_DQBalt_inv<-get_overlap_merged_fix(readed$IIamp2_DQBalt_inv, shift = 205)
  readed$IIamp2_DRBalt<-get_overlap_merged_fix(readed$IIamp2_DRBalt, shift=166)
  readed$IIamp2_DRBalt_inv<-get_overlap_merged_fix(readed$IIamp2_DRBalt_inv, shift =166)
  
  readed$Iamp1alt<-get_nonoverlap_merged(readed$Iamp1alt)
  readed$Iamp1alt_inv<-get_nonoverlap_merged(readed$Iamp1alt_inv)
  print("Reads_assembled")
  print(format(Sys.time(), "%a %b %d %X %Y"))
  readed<-lapply(readed,get_exact)
  readed<-lapply(readed,get_exact_rev) 
  print("Exact matches read1 found")
  print(format(Sys.time(), "%a %b %d %X %Y"))
  readed<-lapply(readed,fathers_and_children)
  readed[12:22]<-lapply(readed[12:22],straight_inverse)
  # readed<-amplist
  readed_intersect<-lapply(readed,function(x){x[(x$parents==1|(x$parents==2&x$freq>0.05)),]})
  print("Graph made. Ready for papas intersection")
  print(format(Sys.time(), "%a %b %d %X %Y"))
  npapas<-sapply(readed_intersect,nrow)
  intersect_reads<-sapply(readed_intersect, function(x) {sum(x$readnumber)})
  print(npapas)
  non_used<-sapply(readed_intersect, function(x) {sum(x[x$Exact=="", ]$readnumber)})
  statistics<-list(n_reads=allreads, amps_reads=ampreads, n_champs=nchamp, champ_reads=champ_reads, n_intersected=npapas, reads_intersected=intersect_reads, non_used=non_used)
  res<-list(safety1=readed,
            safety2=readed_intersect, statistics=statistics
            #,Iclass=intersect_amplicones3(rbind(readed_intersect$Iamp1,readed_intersect$Iamp1_inv),rbind(readed_intersect$Iamp2,readed_intersect$Iamp2_inv)),
            #IIclassDQB=intersect_amplicones3(rbind(readed_intersect$IIamp1_DQB,readed_intersect$IIamp1_DQB_inv),rbind(readed_intersect$IIamp2,readed_intersect$IIamp2_inv)),
            #IIclassOthers=intersect_amplicones3(rbind(readed_intersect$IIamp1_Others,readed_intersect$IIamp1_Others_inv),rbind(readed_intersect$IIamp2,readed_intersect$IIamp2_inv)),
            #Iclass_inv=intersect_amplicones3(readed_intersect$Iamp1_inv,readed_intersect$Iamp2_inv),
            #IIclassDQB_inv=intersect_amplicones3(readed_intersect$IIamp1_DQB_inv,readed_intersect$IIamp2_inv),
            #IIclassOthers_inv=intersect_amplicones3(readed_intersect$IIamp1_Others_inv,readed_intersect$IIamp2_inv),
            
            #Iclass_alt=intersect_amplicones3(rbind(readed_intersect$Iamp1alt,readed_intersect$Iamp1alt_inv),rbind(readed_intersect$Iamp2alt,readed_intersect$Iamp2alt_inv)),
            #IIclassDQB_alt=intersect_amplicones3(rbind(readed_intersect$IIamp1_DQBalt,readed_intersect$IIamp1_DQBalt_inv),rbind(readed_intersect$IIamp2_DQBalt,readed_intersect$IIamp2_DQBalt_inv)),
            #IIclassOthers_alt=intersect_amplicones3(rbind(readed_intersect$IIamp1_DRBalt,readed_intersect$IIamp1_DRBalt_inv),rbind(readed_intersect$IIamp2_DRBalt,readed_intersect$IIamp2_DRBalt_inv))
            # Iclass_inv_alt=intersect_amplicones3(readed_intersect$Iamp1alt_inv,readed_intersect$Iamp2alt_inv),
            # IIclassDQB_inv_alt=intersect_amplicones3(readed_intersect$IIamp1_DQBalt_inv,readed_intersect$IIamp2_DQBalt_inv),
            #IIclassOthers_inv_alt=intersect_amplicones3(readed_intersect$IIamp1_DRBalt_inv,readed_intersect$IIamp2_DRBalt_inv)
  )
  print("Amplicones intersected. Report generated.")
  print(format(Sys.time(), "%a %b %d %X %Y"))
  res
}  



