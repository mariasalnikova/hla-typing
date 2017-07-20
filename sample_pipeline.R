#sample pipeline
#single file
#tst<-HLA_amplicones_full("230617_MiSeq/9_S4_L001_R1_001.fastq.gz","230617_MiSeq/9_S4_L001_R2_001.fastq.gz",threshold=100,read_length=250)
#tstgz<-gzreader("230617_MiSeq/Mstr_S8_L001_R1_001.fastq.gz","230617_MiSeq/Mstr_S8_L001_R2_001.fastq.gz")
#table(nchar(tstgz[,1]))
File_list_pipeline_amplicones<-function(filelist,read_length=250,threshold=100){
  resHLA<-lapply(filelist[,1],"[",1)
  names(resHLA)<-filelist[,1]
  for (i in 1:nrow(filelist)){
    print(names(resHLA)[i])
    resHLA[[i]]<-HLA_amplicones_full(filelist[i,2],filelist[i,3],read_length = read_length,threshold=threshold)
  }
  resHLA
}
# easy13<-HLA_amplicones_full("MiSeq100717/13_S7_L001_R1_001.fastq.gz","MiSeq100717/13_S7_L001_R2_001.fastq.gz",threshold=100,read_length=250)
# easy13_hla<-get_HLA_result(easy13)
miseq1007_safety2_compl<-File_list_pipeline_amplicones(miseq100717)
save(miseq1007_safety2_compl, file="miseq1007_safety2_compl.rda")
TypingMiseq100717<-lapply(Nmiseq1007_safety2_compl, get_HLA_result)
save(TypingMiseq100717, file="TypingMiseq100717.rda")
Mega_table_100717<-do.call(rbind, lapply(TypingMiseq100717, function(x){x$safety5}))
Mega_table_100717<-rownames_to_column(df = Mega_table_100717)
Mega_table_100717$donor<-sapply(str_split(Mega_table_100717$rowname, fixed(".")), function(x) {x[[1]]})
Mega_table_100717<-select(Mega_table_100717, donor, Allele, Score, no_amps, sumreads, meanreads, medianreads)
Mega_table_100717$tidyAllele<-sapply(Mega_table_100717$Allele, gettidyHLA)
Mega_table_100717$tidyAllele<-gsub(Mega_table_100717$tidyAllele, pattern = ":NA", replacement = "", fixed = T)
#Mega_table_100717<-Mega_table_100717[Mega_table_100717$Score>1e-5, ]
Mega_table_100717<-select(Mega_table_100717, donor, tidyAllele, Allele, Score, no_amps, sumreads, meanreads, medianreads)
save(Mega_table_100717, file="Mega_table_100717.rda")

# View(tst$Iclass)
# View(tst$IIclassDQB)
# View(tst$IIclassOthers)
# View(tst$Iclass_alt)
# View(tst$IIclassDQB_alt)
# View(tst$IIclassOthers_alt)
# View(tst$safety1$Iamp1)
# View(tst$safety1$Iamp2)
# View(tst$safety1$Iamp1alt)
# View(tst$safety1$Iamp2alt)
# View(tst$safety1$Iamp1_inv)
# View(tst$safety1$Iamp2_inv)
# View(tst$safety1$Iamp1alt_inv)
# View(tst$safety1$Iamp2alt_inv)
# 
# tst2<-Alleles_bayes_container(tst$safety2)
# View(do.call(rbind,lapply(tst2,function(x)x[rel_fibi>0.01,,][order(-rel_fibi)])))
# View(tst2$A[rel_fibi>0.01,,][order(-rel_fibi)])
# View(tst2$C[rel_fibi>0.01,,][order(-rel_fibi)])
# View(tst2$B[rel_fibi>0.01,,][order(-rel_fibi)])
# View(tst2$DQB[rel_fibi>0.01,,][order(-rel_fibi)])
# View(tst2$DRB[rel_fibi>0.01,,][order(-rel_fibi)])
# 
# View(tst$safety2$Iamp1)
# DT<-mutate(DT, rank=rank(DT$degree), group_by(DT$cluster_id))
# View(DT)
