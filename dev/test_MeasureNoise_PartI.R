#!/usr/bin/env Rscript

###########################
## Test script to test the effect of noise
## Part I (before dtw)
###########################
#Created 2021-01-08

setwd("~/Dropbox/Arabidopsis-eval/R_simul-eval/example_data")
source("~/Dropbox/Arabidopsis-eval/R_simul-eval/source/sim_phyllo_sources.R")
source('~/Dropbox/Arabidopsis-eval/R_simul-eval/source/eval_dtw_sources.R')

#####
#1. simul
####
#Initial (biological) sequence
N1=25
alpha=137.5
a_sd=18.5
i_Gsd=0.8
i_noise_pct=75
seq=make_refseq(N1, alpha, a_sd, i_Gsd, i_noise_pct)
listN1=make_align_list(N1)
#segmentation errors
GAIN=c(10,15)
LOSS=c(1,2,3,8,18,21:26)

#Measures or Noise
measure=FALSE
noise=TRUE
meanA=mean(seq$angles)
meanI=mean(seq$internodes)

#Noise levels
Noise_levels=seq(0,1,0.05)

#Generate tables of refseq and testseq
#Initialize output tables:
refseq.table=data.frame(plantID=NULL,
                  angles=NULL,
                  internodes=NULL)
testseq.table=data.frame(plantID=NULL,
                   angles=NULL,
                   internodes=NULL)
align.intervals=data.frame(plantID=NULL,
                           reference=NULL,
                           modified=NULL,
                           dtw=NULL)

for (i in 1:length(Noise_levels)){
  if (measure){
    #With a measure
    manual_anoise_sd=6 #(in degree, noise with gaussian distrib. of zero mean)
    manual_inoise_sd=0.5 #(in mm, noise with gaussian distrib. of zero mean)
    aut_anoise_sd=60 #(in degree, noise with gaussian distrib. of zero mean) -> x10 manual
    aut_inoise_sd=10 # -> x 10 manual
    
    seq.ref=make_measure(seq, manual_anoise_sd, manual_inoise_sd)
    seq.aut=make_measure(seq, aut_anoise_sd, aut_inoise_sd)
    seq.test=segmentation_errors(seq.aut,listN1,
                                 organ_gain=GAIN,
                                 organ_loss=LOSS) 
  } else if (noise) {
    anoise_sd=Noise_levels[i]*meanA #(in degree, noise with gaussian distrib. of zero mean)
    inoise_sd=Noise_levels[i]*meanI #(in mm, noise with gaussian distrib. of zero mean)
    seq.noise=make_measure(seq, anoise_sd, inoise_sd)
    seq.test=segmentation_errors(seq.noise,listN1,
                                 organ_gain=GAIN,
                                 organ_loss=LOSS)
  } else {
    #Without noise or measure  
    seq.ref=seq
    seq.test=segmentation_errors(seq.ref,listN1,
                                 organ_gain=GAIN,
                                 organ_loss=LOSS)
  }
  
  refseq.table=rbind.data.frame(refseq.table,
                                cbind.data.frame(plantID=rep(paste0("Plant_",Noise_levels[i]), nrow(seq)), 
                                           seq[,c(2,3)])) #drop the "intervals" column
  testseq.table=rbind.data.frame(testseq.table,
                                 cbind.data.frame(plantID=rep(paste0("Plant_",Noise_levels[i]), nrow(seq.test$values)), 
                                                  seq.test$values[,c(2,3)])) #drop the "intervals" column
  align.intervals=rbind.data.frame(align.intervals,
                                   cbind.data.frame(plantID=rep(paste0("Plant_",Noise_levels[i]), nrow(seq.test$I)), 
                                                    seq.test$I))
}

#####
#2. write out tables
####
colnames(refseq.table)=c("Plant ID", "angles", "Internodes")
write.csv(refseq.table, file="reference_sequences.csv", row.names = FALSE)
colnames(testseq.table)=c("Plant ID", "angles", "Internodes")
write.csv(testseq.table, file="test_sequences.csv",row.names = FALSE)
write.csv(align.intervals, file="align_intervals.csv",row.names = FALSE)

######
## Plot sequences
######
run=FALSE
if (run){
  #Plantid & noise level
  seq.id=levels(testseq.table$`Plant ID`)
  Noise_levels
  #Noise_levels=as.numeric(gsub("Plant_","", seq.id)) #alternative way to retrieve noise levels
  
  #select the noise level
  #   level/id    idx
  #      0         1
  #      0.05      2
  #      0.1       3
  #      0.20      5
  #      0.5       11
  NOISE=11

  #Select sequence depending on noise level:
  test.plot=testseq.table[testseq.table$`Plant ID`==seq.id[NOISE],-1]
  test.plot=cbind.data.frame(intervals=seq(1:nrow(test.plot)),
                             test.plot)
  colnames(test.plot)=colnames(seq)
  test.align=align.intervals[align.intervals$plantID==seq.id[NOISE],-1]

  multiseq_plot(list(seq, test.plot), 
                align.df = test.align, 
                title=paste("Noise level = ", Noise_levels[NOISE]))
}



cat("data generated - end of partI \n")