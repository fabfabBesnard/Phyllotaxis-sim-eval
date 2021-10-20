#!/usr/bin/env Rscript

###########################
## Test script to test the effect of noise
## Part II (before dtw)
###########################
#Created 2021-01-08

setwd("~/Dropbox/Arabidopsis-eval/R_simul-eval/example_data")
source("~/Dropbox/Arabidopsis-eval/R_simul-eval/source/sim_phyllo_sources.R")
source('~/Dropbox/Arabidopsis-eval/R_simul-eval/source/eval_dtw_sources.R')

#####
#1.Import dtw results
####
result.files=list.files(pattern="Plant_")
seq.id=gsub("_result.csv", "", result.files)
dtw.data=read.csv(result.files[1], header=TRUE)
dtw.data=cbind.data.frame(plantID=rep(seq.id[1], nrow(dtw.data)),
                          dtw.data)

if (length(result.files)>1){
  for (i in 2:length(result.files)){
    new.data=read.csv(result.files[i], header=TRUE)
    new.data=cbind.data.frame(plantID=rep(seq.id[i], nrow(new.data)),
                              new.data)
    dtw.data=rbind(dtw.data,new.data)
  }  
}
true.data=read.csv("align_intervals.csv", header=TRUE)

#####
#2. convert data
####
dtw.sub=dtw.data[dtw.data$plantID==seq.id[1],]
sub_dtw_results=convert_dtw_results(dtw.sub[,-1]) #drop the plantID column
dtw_results=cbind.data.frame(plantID=rep(seq.id[1], nrow(sub_dtw_results)), #restore a plantID column in the new df
                                 sub_dtw_results)

if (length(seq.id)>1){
  for (i in 2:length(seq.id)){
    dtw.sub=dtw.data[dtw.data$plantID==seq.id[i],]
    sub_dtw_results=convert_dtw_results(dtw.sub[,-1])
    sub_dtw_results=cbind.data.frame(plantID=rep(seq.id[i], nrow(sub_dtw_results)), 
                                     sub_dtw_results)
    dtw_results=rbind(dtw_results, sub_dtw_results)
  }
}

#####
#3. Assess data
####
eval.subtab=compare_align(dtw_results[dtw_results$plantID==seq.id[1], -1], 
                   true.data[true.data$plantID==seq.id[1],-1])
eval.table=cbind.data.frame(plantID=rep(seq.id[1], nrow(eval.subtab)),
                            eval.subtab)
eval.subsum=eval_summarize(true.data[true.data$plantID==seq.id[1],-1], eval.subtab)
eval.summary=cbind.data.frame(plantID=seq.id[1],
                              eval.subsum)

if (length(seq.id)>1){
  for (i in 2:length(seq.id)){
    #print(paste("loop is at step n°", i))
    #print(seq.id[i])
    #full table
    eval.subtab=compare_align(dtw_results[dtw_results$plantID==seq.id[i], -1], 
                           true.data[true.data$plantID==seq.id[i],-1])
    eval.subtab=cbind.data.frame(plantID=rep(seq.id[i], nrow(eval.subtab)),
                                eval.subtab)
    eval.table=rbind(eval.table, eval.subtab)
    #summary
    eval.subsum=eval_summarize(true.data[true.data$plantID==seq.id[i],-1], eval.subtab)
    eval.subsum=cbind.data.frame(plantID=seq.id[i],
                                  eval.subsum)
    eval.summary=rbind(eval.summary, eval.subsum)
  }
}

#####
#4. plots
####
plots=FALSE
if (plots){
  require(ggplot2)
  eval.summary$noise.level=as.numeric(gsub("Plant_","", seq.id))
  ggplot(data=eval.summary, aes(x=noise.level, y=Correct.ratio))+
    geom_line()+geom_point()+
    ylim(0,1)
  
  plant.group1=seq.id[seq.id %in% eval.summary[eval.summary$Correct.ratio==1,]$plantID]
  plant.group2=seq.id[! seq.id %in% eval.summary[eval.summary$Correct.ratio==1,]$plantID]
  plant.group3=seq.id[6:13]
  ggplot(data=eval.table[eval.table$plantID %in% plant.group1,], aes(dtw.cost, fill=eval))+
    geom_histogram(position = "dodge", binwidth = 0.05, color="black")+
    facet_wrap(~plantID)
  ggplot(data=eval.table[eval.table$plantID %in% plant.group3,], aes(dtw.cost, fill=eval))+
    geom_histogram(position = "dodge", binwidth = 0.05)+
    facet_wrap(~plantID)
}

