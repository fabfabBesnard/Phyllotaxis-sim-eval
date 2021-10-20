#!/usr/bin/env Rscript

###########################################################################
#### Copyright (C) 2020 - INRAe (Fabrice Besnard, RDP)
#### This script is a free software: you can redistribute it
#### and/or modify it under the terms of the GNU General Public
#### License (GNU GPL) as published by the Free Software Foundation, either
#### version 3 of the License, or (at your option) any later version.
#### Distributed without any warranty.
###########################################################################
#started 2020-12
# last edit: 2020-03-05
#Version v0

###############
##   Usage   ##
###############
#Rscript test_dtw.R -f input_table.csv -p
Program_Description="
  Script test_dtw.R \n
Description: \n
This scripts generates
-a raw artifical pyllotaxis sequences (angles + internodes) \n
-a derived perturbed sequences because of segmentation errors, given by an input file \n
-the alignment/correspondance between the raw and derived sequences (as tables and optionnaly plots) \n
The input file contains the segmentation error to introduce in the alignment. 
It is a table, on plant per row, precising the total number of intervals, id of organs gains and ids of organs lost. \n
"

################################
####   INPUTS / Arguments    ###
################################
suppressPackageStartupMessages(require(optparse))

option_list = list(
  make_option(c("-f", "--file"), type="character", default=NULL, 
              help="input table indicating the sequences to simulate (N and segmentation errors)", metavar="character"),
  make_option(c("-e", "--permutation"),type="integer", default=0,
              help="Likelihood that 2 close organs will be permuted (range [0,1]). 0 deactivate the function. [default %default]", metavar="number"),
  make_option(c("-c", "--close_organs"),type="integer", default=2,
              help="Internode value below which the corresponding pair of organs can be permuted. [default %default]", metavar="number"),
  make_option(c("-n", "--noplots"), action="store_true", default=TRUE,
              help="do not print plots [default]"),
  make_option(c("-p", "--plots"), action="store_false", 
              dest="noplots", help="Print plots")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if (is.null(opt$file)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}

if (opt$help){
  cat(Program_Description)
}

#lines for debug (uncomment to run this script from Rconsole)
# setwd("~/Dropbox/Arabidopsis-eval/R_simul-eval/example_data")
# data=read.delim2("dtw_input_plants.csv", header=TRUE)
# opt=list()
# opt$permutation=1
# opt$close_organs=2
# opt$noplots=FALSE

##################
##  PARAMETERS  ##
##################
###I. PHYLLOTAXIS
#####
## Divergence angles
#Canonical angle
alpha=137.5
#angle_sd (correspond to real biological variation)
a_sd=18.5
#-> cf Guedon et al. JTB 2013: standard deviation a_sd=18.5

## Internodes
#Internode_noise: gaussian noise: mean=0, sd=i_noise
i_Gsd=0.8
#ratio of the biological noise/variation compared to the value of the internode, expressed in pct
i_noise_pct=75

####################
#### Body Run  #####
####################
source("~/Dropbox/Arabidopsis-eval/R_simul-eval/source/sim_phyllo_sources.R")
setwd("~/Dropbox/Arabidopsis-eval/R_simul-eval/example_data")

data=read.delim2(opt$file, header=TRUE)

#Initialize output tables:
refseq=data.frame(PlantID=NULL,
                  angles=NULL,
                  internodes=NULL)
testseq=data.frame(PlantID=NULL,
                   angles=NULL,
                   internodes=NULL)
align.intervals=data.frame(PlantID=NULL,
                           reference=NULL,
                           modified=NULL,
                           dtw=NULL)
align.organs=data.frame(PlantID=NULL,
                        reference=NULL,
                        modified=NULL,
                        segmentation=NULL)

#Loop over plant in the tables to generate the data
for (i in 1:nrow(data)){
  print(paste("processing data for plant", data[i,]$PlantID ))
  N=data[i,]$N_interval
  
  #Generate data
  #make a reference sequence
  ref=make_refseq(N, alpha, a_sd, i_Gsd, i_noise_pct)
  #make default alignments (organs / intervals)
  align.list=make_align_list(N)
  #introduce segmentation errors
  gains=as.numeric(unlist(strsplit(as.character(data[i,]$organ_gain), ",")))
  losses=as.numeric(unlist(strsplit(as.character(data[i,]$organ_loss), ",")))
  test=segmentation_errors(ref, align.list = align.list,
                           organ_gain=gains,
                           organ_loss=losses)
  #If measure-permutations are activated, run permutation function:
  if (opt$permutation > 0){
    test=simple_measure_permutation(test$values, align.list = list(test$I, test$O), 
                                    i_threshold = opt$close_organs, proba=opt$permutation, verbose = FALSE)
  }
  
  
  #Print the alignment if asked by user
  if (!opt$noplots) {
    multiseq_plot(list(ref, test$values), 
                  align.df = test$I,
                  id.names = c("ref", "test"), 
                  title=as.character(data[i,]$PlantID),
                  ref.first = TRUE)
  }
  
  #Fill up output tables
  refseq=rbind.data.frame(refseq,
                          cbind.data.frame(PlantID=rep(data[i,]$PlantID, nrow(ref)), 
                                           ref[,c(2,3)])) #drop the "intervals" column
  testseq=rbind.data.frame(testseq,
                           cbind.data.frame(PlantID=rep(data[i,]$PlantID, nrow(test$values)), 
                                            test$values[,c(2,3)])) #drop the "intervals" column
  align.intervals=rbind.data.frame(align.intervals,
                                   cbind.data.frame(PlantID=rep(data[i,]$PlantID, nrow(test$I)), 
                                                    test$I))
  align.organs=rbind.data.frame(align.organs,
                                   cbind.data.frame(PlantID=rep(data[i,]$PlantID, nrow(test$O)), 
                                                    test$O))
}

#Write out the data
colnames(refseq)=c("PlantID", "angles", "Internodes")
write.csv(refseq, file="reference_sequences.csv", row.names = FALSE)
colnames(testseq)=c("PlantID", "angles", "Internodes")
write.csv(testseq, file="test_sequences.csv",row.names = FALSE)
write.csv(align.intervals, file="align_intervals.csv",row.names = FALSE)
write.csv(align.organs, file="align_organs.csv",row.names = FALSE)

cat("data generated - end of script \n")
