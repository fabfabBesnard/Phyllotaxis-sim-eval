#!/usr/bin/env Rscript

###########################################################################
#### Copyright (C) 2021 - INRAe (Fabrice Besnard, RDP)
#### This script is a free software: you can redistribute it
#### and/or modify it under the terms of the GNU General Public
#### License (GNU GPL) as published by the Free Software Foundation, either
#### version 3 of the License, or (at your option) any later version.
#### Distributed without any warranty.
###########################################################################
#started 2021-04
# last edit: 2021-11-08
#Version v0

###############
##   Usage   ##
###############
#Rscript simul_data.R -f input_table.csv -p -o data1
# mandatory options:
            # -f (--file) table file (.csv) describing the properties of the paired sequences to generate
# + facultative options:
            # -p (--plots): print plots
            # -o (--output_prefix): prefix for all outputs
            # -s (--setseed): provides a seed for the random values generated in simulated sequences
            # -R (--repository) path/to/Phyllotaxis-sim-eval/ (default is ~/Dropbox/Arabidopsis-eval/Phyllotaxis-sim-eval/)
            # -D (--destination) path/to/dest (default is current working directory)
            # -v (--verbose)

Program_Description="
  Script simul_data.R \n
Description: \n
The program ‘simul_data.R’ generates pairs of simulated phyllotaxis sequences representing a reference sequence and a test sequence
Divergences between ref & test come either from segmentation errors, permutations or noise on values.
Two different designs can be simulated: either a simple divergence from an initial 'reference' sequence or a comparison of two parallel measurements (e.g. manual/computer) 

The program ouputs 4 tables (.csv format) and (optionaly) 1 plot file (.pdf):
- reference_sequences.csv : table with the values of angles and internodes of the reference sequence
- test_sequences.csv : table with the values of angles and internodes of the test sequence
- align_intervals.csv : table aligning each intervals of the reference and test sequences for all plants, with indications of the nature of the match/mismatch
- align_organs.csv : table aligning each intervals of the reference and test sequences for all plants, with indications of the nature of the match/mismatch
- Rplots.pdf : (optional) plots of the alignment of the test sequence against the reference sequence

The input file details the segmentation errors that will be introduced in the alignment. 
Its format is fixed and must be respected for the program to run correctly.
It is a table with one plant per row, precising the characteristics of the reference sequence and the scenario and parameters that modify it into the test sequence.
Please refer to the readme section of the template input file.
"

################################
####   INPUTS / Arguments    ###
################################
suppressPackageStartupMessages(require(optparse))

option_list = list(
  make_option(c("-f", "--file"), type="character", default=NULL, 
              help="input table that configures the sequences to simulate (ex: Number, length, segmentation errors, ...)", metavar="character"),
  make_option(c("-n", "--noplots"), action="store_true", default=TRUE,
              help="do not print plots [default]"),
  make_option(c("-p", "--plots"), action="store_false", 
              dest="noplots", help="Print plots"),
  make_option(c("-o", "--output_prefix"), type="character", default=NULL, 
              help="prefix for all outputs", metavar="character"),
  make_option(c("-s", "--setseed"), type="numeric", default=NULL, 
              help="give a number for seed", metavar="numeric"),
  make_option(c("-D", "--destination"), type="character", default=NULL, 
             help="destination folder", metavar="character"),
  make_option(c("-R", "--repository"), type="character", default="~/Phyllotaxis-sim-eval/", 
              help="local path to 'Phyllotaxis-sim-eval' repository", metavar="character"),
  make_option(c("-v", "--verbose"), action="store_true", default=FALSE,
              help="increase verbosity")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if (is.null(opt$file)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file).", call.=FALSE)
}

if (opt$help){
  cat(Program_Description)
}

## lines for Rconsole debug (uncomment to run this script from Rconsole)
# # setwd("~/Documents/RDP/MyProjects/ROMI/Data/Eval_AnglesAndInternodes/experiments/Exp1/")
# setwd("~/Documents/RDP/MyProjects/ROMI/Data/Eval_AnglesAndInternodes/Phyllotaxis-sim-eval/example_data/Notebook_tests/")
# opt=list()
# #opt$file="E1_INPUT_test-noiselevels.csv"
# opt$file="simulation_plants_nb.csv"
# opt$noplots=FALSE
# opt$repository="~/Documents/RDP/MyProjects/ROMI/Data/Eval_AnglesAndInternodes/Phyllotaxis-sim-eval/"
# opt$destination="~/Documents/RDP/MyProjects/ROMI/Data/Eval_AnglesAndInternodes/tests/"
# opt$output_prefix="debug"
# opt$setseed=NULL
# opt$verbose=TRUE

#############################
##  Hard-coded PARAMETERS  ##
#############################
###I. PHYLLOTAXIS
#####
## Divergence angles
#Canonical angle:
alpha=137.5
#angle_sd (correspond to real biological variation)
a_sd=18.5 #-> cf Guedon et al. JTB 2013: standard deviation a_sd=18.5

## Internodes
#Internode_noise: gaussian noise: mean=0, sd=i_noise
i_Gsd=0.8
#ratio of the biological noise/variation compared to the value of the internode, expressed in pct
i_noise_pct=75

cat("Starting script to simulate paired sequences of phyllotaxis \n")
########################
## up-load input data ##
########################
data=read.delim(opt$file, header=TRUE)

#################################
#### Set-up in/out options  #####
#################################
#Path to local code repository
local.repo=opt$repository #must end by '/Phyllotaxis-sim-eval/'
if (!grepl("/$", local.repo)){#add an ending / if missing
  local.repo=paste0(local.repo, "/")
}
source(paste0(local.repo, "source/sim_phyllo_sources.R"))
source(paste0(local.repo, "source/plot_sequences_sources.R"))

#Set-up the destination folder for the outputs
if (is.null(opt$destination)){ opt$destination=getwd() }
setwd(opt$destination)

####################
#### Body Run  #####
####################
#reminder: expected fields in the input config table:
# PlantID
# N_interval
# organ_gain
# organ_loss
# Noise_or_Measures
# measure1_angle_sd
# measure1_internode_sd
# measure2_angle_sd
# measure2_internode_sd
# sd_noise_level
# sd_noise_scale
# mean_noise_bias
# permutation_length
# permutation_likelihood
# Note

if (!is.null(opt$setseed)){
  if (opt$verbose){ print(paste("A seed will be given for the simulation: seed =", opt$setseed))}
  set.seed(opt$setseed)
}

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

convert_tuple=function(input) {
  #force conversion of "" to NA
  input=ifelse(input=="", 0, input)
  if (!is.numeric(input)){
    #remove brackets
    temp=gsub(")","",gsub("(", "", input, fixed=TRUE), fixed = TRUE)
    output=as.vector(sapply(strsplit(temp, ",")[[1]], as.numeric))
  } else (output=input)
  return(output)
}

#LOOP over plants (= rows) in the tables to generate the data
for (i in 1:nrow(data)){
  if (opt$verbose){print(paste("processing data for plant", data[i,]$PlantID ))}
  
  #####
  # Retrieve data from the input file and build the scenario (true/false for seg_errors, measure, noise, permutation)
  #####
  N=data[i,]$N_interval #Nber of intervals
  # Segmentation errors
  ## Gains
  if(length(grep("random",data[i,]$organ_gain))>0){
    nb_gain=as.numeric(unlist(strsplit(as.character(data[i,]$organ_gain), ","))[2]) # get the nber of organs to insert randomly
    GAIN=round(runif(nb_gain,min=1, max=(N+2)), digits = 0)
    GAIN=GAIN[order(GAIN)]
  } else {
    GAIN=as.numeric(unlist(strsplit(as.character(data[i,]$organ_gain), ",")))
  }
  if (length(GAIN)==0 || is.na(GAIN)){GAIN=NULL}
  ## Losses
  if(length(grep("random",data[i,]$organ_loss))>0){
    nb_loss=as.numeric(unlist(strsplit(as.character(data[i,]$organ_loss), ","))[2]) # get the nber of organs to insert randomly
    LOSS=round(runif(nb_loss,min=1, max=(N+1)), digits = 0)
    LOSS=LOSS[order(LOSS)]
  } else {
    LOSS=as.numeric(unlist(strsplit(as.character(data[i,]$organ_loss), ",")))
  }
  if (length(LOSS)==0 || is.na(LOSS)){LOSS=NULL}
  if (is.null(GAIN) && is.null(LOSS)){seg_errors=FALSE} else {seg_errors=TRUE}
  
  # two measures versus noise
  NvsM=data[i,]$Noise_or_Measures
  NvsM=tolower(NvsM) #force conversion to lower case
  if (NvsM != "measures" & NvsM != "noise" ){
    stop(paste("for plant", data[i,]$PlantID, ": error in 'Noise_or_Measures' content: please choose among 'measures' or 'noise' only."))}
  if (NvsM == "measures"){
    manual_anoise_sd=data[i,]$measure1_angle_sd
    manual_inoise_sd=data[i,]$measure1_internode_sd
    aut_anoise_sd=data[i,]$measure2_angle_sd
    aut_inoise_sd=data[i,]$measure2_internode_sd } 
  if (NvsM == "noise") {
    sd_noise_level = convert_tuple(data[i,]$sd_noise_level)
    sd_noise_scale = ifelse(data[i,]$sd_noise_scale=="", "mean", data[i,]$sd_noise_scale) #default value is "mean" if the field is empty
    mean_noise_bias = convert_tuple(data[i,]$mean_noise_bias)
  }
  
  # Permutations
  permut_length=data[i,]$permutation_length
  if (is.na(permut_length)){permut_length=0}
  permut_proba=data[i,]$permutation_likelihood
  if (is.na(permut_proba)){permut_proba=0}
  if(permut_length==0 || permut_proba==0){permutation=FALSE} else {permutation=TRUE}
  
  #####
  # Print scenario
  #####
  if (opt$verbose) { print_info(seg_errors=seg_errors, permutation = permutation, Noise_or_Measures = NvsM) }

  #####
  #Run the scenario and generate corresponding data
  #####
  #Start by generating a sequence of angles and internodes with the same length desired for the reference sequence
  seq=make_refseq(N, alpha, a_sd, i_Gsd, i_noise_pct)
  default.align=make_align_list(N) #create a default alignment with only matches between ref and test sequence
  if (NvsM == "measures"){
    #Modify seq by two independent measures:
    seq.ref=make_measure(seq, anoise_sd = manual_anoise_sd, inoise_sd = manual_inoise_sd,
                         noise.scale = "absolute", verbose=opt$verbose)
    seq.aut=make_measure(seq, anoise_sd = aut_anoise_sd, inoise_sd = aut_inoise_sd,
                         noise.scale = "absolute", verbose=opt$verbose)
    
    if (seg_errors){ #segmentation errors only affect the test sequence, not the reference sequence
      seq.test=segmentation_errors(seq.aut,default.align,
                                   organ_gain=GAIN,
                                   organ_loss=LOSS) }
    else {
      seq.test=list(seq.aut, default.align$Ialign, default.align$Oalign) 
      names(seq.test)=c("values", "I", "O") }
    
    if (permutation){
      seq.test=simple_measure_permutation(seq.test$values, align.list = list(seq.test$I, seq.test$O),
                                          i_threshold = permut_length, proba = permut_proba, verbose = opt$verbose )
    }
    
  } else if (NvsM == "noise") {
    seq.ref=seq
    if (length(sd_noise_level) == 1 ) {
      sd_noise_level=ifelse(is.na(sd_noise_level),0,sd_noise_level) #default value is (0,0) if the field is empty
      #duplicate the value for angles, internodes, respectively
      sd_noise_level=c(sd_noise_level, sd_noise_level)
    }
    if (length(mean_noise_bias) == 1 ) {
      mean_noise_bias=ifelse(is.na(mean_noise_bias),0,mean_noise_bias) #default value is (0,0) if the field is empty
      #duplicate the value for angles, internodes, respectively
      mean_noise_bias=c(mean_noise_bias, mean_noise_bias)
    }
    #Apply noise:
    seq.noise=make_measure(seq, anoise_sd = sd_noise_level[1], inoise_sd = sd_noise_level[2],
                           noise.scale = sd_noise_scale, 
                           anoise.mean = mean_noise_bias[1], inoise.mean = mean_noise_bias[2], 
                           verbose=opt$verbose)
    
    if (seg_errors){
      seq.test=segmentation_errors(seq.noise,default.align,
                                   organ_gain=GAIN,
                                   organ_loss=LOSS) }
    else {
      seq.test=list(seq.noise, default.align$Ialign, default.align$Oalign)
      names(seq.test)=c("values", "I", "O")}
    
    if (permutation){
      seq.test=simple_measure_permutation(seq.test$values, align.list = list(seq.test$I, seq.test$O),
                                          i_threshold = permut_length, proba = permut_proba, verbose = opt$verbose )
    }
  } else {
    #Without noise nor measure  
    seq.ref=seq
    
    if (seg_errors){
      seq.test=segmentation_errors(seq.ref,default.align,
                                   organ_gain=GAIN,
                                   organ_loss=LOSS) }
    else {
      seq.test=list(seq.ref, default.align$Ialign, default.align$Oalign) 
      names(seq.test)=c("values", "I", "O") }
    
    if (permutation){
      seq.test=simple_measure_permutation(seq.test$values, align.list = list(seq.test$I, seq.test$O),
                                          i_threshold = permut_length, proba = permut_proba, verbose = opt$verbose )
    }
  }
  
  #####
  #Fill up output tables
  #####
  refseq=rbind.data.frame(refseq,
                          cbind.data.frame(PlantID=rep(data[i,]$PlantID, nrow(seq.ref)), 
                                           seq.ref[,c(2,3)])) #drop the "intervals" column
  testseq=rbind.data.frame(testseq,
                           cbind.data.frame(PlantID=rep(data[i,]$PlantID, nrow(seq.test$values)), 
                                            seq.test$values[,c(2,3)])) #drop the "intervals" column
  align.intervals=rbind.data.frame(align.intervals,
                                   cbind.data.frame(PlantID=rep(data[i,]$PlantID, nrow(seq.test$I)), 
                                                    seq.test$I))
  align.organs=rbind.data.frame(align.organs,
                                   cbind.data.frame(PlantID=rep(data[i,]$PlantID, nrow(seq.test$O)), 
                                                    seq.test$O))
}

########################
## Write out the data ##
########################
if (is.null(opt$output_prefix)){
  opt$output_prefix="" } else { opt$output_prefix=paste0(opt$output_prefix,"_")}

# plots:
if (!opt$noplots) {
  multiseq_plot_pdf(seq.ref=refseq, seq.test=testseq, 
                true.align=align.intervals,
                id.names=c("reference", "test"),
                pdf.name=paste0(opt$output_prefix,"SimulatedPairedSequences.pdf"))
}

# Tables:
#Note: destination folder has been set-up above in the loop
colnames(refseq)=c("PlantID", "angles", "Internodes")
write.csv(refseq, file=paste0(opt$output_prefix,"reference_sequences.csv"), row.names = FALSE)
colnames(testseq)=c("PlantID", "angles", "Internodes")
write.csv(testseq, file=paste0(opt$output_prefix,"test_sequences.csv"),row.names = FALSE)
write.csv(align.intervals, file=paste0(opt$output_prefix,"align_intervals.csv"),row.names = FALSE)
write.csv(align.organs, file=paste0(opt$output_prefix,"align_organs.csv"),row.names = FALSE)

cat("simulated data generated - end of script \n")
