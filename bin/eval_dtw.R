#!/usr/bin/env Rscript

###########################################################################
#### Copyright (C) 2021 - INRAe (Fabrice Besnard, RDP)
#### This script is a free software: you can redistribute it
#### and/or modify it under the terms of the GNU General Public
#### License (GNU GPL) as published by the Free Software Foundation, either
#### version 3 of the License, or (at your option) any later version.
#### Distributed without any warranty.
###########################################################################
# started: 2021-01-05
# last edit: 2021-11-04
#Version v0

###############
##   Usage   ##
###############
#Rscript eval_dtw.R -a dtw_res.csv

################################
####   INPUTS / Arguments    ###
################################
require("optparse")
option_list = list(
  make_option(c("-a", "--alignment_dtw"), type="character", default=NULL, 
              help="dtw result file name", metavar="character"),
  make_option(c("-r", "--reference_seq"), type="character", default=NULL, 
              help="dataset file name", metavar="character"),
  make_option(c("-t", "--test_seq"), type="character", default=NULL, 
              help="dataset file name", metavar="character"),
  make_option(c("-i", "--intervals_truealign"), type="character", default=NULL, 
              help="dataset file name", metavar="character"),
  make_option(c("-o", "--organs_truealign"), type="character", default=NULL, 
              help="dataset file name", metavar="character"),
  make_option(c("-R", "--repository"), type="character", default="~/Dropbox/Arabidopsis-eval/Phyllotaxis-sim-eval/", 
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

#lines for debug (run from Rconsole)
setwd("~/Dropbox/Arabidopsis-eval/Phyllotaxis-sim-eval/example_data/Notebook_tests")
opt=list()
opt$repository="~/Dropbox/Arabidopsis-eval/Phyllotaxis-sim-eval"
opt$alignment_dtw="ploufplouf_result.csv"
opt$reference_seq="reference_sequences.csv"
opt$test_seq="test_sequences.csv"
opt$intervals_truealign="align_intervals.csv"
opt$organs_truealign="align_organs.csv"
opt$verbose=TRUE

####################
#### Body Run  #####
####################
#up-load libraries
local.repo=opt$repository #must end by '/Phyllotaxis-sim-eval/'
if (!grepl("/$", local.repo)){#add an ending / if missing
  local.repo=paste0(local.repo, "/")
}
source(paste0(local.repo, "source/eval_dtw_sources.R"))
source(paste0(local.repo, "source/plot_sequences_sources.R"))

#up-load input data:
raw.results=read.csv(opt$alignment_dtw, header=TRUE)
seqs=read.csv(opt$reference_seq, col.names = c("PlantID", "angles", "internodes"))
tests=read.csv(opt$test_seq, col.names = c("PlantID", "angles", "internodes"))
Ialign=read.csv(opt$intervals_truealign)

#Convert raw results from dtw:
dtw_results=convert_dtw_results(raw.results, 
                                seq.ref = seqs, seq.test = tests, 
                                verbose=opt$verbose)

#Compare dtw alignment with the ground truth alignment
prediction_eval=evaluate_align_prediction(dtw_results = dtw_results, 
                                          true_align = Ialign, 
                                          verbose = opt$verbose)
#Summarize the assessment by plant
summarize_prediction_eval(comparison.df = prediction_eval, 
                          true_align.df = Ialign, 
                          verbose=opt$verbose)
