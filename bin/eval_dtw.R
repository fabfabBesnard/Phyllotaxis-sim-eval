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
# last edit: 2021-11-05
#Version v0

###############
##   Usage   ##
###############
#Requirements: the four input csv files.
#Rscript eval_dtw.R -a dtw_res.csv -r reference_seq.csv -t test_seq.csv -i true_interval_align.csv 

#Compulsory options:
                  # -a (--alignment_dtw)
                  # -r (--reference_seq)
                  # -t (--test_seq)
                  # -i (--intervals_truealign)
# + facultative options:
#                 -p (--plots): print plots
#                 -d (--detail): output detailed table
#                 -o (--output_prefix): prefix for all outputs
#                 -R (--repository) path/to/Phyllotaxis-sim-eval/ (default is ~/Dropbox/Arabidopsis-eval/Phyllotaxis-sim-eval/)
#                 -D (--destination) path/to/dest (default is current working directory)
#                 --verbose

cat("Starting script to evaluate dtw alignment prediction \n")

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
  make_option(c("-n", "--noplots"), action="store_true", default=TRUE,
              help="do not print plots [default]"),
  make_option(c("-p", "--plots"), action="store_false", 
              dest="noplots", help="Print plots"),
  make_option(c("-d", "--detail"), action="store_true", default=FALSE,
              help="outputs the detailed evaluation interval by interval"),
  make_option(c("-o", "--output_prefix"), type="character", default=NULL, 
              help="prefix for all outputs", metavar="character"),
  make_option(c("-D", "--destination"), type="character", default=NULL, 
              help="destination folder", metavar="character"),
  make_option(c("-R", "--repository"), type="character", default="~/Dropbox/Arabidopsis-eval/Phyllotaxis-sim-eval/", 
              help="local path to 'Phyllotaxis-sim-eval' repository", metavar="character"),
  make_option(c("-v", "--verbose"), action="store_true", default=FALSE,
              help="increase verbosity")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

#lines for debug (run from Rconsole)
# setwd("~/Dropbox/Arabidopsis-eval/Phyllotaxis-sim-eval/example_data/Notebook_tests")
# opt=list()
# opt$repository="~/Dropbox/Arabidopsis-eval/Phyllotaxis-sim-eval"
# opt$alignment_dtw="ploufplouf_result.csv"
# opt$reference_seq="reference_sequences.csv"
# opt$test_seq="test_sequences.csv"
# opt$intervals_truealign="align_intervals.csv"
# opt$noplots=FALSE
# opt$detail=TRUE
# opt$output_prefix=NULL
# opt$verbose=TRUE

########################
## up-load input data ##
########################
raw.results=read.csv(opt$alignment_dtw, header=TRUE)
seqs=read.csv(opt$reference_seq, col.names = c("PlantID", "angles", "internodes"))
tests=read.csv(opt$test_seq, col.names = c("PlantID", "angles", "internodes"))
Ialign=read.csv(opt$intervals_truealign)

#################################
#### Set-up in/out options  #####
#################################
#Path to local code and libraries repository
local.repo=opt$repository #must end by '/Phyllotaxis-sim-eval/'
if (!grepl("/$", local.repo)){#add an ending / if missing
  local.repo=paste0(local.repo, "/")
}
source(paste0(local.repo, "source/eval_dtw_sources.R"))
source(paste0(local.repo, "source/plot_sequences_sources.R"))

#Set-up the destination folder for the outputs
if (is.null(opt$destination)){ opt$destination=getwd() }
setwd(opt$destination)

####################
#### Body Run  #####
####################
#Convert raw results from dtw:
dtw_results=convert_dtw_results(raw.results, 
                                seq.ref = seqs, seq.test = tests, 
                                verbose=opt$verbose)

#Compare dtw alignment with the ground truth alignment
prediction_eval=evaluate_align_prediction(dtw_results = dtw_results, 
                                          true_align = Ialign, 
                                          verbose = opt$verbose)
#Summarize the assessment by plant
summary=summarize_prediction_eval(comparison.df = prediction_eval, 
                                  true_align.df = Ialign, 
                                  verbose=opt$verbose)

#Write out the data
if (is.null(opt$output_prefix)){
  opt$output_prefix="" } else { opt$output_prefix=paste0(opt$output_prefix,"_")}
write.csv(summary, file=paste0(opt$output_prefix,"PredictionEval_summary.csv"), row.names = FALSE)
if (opt$detail){
  write.csv(prediction_eval, file=paste0(opt$output_prefix,"PredictionEval_detail.csv"), row.names = FALSE)
}
if (!opt$noplots) {
  compare_plots(seq.ref=seqs, seq.test=tests, 
                true.align=Ialign, dtw.results = dtw_results, prediction.eval= prediction_eval,
                id.names=c("reference", "test"),
                PDF=TRUE, pdf.name=paste0(opt$output_prefix,"PredictionEval_plots.pdf"))
}

cat("data generated - end of script \n")