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
# last edit: 2021-01-07
#Version v0

###############
##   Usage   ##
###############
#Rscript eval_dtw.R -f dtw_res.csv

################################
####   INPUTS / Arguments    ###
################################
require("optparse")
option_list = list(
  make_option(c("-f", "--file"), type="character", default=NULL, 
              help="dataset file name", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if (is.null(opt$file)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}
data=read.csv(opt$file, header=TRUE)

#lines for debug (run from Rconsole)
setwd("~/Dropbox/Arabidopsis-eval/R_simul-eval/example_data")
data=read.csv("Plant#2_result.csv", header=TRUE)

