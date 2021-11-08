#!/bin/bash

###########################################################################
#### Copyright (C) 2021 - INRAe (Fabrice Besnard, RDP)
#### This script is a free software: you can redistribute it
#### and/or modify it under the terms of the GNU General Public
#### License (GNU GPL) as published by the Free Software Foundation, either
#### version 3 of the License, or (at your option) any later version.
#### Distributed without any warranty.
###########################################################################
# created: 2021-11-08
version="v0" #version

###########
#  USAGE  #
###########
usage() { 
echo "Script 'simul_predict_and_eval.sh'; version $version"
echo "Usage: $0 [--file <PATH>] [OPTIONS]"
echo "OPTIONS:"
echo "	-f, --file (mandatory)
	Path to the table file (.csv) describing the properties of the paired sequences to generate. Use absolute paths.
"
echo "	-R, --repository (optional)
	Path to the local repository of 'Phyllotaxis-sim-eval'. Note that dtw must be installed through a conda environment.
"
echo "	-D, --destination (optional)
	Path to destination folder for outputs. Use absolute paths. If the folder does not exist, the program will create it at the specified path.
"
echo "	-o, --output_prefix (optional)
	prefix for all outputs files
"
echo "	-p, --plots (optional)
	generate plots
"
echo "	-d, --detail (optional)
	provides a detailed table of the assessment for each interval of the test sequences.
"
echo "	-h, --help: 
	print this help and exit
"
echo "  -v, --verbose: 
	increase verbosity
"
exit 1;}

# Running example
# local=path/to/local/Phyllotaxis-sim-eval
# dest=path/to/output/folder

# bash $local/bin/simul_predict_and_eval.sh \
# --file $local/example_data/Notebook_tests/simulation_plants_nb.csv \
# -R $local \
# -D $dest \
# -o testbash \
# --detail --plots --verbose 

####################
#  Local functions #
####################

function print_verbose {
#funtion to print info only if verbose mode is active. I uses the global $VERBOSE variable defined in options
   local MESSAGE="${@}"
   if [[ "${VERBOSE}" == true ]];then
      echo "${MESSAGE}"
   fi
}

##############################
#  Read & check User inputs  #
##############################

options_step1=""
options_step3=""

while [ "$1" != "" ]; do
  case $1 in
  -f | --file)
    shift
    file=$1
    ;;
  -R | --repository)
    shift
    localrepo=$1
    ;;
  -D | --destination)
    shift
    dest=$1
    ;;
  -o | --output_prefix)
    shift
    output_prefix=$1
    ;;
  -p | --plots)
    plots='true'
    options_step3+=" --plots"
    ;;
  -d | --detail)
    detail_table='true'
    options_step3+=" --detail"
    ;;
  -v | --verbose)
    VERBOSE='true'
    options_step1+=" --verbose"
    options_step3+=" --verbose"
    print_verbose "Verbose mode is ON"
    ;;
  -h | --help)
    usage
    exit
    ;;
  *)
    usage
    exit 1
    ;;
  esac
  shift
done

if [ -z $file ] ; then
	echo "input file -f (--file) is mandatory"
	usage
fi

################
#  script body #
################

if [ ! -d $dest ]; then
	echo "destination folder does not exist yet, so the program will create it."
	mkdir $dest
fi

## Step1: simulate data
print_verbose "Step1: Simulating paired sequences of phyllotaxis from input file"
Rscript $localrepo/bin/simul_data.R \
	--repository $localrepo \
	--file $file \
	--destination $dest \
	--output_prefix $output_prefix \
	$options_step1

## Step2: Predict an alignment for the paired sequences using sm-dtw
print_verbose "Step2: Predict an alignment for the paired sequences using sm-dtw"
source ~/softwares/miniconda3/bin/activate
conda activate romi

align_csv_database.py $dest/${output_prefix}_reference_sequences.csv $dest/${output_prefix}_test_sequences.csv \
$output_prefix --free_ends 0.4

## Step3: Assess the alignment prediction made by sm-dtw
print_verbose "Step3: Assess the alignment prediction made by sm-dtw"

Rscript $localrepo/bin/eval_dtw.R \
--repository $localrepo \
--destination $dest \
--alignment_dtw $dest/${output_prefix}_result.csv \
--reference_seq $dest/${output_prefix}_reference_sequences.csv \
--test_seq $dest/${output_prefix}_test_sequences.csv \
--intervals_truealign $dest/${output_prefix}_align_intervals.csv \
--output_prefix $output_prefix \
$options_step3

##Cleaning
if [[ $plots == 'true' ]]; then
	rm $dest/Rplots.pdf
fi

##end of script
print_verbose "End of script 'simul_predict_and_eval.sh'"
