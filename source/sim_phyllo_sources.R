###########################################################################
#### Copyright (C) 2020 - INRAe (Fabrice Besnard, RDP)
#### This script is a free software: you can redistribute it
#### and/or modify it under the terms of the GNU General Public
#### License (GNU GPL) as published by the Free Software Foundation, either
#### version 3 of the License, or (at your option) any later version.
#### Distributed without any warranty.
###########################################################################
#started 2020-09-20
# last edit: 2021-01-13
#Version v0

##Content:
#Custom functions to generate synthetic phyllotaxis sequences that mimic different real situations:
#-real biological sequences with biological variations (like natural variance of different parameters and permutations)
#-measured sequences without errors
#-measured sequences with errors of segmentations

##################
##   Libraries  ##
##################
suppressPackageStartupMessages(require(ggplot2))
suppressPackageStartupMessages(require(reshape2))
suppressPackageStartupMessages(require(gridExtra))

###################################
## Generate Reference sequences  ##
###################################
make_refseq=function(N, #length of the sequence
                     alpha=137.5, a_sd=30, #canonical angle value / sd of angles (Gaussian dist.)
                     natural.permutation=TRUE, permutation.frequency=0.04,
                     i_Gsd=1.5, i_noise_pct=75, 
                     i_beta=1.5, i_max=100, i_plateau=5,
                     verbose=FALSE){
  #DESCRIPTION: it generates typical "Arabidopsis" phyllotaxis data made of two sequences: 
  #               - the sequence of successive divergence angle 
  #               - and of successive internodes
  #             from the first cauline branch (base of the raceme) to the top of the inflorescence
  #Default values have been optimized to simulate realistic data
  # [output]: a dataframe with three columns: $intervals, $angles and $internodes
  
  #1. Generate angles (only integer values)
  angles=round(rnorm(N, mean=alpha, sd=a_sd), digits = 0) %% 360
  
  #2. Generate internodes (ony integer values)
  # parameters for internodes:
  #Exponential decay -> rate
  #i_beta
  #Maximum internode length
  #i_max
 make_internodes = function(N, i_beta, i_Gsd, i_noise_pct, i_plateau){
   internode.seq=i_max*exp(-seq(1:N)/i_beta)+ #negative exponential as a base
     mean(i_max*exp(-seq(1:N)/i_beta))*rnorm(N,0,sd=i_Gsd)*i_noise_pct/100 + #Gaussian noise proportional to the mean internode length
     i_plateau #Plateau value -> average minimal value of internodes towards the end of the sequence
 }
  internodes=make_internodes(N, i_beta, i_Gsd, i_noise_pct, i_plateau)
  
  #Make sure noise does not introduce negative internode length
  # When the internode is null, re-sample a value from the noised negative exponential
  null.i=which(internodes<0)
  #print(internodes)
  #print(null.i)
  if (length(null.i)>0){
    for (i in 1:length(null.i)){
      while(internodes[null.i[i]] < 0 ){
        new.i=make_internodes(N, i_beta, i_Gsd, i_noise_pct, i_plateau)
        internodes[null.i[i]]=new.i[null.i[i]]
      }
    }
  }
  #print(which(internodes<0))
  
  #Introduce random null internodes: 
  #Likelihood of null internodes: each internode length is modeled as (independant) Bernouilli variable
  null_internode_frequency=1/100 #(from real data: 3/396 intervals)
      # multiply by a random variable that takes 1 or 0 value following Poisson law
      # Use the rbinom function with sample size=1 (complementarity to inverse success/failure)
  internodes=internodes*(1-rbinom(N, 1, null_internode_frequency))
  internodes=round(internodes, digits=0)

  ref.seq=cbind.data.frame(intervals=seq(1:N),
                           angles,
                           internodes)
  
  if (natural.permutation){
    if (verbose){print("Natural permutations can be added to the divergence angle sequence")}
    #MAIN PARAMATER FOR NATURAL PERMUTATION: 'permutation.frequency' is given as input
    #default: 0.1
    
    ## Natural Permutation; Step1: Generate organ permutations
    organ.idx=seq(N+3)-1 #N+1 organs + an organ before (N°0) + an organ after (censored permutations)
    #Natural permutation is a random Bernoulli variable on each organ: =1 means that the organ is permuted with the next organ
    permut.events=rbinom(N+2, 1, permutation.frequency) #idx of perumtation can range from organ idx '0' to 'N+1', so N+2 values are possible
    
    # For debug:
    # if (verbose){
    #   print("following indexes will be permuted")
    #   print(permut.events)}
    
    #The following function allows to separate consecutive versus isolated permutations:
    permut.analysis=get_consecutive_idx(permut.events, value = 1, return.isolated = TRUE) #is a list
    
    #Case of simple isolated permutations
    simple.permut=permut.analysis[[2]]
    if (length(simple.permut>0)){
      if (verbose){print("there are isolated permutations") 
        # For debug
        # print("starting at the following indexes")
        # print(simple.permut)
        }
      organ.idx[simple.permut]=organ.idx[simple.permut]+1
      organ.idx[simple.permut+1]=organ.idx[simple.permut+1]-1
    } else {
      if (verbose){print("No isolated permutations have been drawn in the sequence")}
    }
    
    #Case of consecutive permutations
    consecutive.permut=permut.analysis[[1]]
    if (nrow(consecutive.permut)>0){
      if (verbose){print("consecutive permutations have been drawn in the sequence")
        #For debug:
        # print("They start at the following indexes")
        # print(consecutive.permut)
      }
      #print(organ.idx)
      for (r in nrow(consecutive.permut)){
        new.seq=seq(organ.idx[consecutive.permut$start.idx[r]], 
                    (organ.idx[consecutive.permut$end.idx[r]]+1))
        #print(new.seq)
        while (sum(new.seq==seq(organ.idx[consecutive.permut$start.idx[r]], 
                           (organ.idx[consecutive.permut$end.idx[r]]+1)))==length(new.seq)){
          new.seq=sample(new.seq) #this loop ensures that the shuffling is not by randomness the original order
        }
        organ.idx[seq(consecutive.permut$start.idx[r],(consecutive.permut$end.idx[r]+1))]=new.seq
      }
    } else {
      if (verbose){print("No consecutive permutations have been drawn in the sequence")}
    }
    
    ## Natural Permutation; Step2: Compute new divergence angles
    # Compute new divergence angles based on the previously generated organ permutations
    #print(organ.idx)
    #print(ref.seq$angles)
    before.angles=ref.seq$angles #store the values of the angles before permutation
    successive.gaps=organ.idx[(3:(length(organ.idx)-1))]-organ.idx[(2:(length(organ.idx)-2))]
    organ.idx=organ.idx[c(-1, -length(organ.idx))] #resize the organ idx to the window of the sequence
    # For debug only:
    # if (verbose){
    #   print("The new organ order is")
    #   print(organ.idx)
    #   print("New intervals are computed according to the following gaps")
    #   print(successive.gaps)
    # }
    
    if (length(successive.gaps) != N){stop("error: the final number of intervals has been modified by natural permutations. Consider debugging.")}
    if(verbose){print("Computing new divergence angles after natural permutations")}
    
    for (i in 1:length(successive.gaps)){
      #print("gap:")
      #print(successive.gaps[i])
      if (successive.gaps[i] != 1){
        g=successive.gaps[i]
        Oi=organ.idx[i] #original idx of the first organ of the current interval
        Oii=Oi+g #original idx of the last organ of the current interval
        #Note: [Oi-Oii] is [0, N+2]
        if ( Oi == 0 | Oi == N+2 | Oii == 0 | Oii == N+2 ){
          #All possible cases of censored permutations: a divergence angle is missing to compute the new divergence angle
          missing.angle=round(rnorm(1, mean=alpha, sd=a_sd), digits = 0) %% 360 }
        # For debug only:
        # print(c("interval, gap, first organ, second organ:"))
        # print(c(i, g, Oi, Oii))
        # print(before.angles[i])
        
        if (g>0){ 
          if (Oi == 0 & Oii== 27 ){
            ref.seq$angles[i] = (missing.angle + 
                                   sum(before.angles) +
                                   round(rnorm(1, mean=alpha, sd=a_sd), digits = 0) %% 360) %% 360
          } 
          else if (Oi == 0){#Censored permutations at the beginning: interval starting at organ '0'
            #Note: gap is necessarily positive, '0' is the minimum possible organ index
            ref.seq$angles[i] = (missing.angle + sum(before.angles[1:(g-1)]) ) %% 360
          } else if (Oii ==  N+2 ) {#Censored permutations at the end: interval ending at organ 'N+2'
            #Note: gap is necessarily positive, 'N+2' is the minimum possible organ index
            ref.seq$angles[i] = (sum(before.angles[Oi:length(before.angles)]) + missing.angle) %% 360
          } else { #all other cases
            ref.seq$angles[i] = sum(before.angles[Oi:(Oi+g-1)]) %% 360 }
          }
        else {#g<0
          if (Oi == 27 & Oii== 0 ){
            ref.seq$angles[i] = (- missing.angle - 
                                   sum(before.angles) -
                                   round(rnorm(1, mean=alpha, sd=a_sd), digits = 0) %% 360) %% 360 } 
          else if (Oii == 0){#Censored permutations at the beginning: interval ending at organ '0'
            ref.seq$angles[i] = (- missing.angle - sum(before.angles[1:Oi]) ) %% 360 }
          else if (Oi == N+2 ){#Censored permutations at the end: interval starting at organ 'N+2'
            ref.seq$angles[i] = (- ifelse(is.na(sum(before.angles[Oii:(length(before.angles))])), 0, sum(before.angles[Oii:(length(before.angles))]) ) 
                                 - missing.angle ) %% 360} 
          else {
            ref.seq$angles[i] = -sum(before.angles[(Oi+g):(Oi-1)]) %%360 }
        }
        #print(ref.seq$angles[i])
        }
    }
  }
  return(ref.seq)
}

get_consecutive_idx=function(myvector, value=value, return.isolated=TRUE){
  #DESCRIPTION: identify segments made by consecutive vector elements taking the same value given by the input [value]
  #if [input] 'return.isolated=TRUE, a second dataframe with the isolated values is returned.
  #the two dataframe are return as a list
  
  idx=which(myvector==value)
  #print(idx)
  consecutive.segments=data.frame(segment.no=integer(),
                                 start.idx=integer(),
                                 end.idx=integer(),
                                 length=integer())
  if (return.isolated){isolated.idx=c()}
  
  if (length(idx)>0){
    idx.shift=idx[-length(idx)]+1  
    breaks=c(which(idx[-1] != idx.shift), length(idx))
    #print(breaks)
    pointer=0
    segment.no=1
    for (b in breaks){
      l=b-pointer
      if (l>1){
        #print(idx[pointer+1])
        consecutive.segments=rbind.data.frame(consecutive.segments,
                                             data.frame(segment.no=segment.no, 
                                                        start.idx=idx[pointer+1], 
                                                        end.idx=idx[b],
                                                        length=l))
      segment.no=segment.no+1
      } else if (return.isolated){
        isolated.idx=c(isolated.idx, idx[b])
      }
      pointer=b
    }
  } 
  if (return.isolated){return(list(consecutive.segments, isolated.idx))}
  else {return(consecutive.segments)}
}

##################################
## Generate Measured sequences   ##
###################################
#The measure can be manual (physical device) or derived from computer-vision algorithms in a phenotyping pipeline
#Here the noise is technical and comes from the measurement procedure: it is inversely proportional to the accuracy of the measure
#In other words, while biological variations are real, technical noise is artefactual.

make_measure=function(in.seq, anoise_sd, inoise_sd,
                      noise.scale=c("mean", "sd", "absolute"),
                      anoise.mean=0, inoise.mean=0, 
                      verbose=FALSE){
  #DESCRIPTION: it simulates measured sequences from an input reference sequence by adding an independant random gaussian noise to each values of the input sequence
  #input [in.seq]: the reference sequence = dataframe with three fields: interval ($interval), angle values ($angles) and internode length ($internodes)
  #input [anoise_sd]: sd of the gaussian noise applied to angles (see noise.scale and examples to see the meaning of the exact value)
  #input [inoise_sd]: sd of the gaussian noise applied to internodes
  #input [noise.scale]: indicates what is the scale of the angle/internode noise sd provided by the user.
  #                     choice among "mean", "sd", "absolute"
  #input [anoise.mean]: the mean value of the gaussian noise for angles (default is 0: unbiased centered noise)
  #input [inoise.mean]: the mean value of the gaussian noise (default is 0: unbiased centered noise)
  #output [meas.seq]: the modified reference sequence with measurement noise added to both angles and internodes
  
  #examples:
  # if scale = "mean": anoise_sd (or inoise_sd) ranges from [0 +Inf[ and is the ratio compared to the mean of the input values
  # if scale = "sd": anoise_sd (or inoise_sd) ranges from [0 +Inf[ and is the ratio compared to the sd of the input values
  # if scale = "absolute": anoise_sd (or inoise_sd) ranges from [0 +Inf[ is the absolute value of sd (in degree for angles and in mm for internodes)                                                           
  
  #Checking inputs:
  noise.scale=match.arg(noise.scale)
  
  #Function body:
  #Notes: Technical noise has been measured for manual device: 
  #it looks Gaussian, mean error=0 sd ~ 10 / 0.7 for angles/internodes
  if (verbose){print(paste("the sd of the gaussian noise applied to input values will be scaled to", noise.scale))}
  if (noise.scale == "mean"){
    scaled.anoise_sd=mean(in.seq$angles)*anoise_sd
    angle_noise=round(rnorm(nrow(in.seq), mean=anoise.mean, sd=scaled.anoise_sd), 
                      digits = 0)
    scaled.inoise_sd=mean(in.seq$internodes)*anoise_sd
    internode_noise=round(rnorm(nrow(in.seq),mean=inoise.mean, sd=scaled.inoise_sd), 
                          digits=0)
  } else if (noise.scale == "sd"){
    scaled.anoise_sd=sd(in.seq$angles)*anoise_sd
    angle_noise=round(rnorm(nrow(in.seq), mean=anoise.mean, sd=scaled.anoise_sd), 
                      digits = 0)
    scaled.inoise_sd=sd(in.seq$internodes)*anoise_sd
    internode_noise=round(rnorm(nrow(in.seq),mean=inoise.mean, sd=scaled.inoise_sd), 
                          digits=0)
  } else {#absolute values given for noise_sd
    angle_noise=round(rnorm(nrow(in.seq), mean=anoise.mean, sd=anoise_sd), 
                      digits = 0)
    internode_noise=round(rnorm(nrow(in.seq),mean=inoise.mean, sd=inoise_sd), 
                          digits=0)
  }
  
  meas.seq=in.seq
  meas.seq$angles=(in.seq$angles+angle_noise) %% 360
  meas.seq$internodes=in.seq$internodes+internode_noise
  
  #Correct internodes sequence to avoid negative values:
  #meas.seq$internodes[meas.seq$internodes<0]=0
  #Make sure noise does not introduce negative internode length
  # When the internode is null, re-sample a value from the noised negative exponential
  null.i=which(meas.seq$internodes<0)
  if (length(null.i)>0){
    # print("correction negative values")
    # print(meas.seq$internodes)
    # print(null.i)
    if (noise.scale == "absolute"){ inoise_sd_2=inoise_sd} else {inoise_sd_2=scaled.inoise_sd} #use the right scale for the noise
    for (i in 1:length(null.i)){
      while(meas.seq$internodes[null.i[i]] < 0 ){
        new.inoise=round(rnorm(nrow(in.seq),mean=inoise.mean, sd=inoise_sd_2), 
                         digits=0)
        meas.seq$internodes[null.i[i]]=in.seq$internodes[null.i[i]]+new.inoise[null.i[i]]
      }
    }
  }
  #print(which(internodes<0))
  
  return(meas.seq)
}

#######################################################
## VALUE SEQUENCE OPERATIONS (SPLIT/MERGE EVENTS)   ###
#######################################################
#Following functions require a input dataframe [in.seq], which is a dataframe with specific format requirements, as follows:
# the dataframe contains three fields: 
#             intervals ($intervals), 
#             angle values ($angles) 
#             internode length ($internodes)

#1. MAKE FREE ENDS
value_chops=function(in.seq, base, top){
  #DESCRIPTION: this function removes values at the base/beginning and top/end of a sequence
  #input [in.seq]: the reference sequence=dataframe with three fields: intervals ($intervals), angle values ($angles) and internode length ($internodes)
  #input [base]: the number of intervals to remove at the beginning/base of the sequences
  #note: this is equivalent to the number of organs to remove at the beginning of the sequence
  #input [top]: the number of intervals to remove at the top/end of the sequences
  #note: this is equivalent to the number of organs to remove at the end of the sequence
  #output: add a column tracking the interval matching to the reference/input sequence
  
  #Input checks
  if (top+base>nrow(in.seq)){
    stop("deleting more rows than existing")
  }
  if (top+base == nrow(in.seq)){
    return(NULL)
  }
  
  #Chop free ends at start/end of the sequence
  out.seq=in.seq[(base+1):(nrow(in.seq)-top),]
  
  #Compute new intervals
  out.seq$intervals=seq(1:nrow(out.seq))

  return(out.seq)
}

value_tails=function(in.seq, base, top){
  #DESCRIPTION: this function creates values at the base/beginning and top/end of a sequence
  #input [in.seq]: the reference sequence=dataframe with three fields: intervals ($intervals), angle values ($angles) and internode length ($internodes)
  #input [base]: the number of intervals to add at the beginning/base of the sequences
  #note: this is equivalent to the number of organs to add at the beginning of the sequence
  #input [top]: the number of intervals to add at the top/end of the sequences
  #note: this is equivalent to the number of organs to add at the end of the sequence
  #output: add a column tracking the interval matching to the reference/input sequence
  
  amplif=1 #parameter to amplify possible values taken by the new values, range [0 ; +inf]
  local=5 #Parameter that set the Nber of organs taken into account to compute the values that will be added (range: [1 ; nrow(in.seq)])
  base_tail_angles=round(rnorm(base, 
                               mean=mean(in.seq$angles[1:local]),
                               sd=amplif*sd(in.seq$angles[1:local])),
                         digits = 0) %% 360
  base_tail_internodes=round(rnorm(base,
                                   mean=mean(in.seq$internodes[1:local]),
                                   sd=amplif*sd(in.seq$internodes[1:local])),
                              digits = 0)
  #correct possible negative values for internode length:
  base_tail_internodes[base_tail_internodes<0]=0
  
  top_tail_angles=round(rnorm(top, 
                              mean=mean(in.seq$angles[(nrow(in.seq)-local):nrow(in.seq)]),
                              sd=amplif*sd(in.seq$angles[(nrow(in.seq)-local):nrow(in.seq)])),
                         digits = 0) %% 360
  top_tail_internodes=round(rnorm(top, 
                                  mean=mean(in.seq$internodes[(nrow(in.seq)-local):nrow(in.seq)]),
                                  sd=amplif*sd(in.seq$internodes[(nrow(in.seq)-local):nrow(in.seq)])),
                              digits = 0)
  #correct possible negative values for internode length:
  top_tail_internodes[top_tail_internodes<0]=0
  
  #add values in the dataframe:
  out.seq=data.frame(intervals=seq(1,(base+nrow(in.seq)+top)),
                     angles=c(base_tail_angles, in.seq$angles, top_tail_angles),
                     internodes=c(base_tail_internodes, in.seq$internodes, top_tail_internodes))
  return(out.seq)
}

#2. MERGE EVENTS (create two values that should be merged to the initial reference value)
create_merge_angle=function(in.seq){
  #DESCRIPTION: this function creates the first additional angle value of typical MERGE cases. 
  #It has no constraints but the second value is fully determined by the initial value and this new created one.
  #Note for improvement: more sophisticated values could take into account local values only around the index impacted by the MERGE (hence this index should be given as input)
  
  amplif=0 #parameter to amplify possible values taken by the new angle, range [-Inf (increase); 0-1 (decrease)]
  round(runif(1, 
              min=(1-amplif)*min(in.seq$angles), #cannot be <0
              max=min(360,(1+amplif)*max(in.seq$angles))), #cannot be more than 360°
        digits = 0)
}
create_merge_internode=function(in.seq,idx){
  #DESCRIPTION: the function that creates the first additional internode value like in MERGE cases (two values are added by a merge)
  #Note: more sophisticated values could take into account local values only around the index impacted by the MERGE (hence this index should be given as input)
  if (idx==1 || idx==(nrow(in.seq)+1)){
    #No constraints on the new value
    round(runif(1, min=0, max=in.seq$internodes),
          digits = 0)
  } else {
    #The new value cannot be larger than the existing interval !
    #take a random value in the existing interval
    round(runif(1, min=0, max=in.seq$internodes[idx]),
          digits = 0)
  }
  
}
value_merge=function(in.seq, i, verbose=FALSE){
  #DESCRIPTION: this function creates TWO values that should be MERGED to a single value when aligned the reference sequence
  #a MERGE corresponds to an ADDITION of a value in the test sequence compared to the reference sequence
  #It impacts one value of the reference sequence and add a new one.
  #Biologically, it corresponds to an ADDITION of an organ (over-segmentation) 
  #input [in.seq]: the reference sequence=dataframe with three fields: intervals ($intervals), angle values ($angles) and internode length ($internodes)
  #input [i]: the (new) idx of the FALSE ORGAN which is ADDED (ie the idx taken by the false organ in the new sequence). All true idx downstream will be shifted +1.
  #Note: for MERGE, [i] can range from 1 (1: organ added before the first ref organ) to the end of seq +2 (organ added after the LAST ref organ)
  #output [out.seq]: the test sequence with ONE MORE value (the (i-1)th and i.th value must be merged to the (i-1)th value of the aligned ref sequence). 
  
  if (i<1 || i>(nrow(in.seq)+2)){
    stop("i out of range")
  }
  
  if (i==1){
    # use the function value_tails
    out.seq=value_tails(in.seq, 1, 0)
    if (verbose){ cat("This simply adds values at the beginning of sequences (not a proper merge case). The function 'value_tails' was used instead.\n") }
  } 
  else if (i<=nrow(in.seq)+1){
    #Generate a new dataframe, from start up to the involved interval i-1 
    #and introduce an additional interval created by the MERGE by duplicating the (i-1)th row
    out.seq=rbind(in.seq[1:i-1,],in.seq[i-1,])
    #erase previous values for the two last intervals corresponding to idx (i-1) and i of new output seq:
    out.seq[i-1,c(2,3)]=c(NA, NA) #only modify angles + internodes
    out.seq[i,c(1,2,3)]=c(i,NA, NA) #modify index + angles + internodes
    #if needed, retrieve the end of the reference sequence:
    if (i<=nrow(in.seq)){#if i=nrow+1, all the sequence is already loaded in out.seq
      end.seq=in.seq[(i):nrow(in.seq),]
      #modify interval numbering after the MERGE (add +1):
      end.seq$intervals=end.seq$intervals+1
      out.seq=rbind(out.seq, end.seq)
      }
    #Create MERGE values
    #Angle values
      #The two new values should sum-up to previous value
      out.seq[i-1,]$angles=create_merge_angle(in.seq)
      #name this new angle new_alpha
      new_alpha=out.seq[i-1,]$angles
      if (new_alpha < in.seq[i-1,]$angles){
        out.seq[i,]$angles=in.seq[i-1,]$angles-new_alpha} 
      else {
        out.seq[i,]$angles=360-new_alpha+in.seq[i-1,]$angles
      }
    #Internode values
      #The two new values should sum-up to previous value
      out.seq[i-1,]$internodes=create_merge_internode(in.seq,i-1)  
      out.seq[i,]$internodes=in.seq$internodes[i-1]-out.seq[i-1,]$internodes
  } 
  else { #meaning if i==nrow(in.seq)+2
    #use the function value_tails
    out.seq=value_tails(in.seq, 0, 1)
    if (verbose){ cat("This simply adds values at the end of sequences (not a proper merge case). The function 'value_tails' was used instead.\n")}
    
  }
  
  return(out.seq)
}

#3. SPLIT EVENTS (create a new value that should be split into two initial reference values)
value_split=function(in.seq, i){
  #DESCRIPTION: this function creates a NEW value that should be SPLIT into two reference values when aligned the reference sequence
  #a SPLIT corresponds to a LOSS of a value in the test sequence compared to the reference sequence
  #it impacts two values of the reference sequence
  #Biologically, it corresponds to an ABSENCE of an organ (under-segmentation)
  #input [in.seq]: the reference sequence=dataframe with three fields: intervals ($intervals), angle values ($angles) and internode length ($internodes)
  #input [i]: idx of the ORGAN which is missed (impacts the (i-1)th and the i.th interval)
  #Note: for SPLIT, [i] can range from 1 to the end of seq +1 (ie. nrow(in.seq)+1 )
  #output [out.seq]: the test sequence with one less value, the (i-1)th of the test seq must be split to the (i-1)th and (i)th position in aligned ref seq
  
  if (i<1 || i>(nrow(in.seq))+1){
    stop("i out of range")
  }
  
  if (i==1){
    #case where the first organ is missed, so only the first values of the reference seq is affected
    out.seq=in.seq[2:nrow(in.seq),]
    #Offset of -1 to all previous intervals
    out.seq$intervals=out.seq$intervals-1
  } else if (i<=nrow(in.seq)-1){
    #remove the i-th interval:
    out.seq=rbind(in.seq[1:(i-1),], in.seq[(i+1):nrow(in.seq),])
    #rename intervals
    out.seq$intervals=1:(nrow(in.seq)-1)
    #change value of the (i-1)th: the sum of previous (i-1)th and i-th:
    out.seq$angles[i-1]=(in.seq$angles[i-1]+in.seq$angles[i]) %% 360
    out.seq$internodes[i-1]=(in.seq$internodes[i-1]+in.seq$internodes[i])
  } else if (i==nrow(in.seq)) {
    #remove the i-th and last interval:
    out.seq=in.seq[1:(i-1),]
    #change value of the (i-1)th: the sum of previous (i-1)th and i-th:
    out.seq$angles[i-1]=(in.seq$angles[i-1]+in.seq$angles[i]) %% 360
    out.seq$internodes[i-1]=(in.seq$internodes[i-1]+in.seq$internodes[i])
  } else {
    #Necessary, i==nrow(in.seq)+1
    #case where the last organ has been missed, so only the last values of the reference seq. are affected
    out.seq=in.seq[1:(nrow(in.seq)-1),]
  }
  
  return(out.seq)
}

#4. PERMUTATIONS
value_simple_permut=function(in.seq, idx){
  #DESCRIPTION: this function computes the new angle values affected by a single 2-permutation
  #Note: here, permutations are only allowed between close organs (short internodes). 
  #     Internode values are unchanged, the longitudinal position of organs along the stem are just swapped
  #input [idx]: integer, index of the central interval containing the two permuted organs
  #output [out.seq]: in.seq (angles) values modified by the permutation
  
  #check input arguments
  if (idx<1 || idx>nrow(in.seq)){
    stop(paste0("input index out of range: it should be in [1-", nrow(in.seq),"]"))
  }
  
  #Main body
  out.seq=in.seq
  #Modify the angle values:
  if (idx==1){
    out.seq$angles[idx]=(360-in.seq$angles[idx]) %% 360
    out.seq$angles[idx+1]=(in.seq$angles[idx]+in.seq$angles[idx+1]) %% 360
  } else if (idx==nrow(in.seq)){
    out.seq$angles[idx-1]=(in.seq$angles[idx-1]+in.seq$angles[idx]) %% 360
    out.seq$angles[idx]=(360-in.seq$angles[idx]) %% 360
  } else {
    out.seq$angles[idx-1]=(in.seq$angles[idx-1]+in.seq$angles[idx]) %% 360
    out.seq$angles[idx]=(360-in.seq$angles[idx]) %% 360
    out.seq$angles[idx+1]=(in.seq$angles[idx]+in.seq$angles[idx+1]) %% 360
  }
  
  #Return
  return(out.seq)
}

#############################################
## Intervals/Organs SEQUENCE OPERATIONS   ###
#############################################
#The following functions use a list as input, which should respect a particular format. Here is the specification:
#[align.list]: a list made of two dataframes as element:
  #Dataframe n°1 [Ialign]: contains the alignment of intervals between reference and modified sequence
  #                         Formatted with three invariant fields: $reference, $modified, $dtw
  #                         $reference and $modified contain either the index/rank of the interval or NA
  #                         $dtw contains a code indicating which operation(s) are needed to get modified interval values back to reference values. 
  #                         $dtw: the code is a combination of the following basic operations: ~=identity / M=merge / S=split / C=chop / T=tail
  #Dataframe n°2 [Oalign]: contains the alignment of organs between reference and modified sequence
  #                         Formatted with three invariant fields: $reference, $modified, $segmentation
  #                         $reference and $modified contain either the index/rank of the interval or NA
  #                         $segmentation contains a factor indicating what happen to the reference organ in the modified sequence
  #                         $segmentation factor levels are given by the function 'define_seg_levels'

#accessory functions
define_seg_levels=function(vec){
  ##DESCRIPTION: define the authorized levels expected for organ segmentation
  #input [vector]: typically align.list$Oalign$segmentation (cf 'align.list' )
  
  #input checks
  if ( dim.data.frame(vec)[1] != 0 ){stop("input is not 1D vector")}
  
  #definition of factor levels
  organ.segmentation.level.set=c("~", "over", "under", 
                                 "chop", "tail", 
                                 "chop/tail","chop/over","under/over",
                                 "perm")
  
  #function's body
  vector=factor(vec, levels= organ.segmentation.level.set)
  
  return(vec)
}

dtw_code_compute=function(vec, idx, code){
  #DESCRIPTION: vec is list.align$Ialign$dtw in the functions below. 
  #This function creates a new dtw code by concatenating the previous record with the new operation (given by [code] input value)
  #input [idx]: indexes of vec to be changed. idx can be a vector of several indexes
  #input [code]: code to add in the rows given by idx
  #output: vec with updated code
  
  #Internal function: allows to erase "~" from code concatenation
  clean_code=function(car){ifelse(car=="~", "", car)}
  
  #Main body
  #vec data must be factors
  vec=as.factor(vec)
  #Order idx
  ordered_Idx=idx[order(idx)]
  #Store previous code(s) at each row
  code_record=vec[ordered_Idx]
  #Compute the new code: concatenation of input [code] with previous content:
  if (length(code_record) == 1 && is.na(code_record)){ new_code = code } #typical case with only one row to change
  else { new_code = paste0(clean_code(as.character(code_record)), code) }
  #Manage levels
  if (sum(new_code %in% levels(vec)) == length(idx)){
    #New code(s) was(were) already present and hence is(were) already an existing level.
    #Just replace the content of the dtw cells in the df:
    vec[ordered_Idx]=new_code } 
  else {
    #the level set needs to be updated to include new code(s)
    cat_level=c(levels(vec), new_code)
    #Note: duplicated levels must be eliminated while defining new levels with 'factor' function
    vec=factor(vec, levels=cat_level[!duplicated(cat_level)])
    #Cell content can now be replaced
    vec[ordered_Idx]=new_code }
  return(vec)
}

make_align_list=function(nb, type=c("interval", "organ")){
  #DESCRIPTION: generate rapidly a list called 'align.list' with the desired length. $ref and $modified are all aligned
  #input [nb]: nb of intervals or organs
  #input [type]: if nb refers to interval or organs
  #ouput: list with the expected correct format for subsequent operations
  
  #check input
  type=match.arg(type)
  
  #main body
  if (type == "interval"){N=nb+1} else {N=nb}
  O=data.frame(reference=c(1:N),
                modified=c(1:N),
                segmentation=c(rep("~",N)))
  I=data.frame(reference=c(1:(N-1)),
                modified=c(1:(N-1)),
                dtw=c(rep("~",(N-1))))
  align.list=list(Ialign=I, Oalign=O)
  
  #output
  return(align.list)
}

check_align_list=function(align.list, verbose=FALSE){
  ## DESCRIPTION: function checking that align.list argument is well formatted to be used by other functions of the program
  ##[output]: the list with standard names for list element, ordered as ($Ialign, $Oalign) and standard dataframe columns
  
  #1. check the format = list
  if (!is.list(align.list)){stop("input must be a list")}
  #2. check list length = 2 elements
  if (length(align.list) != 2){
    if (length(align.list) ==3){cat("Debug advice: Check whether the '$value' dataframe was not present in the input list.\n")}
    stop(paste("input data has", length(align.list),"fields, while alignment list should typically contain only 2 fields: interval and organ alignments."))}

  #3. check that each list element is a dataframe
  if (sum(sapply(align.list, FUN=is.data.frame)) != length(align.list)){
    print("Error: input list must contain only 'dataframe' elements")
    for (i in 1:length(align.list)){
      if (!is.data.frame(align.list[[i]])){
        print(paste0("tracking error: element list '", align.list[i], "' (n° ", i, ") is not formated as dataframe"))
      }
    }
    stop("reformat input before using the function")
  }
  
  #4. Identify which elements correspond to Interval/Organ alignment dataframes
  #if no indices: default is Interval first, organ is second.
  #First guess on colnames of the dataframes, then on the column content
  if ("dtw" %in% c(colnames(align.list[[1]]), colnames(align.list[[2]]))){
    if ("dtw" %in% colnames(align.list[[2]])){ 
      align.list=rev(align.list) }}
  else {
    if (verbose){cat("Absence of standard names for alignment data. Guessing assignment from content...\n")}
    #Detect dtw code content in the third column:
    if ((length(grep("[S,M,T,C]", t(align.list[[2]][3]))) > 0 && length(grep("[S,M,T,C]", t(align.list[[1]][3]))) == 0 )){
      align.list=rev(align.list)}
    #Force the dataframe colnames
    if (verbose){cat("Correcting column names of alignments. Double check that the content is appropriate and well formated.\n")}
    colnames(align.list[[1]])=c("reference", "modified" , "dtw")
    colnames(align.list[[2]])=c("reference", "modified", "segmentation")
    }
  names(align.list)=c("Ialign", "Oalign")
  
  #Note for future development/improvement: also force names of each fields of the df ?
  return(align.list)
}

#Functions to manipulate intervals
seq_insert=function(align.list, idx, verbose=FALSE){
      ##DESCRIPTION: re-align intervals and element ranks of two sequences after inserting an element (organ) inside the sequence.
      ##The insertion is performed on the so-called "modified" input sequence and new sequences are aligned against the so-called "reference" sequence
      ##input [align.list]: a list made of two dataframes (Ialign and Olign) as element (see specifications)
      ##input [idx]: the idx of the organ to be added in the input modified sequence-> it is the idx the added element (organ) will have in the new/modified sequence
      ## /!\ Warning: do not use the organ-idx of the reference sequence. All modifications must be designed from the input modified sequence.
      ##output: updated list (Ialign/Oalign) modified with a new interval/organ and corresponding operations 
  
  #####
  #input checks
  #####
  align.list=check_align_list(align.list)
  
  #Get the number of organs for reference and currently modified sequences
  No=max(align.list$Oalign$modified, na.rm = TRUE)
  Noref=max(align.list$Oalign$reference, na.rm = TRUE)
  if (verbose){
    cat(" Reference sequence contains", Noref,"organs.\n", 
        "Currently modified sequence contains", No, "organs.\n")  
  }
  
  #####
  #Main body
  #####
  #idx is in range [2, n+1]/where n is the number of intervals; or [2, No]/No Number of organs defined above
  if (idx<1 || idx>No+1){
    stop(paste0("index out of range: it should range from 2 to ", No," (included)")) }
  else if (idx==1 || idx==No+1){
    if (verbose){
      cat("Appending organs at a sequence tail -> use 'seq_append' function instead. \n") 
    }
    align.list=seq_append(align.list, idx, verbose=verbose) }
  else {
    #####
    #by default, set encoding parameters to FALSE
    gap=FALSE #gets activated -TRUE- only if there are true organs chopped before/after the current start/end of the modified sequence
    skip=FALSE #gets activated -TRUE- only if there are false organs tailing (at start or end) the current the modified sequence
    
    #get the index of 'idx' in the Oalign$modified and Ialign$modified table:
    idx.=which(align.list$Oalign$modified==idx)
    NrO=nrow(align.list$Oalign) #nber of rows in the table of organs 'align.list$Oalign'
    idxi.=max(which(align.list$Ialign$modified==(idx-1))) #index of interval, '.' indicates that the considered object is the table containing the interval indexes: idxi. is the rank of that table cell
    NrI=nrow(align.list$Ialign) #Nber of rows in the intervals dataframe
    
    #First, check whether a true reference organ was removed before in this interval
    #In that case, the new insertion in $modified will re-align with that organ
    if (is.na(align.list$Oalign$modified[(idx.-1)])) {#in the alignment with $ref organs, there is a gap before 
      gap=TRUE
      ##Organs
      #change this NA by the organ index 'idx':
      align.list$Oalign$modified[(idx.-1)]=idx
      #shift indexes following this insertion for $modified organs:
      align.list$Oalign$modified[idx.:NrO]=align.list$Oalign$modified[idx.:NrO]+1
      #define segmentation event for the inserted organ:
      if (align.list$Oalign$segmentation[(idx.-1)] != "under"){
        warning(c("expected segmentation 'under' was not found in the organ table for cell row n° ", idx., "."))}
      align.list$Oalign$segmentation=define_seg_levels(align.list$Oalign$segmentation)
      align.list$Oalign$segmentation[(idx.-1)]="under/over"
      
      ##Intervals
      #No new interval will be created: the existing gap in the alignment will be filled by the new interval
      #change the $modified indexes from there up to the end:
      align.list$Ialign$modified[idxi.:NrI]=align.list$Ialign$modified[idxi.:NrI]+1
      #update dtw code:
      align.list$Ialign$dtw=dtw_code_compute(align.list$Ialign$dtw, idx=c((idxi.-1), idxi.), "M" )
    }
    else {#No existing gap before: a new row will be added in the table to keep organ alignment
      if (idx. < NrO){
        if (align.list$Oalign$segmentation[(idx.-1)]=="tail"){#There is already a tail BEFORE
          align.list=seq_append(align.list, 1, verbose=verbose) #so extend the existing starting tail
          skip=TRUE }
        else {#duplicate a row
          align.list$Oalign=rbind(align.list$Oalign[1:idx.,],align.list$Oalign[idx.,], align.list$Oalign[(idx.+1):NrO,]) }
      } 
      else {#meaning idx. == NrO
        if (align.list$Oalign$segmentation[(idx.)]=="tail"){#the last organ is already a tail
          align.list=seq_append(align.list, NrO+1, verbose=verbose) #so extend the existing ending tail
          skip=TRUE }
        else {
          #just append a row at the end by duplicating the last row
          align.list$Oalign=rbind(align.list$Oalign[,] , align.list$Oalign[idx.,]) }
      }
      if (!skip){ #no existing tail has been detected
        ## Organs
        #Re-index organs in the modified sequence and in the newly modified sequence
        #By definition in this function, the inserted organ gets the index given by the input 'idx' value
        #In the reference sequence:
        align.list$Oalign$reference[idx.]=NA
        #in case all the levels have not been defined for $segmentation
        align.list$Oalign$segmentation=define_seg_levels(align.list$Oalign$segmentation)
        #then flag this inserted organ as a "over" ("over-segmentation")
        align.list$Oalign$segmentation[idx.]="over"
        #In the table, shift all indexes of modified sequences after the inserted organ (ie., up to NrO+1)
        align.list$Oalign$modified[(idx.+1):(NrO+1)]=align.list$Oalign$modified[(idx.+1):(NrO+1)]+1
        
        ## Intervals
        #so insert an new interval in the modified sequence by duplicating the idx-th row
        #debug: print(paste("NrI =",NrI))
        if (idxi. < NrI){
          align.list$Ialign=rbind(align.list$Ialign[1:idxi.,],align.list$Ialign[idxi.,], align.list$Ialign[(idxi.+1):NrI,]) } 
        else {#meaning idx == No: just append a row at the end by duplicating the last row
          align.list$Ialign=rbind(align.list$Ialign[,] , align.list$Ialign[idxi.,]) }
        ##debug: print(align.list$Ialign)
        #Annotate cells of Ialign table:
        #In the reference column, the idx will be duplicated, as expected
        #shift downstream indexes of the modified sequence:
        align.list$Ialign$modified[(idxi.+1):(NrI+1)]=align.list$Ialign$modified[(idxi.+1):(NrI+1)]+1 #Note: NA values stay NA
        #dtw:
        #create a dtw code by concatenating previous code with the new operation:
        #Note: if previous code is just "~", it will be eliminated
        align.list$Ialign$dtw=dtw_code_compute(align.list$Ialign$dtw, idx=c(idxi., (idxi.+1)), "M" )
      }
    }
    #####
    #Outputs
    #User information message
    No=max(align.list$Oalign$modified, na.rm = TRUE)
    if (verbose){ cat(" After edition, the modified sequence now contains", No,"organs.\n") }
  }
  return(align.list)
}
seq_remove=function(align.list, idx, verbose=FALSE){
  ##DESCRIPTION: re-align intervals and element ranks of two sequences after REMOVING an element (organ) inside the sequence.
  ##The removal is performed on the so-called "modified" input sequence and new sequences are aligned against the so-called "reference" sequence
  ##input [align.list]: a list made of two dataframes (Ialign and Olign) as element (see specifications)
  ## /!\ Warning: do not use the organ-idx of the reference sequence. All modifications must be designed from the input modified sequence.
  ##output: updated list (Ialign/Oalign) modified intervals/organs and corresponding operations
  
  #####
  #input checks
  #####
  align.list=check_align_list(align.list)
  
  #Get the number of organs for reference and currently modified sequences
  No=max(align.list$Oalign$modified, na.rm = TRUE)
  Noref=max(align.list$Oalign$reference, na.rm = TRUE)
  if (verbose){
    cat(" Reference sequence contains", Noref,"organs.\n", 
        "Currently modified sequence contains", No, "organs.\n")
  }
  #idx is in range [1, n+1]/where n is the number of intervals; or [2, No+1]/No Number of organs defined above
  if (idx < 1 || idx > No){
    stop(paste0("index out of range: it should range from 1 to ", No," (included)")) }
    
  #####
  #Main body
  #####
  ## Impact on 'Organ' dataframe and alignement.
  #get the index of 'idx' in the Oalign$modified table:
  idx.=which(align.list$Oalign$modified==idx)
  #replace the value of organ idx by NA
  align.list$Oalign$modified[idx.]=NA
  #if possible, shift all downstream values of idx
  NrO=nrow(align.list$Oalign) #nber of rows in the table of organs 'align.list$Oalign'
  if (idx. < NrO) {
    align.list$Oalign$modified[idx.:NrO]=align.list$Oalign$modified[idx.:NrO]-1  
  } #no else needed if idx.==NrO
  #in case all the levels have not been defined for $segmentation
  align.list$Oalign$segmentation=define_seg_levels(align.list$Oalign$segmentation)
  #then flag this inserted organ as a "under" ("under-segmentation")
  if ( idx == 1 || idx == No) {align.list$Oalign$segmentation[idx.]="chop"} 
  else { align.list$Oalign$segmentation[idx.]="under"}
  
  #####
  ## Impact on 'Interval' dataframe and alignement.
  #In modified sequence: if removed organ has index=idx, the impacted interval will be of index (idx) as well, except for idx=NrO (last organ)
  #Note: interval indexes can be duplicated 
  if (idx == No){#select the last interval
    idxi.=which(align.list$Ialign$modified==(No-1)) } 
  else {idxi.=which(align.list$Ialign$modified==(idx)) #idxi. is the rank of the table cell align.list$Ialign$modified 
  }
  NrI=nrow(align.list$Ialign)
  if ( idx == 1 ) {
    align.list$Ialign$modified[idxi.]=NA #if several indexes, all are converted to NA
    align.list$Ialign$modified[(max(idxi.)+1):NrI]=align.list$Ialign$modified[(max(idxi.)+1):NrI]-1
  } else if (idx == No) {
    align.list$Ialign$modified[idxi.]=NA #if several indexes, all are converted to NA
  }
  else {#shift idxi. and all downstream values
    align.list$Ialign$modified[min(idxi.):NrI]=align.list$Ialign$modified[min(idxi.):NrI]-1
  }
  #dtw
  #create a dtw code by concatenating previous code with the new operation:
  #Note: if previous code is just "~", it will be eliminated
  if ( idx == 1 ) {
    align.list$Ialign$dtw=dtw_code_compute(align.list$Ialign$dtw, idx=idxi., code="C" )}
  else if (idx == No){
    if (length(idxi. > 1)){
      if (!("C" %in% levels(align.list$Ialign$dtw))){
        align.list$Ialign$dtw=factor(align.list$Ialign$dtw, levels=c(levels(align.list$Ialign$dtw), "C"))
      }
      align.list$Ialign$dtw[idxi.]="C"
    }
    else {
      align.list$Ialign$dtw=dtw_code_compute(align.list$Ialign$dtw, idx=idxi., code="C" )}  
    }
  else {#"S" is the code for split
    align.list$Ialign$dtw=dtw_code_compute(align.list$Ialign$dtw, idx=c((idxi.-1), idxi.), "S" )
  }

  #####
  #Outputs
  #User information message
  No=max(align.list$Oalign$modified, na.rm = TRUE)
  if (verbose){ cat(" After edition, the modified sequence now contains", No,"organs.\n") }
  return(align.list)
}
seq_append=function(align.list, idx, verbose=FALSE){
  ##DESCRIPTION: append intervals and element (organs) two sequences before/after given sequences.
  ##Additions are performed on the so-called "modified" input sequences and new sequences are aligned against the so-called "reference" sequence
  ##input [align.list]: a list made of two dataframes (Ialign and Olign) as element (see specifications)
  ##input [idx]: the idx of the organ to be added in the input modified sequence-> it is the idx the added element (organ) will have in the new/modified sequence
  ##for this function, idx can then only be 1 or (Nb of organ + 1)
  ## /!\ Warning: do not use the organ-idx of the reference sequence. All modifications must be designed from the input modified sequence.
  ##output: updated list (Ialign/Oalign) modified with a new interval/organ and corresponding operations 
  
  #####
  #input checks
  #####
  align.list=check_align_list(align.list)
  
  #Get the number of organs for reference and currently modified sequences
  No=max(align.list$Oalign$modified, na.rm = TRUE)
  Noref=max(align.list$Oalign$reference, na.rm = TRUE)
  if (verbose){
    cat(" Reference sequence contains", Noref,"organs.\n", 
        "Currently modified sequence contains", No, "organs.\n")
  }
  
  #idx is either == 1 or ==No (Number of organs defined above)
  if (!(idx == 1 || idx == (No+1))){
    stop(paste0("incorrect index value: to append a tail, it should be equal to either 1 (start) or ", No+1," (Nber of organs +1 in modified sequence)")) }
  
  #####
  #Main body
  #####
  NrO=nrow(align.list$Oalign)
  NrI=nrow(align.list$Ialign)
  #default behavior: the new organ is appended just before the current 1st organ of the modified sequence
  if (idx==1){#user want to append organ at the beginning of the modified sequence
    idx.=which(align.list$Oalign$modified==idx) #index of 'idx' in the Oalign$modified table
    if (idx.==1){#the current first organ starts the sequence, ie. there are no chopped reference organs
      ## ORGANS
      #duplicate the row containing the first modified organ
      align.list$Oalign=rbind(align.list$Oalign[1,], align.list$Oalign[,])
      #change reference rank value:
      align.list$Oalign$reference[1]=NA
      #shift idx after addition
      align.list$Oalign$modified[2:(NrO+1)]=align.list$Oalign$modified[2:(NrO+1)]+1
      #change $segmentation:
      #add the "tail" flag (it must be checked that this level belongs to the levels of $segmentation factor)
      if ( "tail" %in% levels(align.list$Oalign$segmentation) ) {
        align.list$Oalign$segmentation[idx]="tail" } 
      else {
        align.list$Oalign$segmentation=define_seg_levels(align.list$Oalign$segmentation)
        align.list$Oalign$segmentation[idx]="tail" }
      
      ## INTERVALS
      #safety check
      if (which(align.list$Ialign$modified == idx) != 1 ){
        print(c("debug: idx is ", idx, "(table idx of the corresponding $modified organ) is ", idx., ". Aligned sequence of organs is:"))
        print(align.list$Ialign)
        print("the first $modified interval should start the table")
        stop("an error has occured: conflict between organ and interval table")
      }
      #duplicate the row containing the first modified interval
      align.list$Ialign=rbind(align.list$Ialign[1,], align.list$Ialign[,])
      #change reference rank value:
      align.list$Ialign$reference[1]=NA
      #shift idx after addition
      align.list$Ialign$modified[2:(NrI+1)]=align.list$Ialign$modified[2:(NrI+1)]+1
      #dtw code:
      #change the duplicated code to avoid errors:
      align.list$Ialign$dtw[idx]=NA
      align.list$Ialign$dtw=dtw_code_compute(align.list$Ialign$dtw, idx=1, code = "T") #code for tail is "T"
    } 
    else { #idx.>1
      #there are (at least 1) chopped reference organ(s) that start the $Oalign df
      ## ORGANS
      #so 'NA' should be aligned to this ref organ in the $modified organ
      if ( !is.na(align.list$Oalign$modified[(idx.-1)]) ){
        print(c("debug: idx is ", idx, "idx. (idx of 1st $modified organ) is ", idx., "aligned sequence of organs is:"))
        print(align.list$Oalign)
        stop(paste("error: a NA is expected in the 'modified' column at row", (idx.-1) ) ) }
      #change the value of NA
      align.list$Oalign$modified[(idx.-1)]=1
      #shift all downstream ranks
      align.list$Oalign$modified[idx.:NrO]=align.list$Oalign$modified[idx.:NrO]+1
      #change $segmentation:
      #add the "chop/tail" flag (it must be checked that this level belongs to the levels of $segmentation factor)
      if ( "chop/tail" %in% levels(align.list$Oalign$segmentation) ) {
        align.list$Oalign$segmentation[(idx.-1)]="chop/tail" } 
      else {
        align.list$Oalign$segmentation=define_seg_levels(align.list$Oalign$segmentation)
        align.list$Oalign$segmentation[(idx.-1)]="chop/tail" }
      
      ## INTERVALS
      #likely 'NA' should be aligned to the previous ref interval in the $modified interval
      if ( !is.na(align.list$Ialign$modified[(idx.-1)]) ){
        print(c("debug: idx is ", idx, "idx. (idx of 1st $modified intervak) is ", idx., "aligned sequence of interval is:"))
        print(align.list$Ialign)
        stop(paste("error: a NA is expected in the 'modified' column at row", (idx.-1) ) ) }
      #change the value of NA
      align.list$Ialign$modified[(idx.-1)]=1
      #shift all downstream ranks
      align.list$Ialign$modified[idx.:NrI]=align.list$Ialign$modified[idx.:NrI]+1
      #dtw code:
      align.list$Ialign$dtw=dtw_code_compute(align.list$Ialign$dtw, idx=(idx.-1), code="T") #code for tail is "T"
    }
  }
  else {#idx==No+1: the user want to append organ at the end of the modified sequence
    idx.=which(align.list$Oalign$modified==(idx-1)) #the rank of the table cells that contain the last $modified organ
    if (idx. == NrO){ #existing last $modified organ is precisely at the end of the table
      
      ## Organs
      #duplicate the last row
      align.list$Oalign=rbind(align.list$Oalign[,], align.list$Oalign[NrO,])
      #change rank values:
      align.list$Oalign$reference[(idx.+1)]=NA
      #shift last $modified row:
      align.list$Oalign$modified[(idx.+1)]=No+1
      #change $segmentation:
      #add the "tail" flag (it must be checked that this level belongs to the levels of $segmentation factor)
      if ( "tail" %in% levels(align.list$Oalign$segmentation) ) {
        align.list$Oalign$segmentation[(idx.+1)]="tail" } 
      else {
        align.list$Oalign$segmentation=define_seg_levels(align.list$Oalign$segmentation)
        align.list$Oalign$segmentation[(idx.+1)]="tail" }
      
      ##Intervals
      #idx is the last organ to append -> there are currently (idx-1) organs in the sequence
      #                                -> there are currently (idx-2) intervals in the sequence
      idxi.=max(which(align.list$Ialign$modified == (idx-2))) #is the table (Ialign) index of the last $modified interval
      #safety check
      if (idxi. != NrI ){
        print(paste("debug: idx is ", idx, ". idxi. (table index of the corresponding $modified interval) is ", idxi., ". Aligned sequence of organs is:"))
        print(align.list$Ialign)
        print("the last $modified interval should end the table")
        stop("an error has occured: conflict between organ and interval table")
      }
      #duplicate the row containing the last $modified interval
      align.list$Ialign=rbind(align.list$Ialign[,], align.list$Ialign[idxi.,])
      #change reference rank value:
      align.list$Ialign$reference[(idxi.+1)]=NA
      #shift idx of the addition
      align.list$Ialign$modified[(idxi.+1)]=align.list$Ialign$modified[(idxi.)]+1
      #dtw code:
      #change the duplicated code to avoid errors:
      align.list$Ialign$dtw[(idxi.+1)]=NA
      align.list$Ialign$dtw=dtw_code_compute(align.list$Ialign$dtw, idx=(idxi.+1), code = "T") #code for tail is "T"
    } 
    else { #idx.< NrO: the existing last $modified organ is not at the end of the table (there are -at least one- chopped $ref organs)
      
      ## Organs
      #so there is a stretch of chopped organs ending the modified sequence. Reference organs are aligned to 'NAs'.
      if ( sum(is.na(align.list$Oalign$modified[(idx.+1):NrO])) != (NrO - idx.) ){
        print(c("debug: idx is ", idx, "idx. (idx of last $modified organ) is ", idx., "aligned sequence of organs is:"))
        print(align.list$Oalign)
        stop(paste("error: a stretch of NAs is expected in the modified field of the ORGAN table at row = ", (idx.+1) ) ) }
      #The new "modified" organ will be aligned to the last reference organ.
      #change the value of NA for the last row of the ORGAN table
      align.list$Oalign$modified[NrO]=idx
      #No shifts necessary
      #change $segmentation:
      #add the "chop/tail" flag (it must be checked that this level belongs to the levels of $segmentation factor)
      if ( "chop/tail" %in% levels(align.list$Oalign$segmentation) ) {
        align.list$Oalign$segmentation[NrO]="chop/tail" } 
      else {
        align.list$Oalign$segmentation=define_seg_levels(align.list$Oalign$segmentation)
        align.list$Oalign$segmentation[NrO]="chop/tail" }
      
      ##Intervals
      #idx is the last organ to append -> there are currently (idx-1) organs in the sequence
      #                                -> there are currently (idx-2) intervals in the sequence
      idxi.=which(align.list$Ialign$modified == (idx-2)) #is the table (Ialign) index of the last $modified interval
      # again, there is a stretch of 'C' (chopped) intervals ending the modified sequence. Reference intervals are aligned to 'NAs'.
      if ( sum(is.na(align.list$Ialign$modified[(idxi.+1):NrI])) != (NrI-idxi.) ){
        print(c("debug: idx is ", idx, ". idxi. (index of the last $modified interval) is ", idxi., ". Aligned sequence of interval is:"))
        print(align.list$Ialign)
        stop(paste("error: a NA is expected in the 'modified' column of the INTERVAL table at row", (idxi.+1) ) ) }
      #change the value of NA for the last row of the INTERVAL table
      align.list$Ialign$modified[NrI]=idx-1 #adding a new interval
      #No shifts necessary
      #dtw code:
      align.list$Ialign$dtw=dtw_code_compute(align.list$Ialign$dtw, idx=NrI, code="T") #code for tail is "T"
    }
  } 
  
 #####
 #Outputs
 #User information message
 No=max(align.list$Oalign$modified, na.rm = TRUE)
  if (verbose){cat(" After edition, the modified sequence now contains", No,"organs.\n")}
 return(align.list)
}
seq_simple_permut=function(align.list, idx, verbose=FALSE){
  ###Description: Create a simple permutation between the two organs corresponding to interval idx of the $modified sequence from input (align.list)
  ##Note: a simple permutation will be performed only between error-free organs (only "~" and "perm" allowed in organ segmentation levels)
  ##input [align.list]: a list made of two dataframes (Ialign and Olign) as element (see specifications)
  ##input [idx]: the index of the INTERVAL whose two organs should be permuted in the $modified sequence.
  
  #####
  #input checks
  #####
  align.list=check_align_list(align.list)
  #Get the number of organs for the input modified sequence
  Ni=max(align.list$Ialign$modified, na.rm = TRUE)
  if (idx < 1 || idx > Ni){
    stop(paste0("input index out of range: it should be in [1-", Ni,"]"))
  }
  
  #####
  #function body
  #####
  #In $modified seq, interval idx is NOT necessarily between organ [idx ; idx+1] because indexes may already be shuffled by permutations
  #Use $Ialign to select $modified intervals only with ~ or P dtw code (possibly multiple P like PPP)
  idx.i=which(align.list$Ialign$modified == idx)
  if (length(idx.i)==0){stop(paste("Error: cannot find interval n°", idx, "in input. Check intergrity of the data."))}
  if (length(idx.i)>1){
    #typical SPLIT case: the modified interval is 'duplicated' in front of two reference intervals
    if (verbose){cat("Organs already affected by segmentation errors, skipping permutation. \n")}
    return(align.list) } 
  else { #idx.1 is necessary length==1
    dtw.code=align.list$Ialign$dtw[idx.i]
    dtw.code=gsub('([[:alpha:]])\\1+', '\\1', dtw.code) #suppress repeated character in the string code
    if (!(dtw.code %in% c("~", "P"))){
      #if dtw code contains something else that just "~" or "P(x)"
      if (verbose){cat("Organs already affected by segmentation errors, skipping permutation. \n")}
      return(align.list) }
    else {
      #Interval: modified the code for the interval and the interval before/after
      if (idx.i==1){#modifying the first interval
        align.list$Ialign$dtw=dtw_code_compute(align.list$Ialign$dtw, idx=c(idx.i, idx.i+1), code = "P") #code for permutation is "P"
      } else if (idx == Ni) { #modifying the last interval of the modified sequence
        align.list$Ialign$dtw=dtw_code_compute(align.list$Ialign$dtw, idx=c(idx.i-1, idx.i), code = "P") #code for permutation is "P"
      } else {
        align.list$Ialign$dtw=dtw_code_compute(align.list$Ialign$dtw, idx=c(idx.i-1, idx.i, idx.i+1), code = "P") #code for permutation is "P"  
      }
      
      #Organs: Permute the indexes of the organ forming this interval idx
      #Row indexes of df $Ialign and $Oalign are linked: 
      #interval idx of row idx.i in $Ialign corresponds to organs [idx;idx+1] of rows [idx.i;idx.i+1] in $Oalign
      before=align.list$Oalign$modified[c(idx.i, idx.i+1)] #store previous data
      align.list$Oalign$modified[c(idx.i, idx.i+1)]=rev(before) #permute
      
      #Organs: mention the operation in the segmentation
      align.list$Oalign$segmentation=define_seg_levels(align.list$Oalign$segmentation)
      align.list$Oalign$segmentation[c(idx.i, idx.i+1)]="perm"
      
      #Return the modified list
      return(align.list)
    }
  }
}

########################################################################
## General wrapper for segmentation errors and measure permutations  ###
########################################################################
segmentation_errors=function(in.seq, align.list, organ_loss=NULL, organ_gain=NULL, verbose=FALSE){
  ## DESCRIPTION: compute the consequence of a of input segmentation errors (organ gains and/or losses) on value sequences and alignment of organs and intervals
  #input [in.seq]: sequences of values for $angles and $internodes
  #input [align.list]: a list made of two dataframes (Ialign and Olign) as element (see specifications) 
  #input [organ_loss]: indexes (vector if several) of the REFERENCE organ(s) to remove
  #input [organ_gain]: indexes (vector if several) of the current REFERENCE organ(s) that will be just AFTER the new insertion.
  #output a list made of three elements: value seq, Ialign and Oalign after the segmentation errors
  
  #Input checks
  align.list=check_align_list(align.list, verbose=TRUE) #Impose verbose here to detect if there are issues with align.list input
  
  #Function's main body
  #compute a table of steps
  O.G=data.frame(idx=organ_gain, type=rep("gain", length(organ_gain)))
  O.L=data.frame(idx=organ_loss, type=rep("loss",length(organ_loss)))
  #a true organ cannot be lost several times ! So remove duplicate if any in 0.L
  O.L=O.L[!duplicated(O.L$idx),]
  steps=rbind(O.G, O.L)
  if(nrow(steps)==0){
    if (verbose){
      warning("no segmentations implemented, no organ gains or losses found from input data")
    }
    return(list(values=in.seq, I=align.list$Ialign, O=align.list$Oalign))
  }
  
  steps=steps[order(steps$idx, steps$type),] #gain is ordered before loss because it was defined before in rbind
  nb_steps=nrow(steps)
  steps=cbind.data.frame(steps$idx, steps$idx, steps$type)
  names(steps)=c("ref.idx", "modif.idx", "type")
  #and what if there are already segmentation errors ? -> conversion needed ?
  
  for (s in 1:nb_steps){
    if (verbose){
      print(paste("step n°", s ,"; type=", steps$type[s], "idx = ", steps$modif.idx[s]))
      print(steps)
    }
    if (steps$type[s] == "gain") {
      align.list=seq_insert(align.list, steps$modif.idx[s], verbose=verbose) 
      in.seq=value_merge(in.seq, steps$modif.idx[s], verbose=verbose) 
      
      #modify steps indexes
      steps$modif.idx=steps$modif.idx + 1
    }
    else {# type == "loss"
      align.list=seq_remove(align.list, steps$modif.idx[s], verbose=verbose) 
      if (steps$modif.idx[s] == 1){
        in.seq=value_chops(in.seq, steps$modif.idx[s], 0)
      }
      else if (steps$modif.idx[s] == max(align.list$Oalign$modified, na.rm = TRUE)) {
        in.seq=value_chops(in.seq, 0, 1)
      }
      else {
        in.seq=value_split(in.seq, steps$modif.idx[s])
      }
      
      #modify steps indexes
      steps$modif.idx=steps$modif.idx - 1
    }
    if (verbose){
      print("end of loop iteration")
      print(align.list)
      print(in.seq)}
  }

  #Outputs
  return(list(values=in.seq, I=align.list$Ialign, O=align.list$Oalign))
}

simple_measure_permutation=function(in.seq, align.list, i_threshold=2, proba=1, verbose=FALSE){
  ##DESCRIPTION: perform simple permutations (= single 2-permutations) on organs of a sequence
  #Created to 'mimic' permutations generated by measurement, manual or automatic (different from more complex biological permutations)
  #it changes values (in.seq) and alignment (align.list)
  #simple permutations is applied with a given likelihood [proba] to organs that met conditions:
  #                   -the internode is below a certain threshold [i_threshold] (measure permutations never occur when organ are well separated on the stem)
  #                   -the two organs involved are not previously affected by segmentation errors
  #It is recommended to use it after 'segmentation_errors' (it will avoid permutations on lost/over segmented organs)
  #input [in.seq]: sequences of values for $angles and $internodes
  #input [align.list]: a list made of two dataframes (Ialign and Olign) as element (see specifications) 
  #input [i_threshold]: numeric, internode length under which two organs can be permuted
  #input [proba]: numerical in range [0-1], likelihood that a pair of organ meeting criteria will be effectively permuted
  #output: a list made of three elements: value seq, Ialign and Oalign after the segmentation errors
  
  #Input checks
  align.list=check_align_list(align.list, verbose=TRUE) #Impose verbose here to detect if there are issues with align.list input
  
  #Function's main body
  #Identify the short internodes that could be permuted
  #1. Criterium 1: short internodes
  idx_short=which(in.seq$internodes <= i_threshold)
  if (length(idx_short) == 0){
    #stop the function, there are no organ pairs to permute
    if (verbose){cat(paste("There are no intervals below the treshold value, no permutations performed, input unchanged.\n"))}
    return(list(values=in.seq, I=align.list$Ialign, O=align.list$Oalign)) } 
  else {
    #2. Criterium 2: no previous segmentation errors 
    idx_perm=c() #initialize
    #For each interval, identify involved organs
    for (i in 1:length(idx_short)){
      idx.o=which(align.list$Oalign$modified == idx_short[i] ) #identify the first organ of the pair corresponding to interval idx_short[i]
      if ( (align.list$Oalign$segmentation[idx.o] %in% c("~", "perm")) && (align.list$Oalign$segmentation[idx.o+1] %in% c("~", "perm")) ){
        #both organs are error-free, the pair can be permuted
        idx_perm=c(idx_perm, idx_short[i]) }
    }
    if (verbose){cat(paste("Number of short error-free internodes that can be permuted =",
                           length(idx_perm),".\n"))}
    if (length(idx_perm)==0){
      #stop the function, there are no organ pairs to permute
      if (verbose){cat(paste("There are no intervals meeting the criteria for permutations, input unchanged.\n"))}
      return(list(values=in.seq, I=align.list$Ialign, O=align.list$Oalign)) }
    else {
      #Loop on these target 'short + error-free intervals'
      count=0 #record the number of effective permutation that will be performed
      for (i in 1:length(idx_perm)){
        #In case of successive permutation, verify that an organ is not drifting away at a larger distance than i_threshold
        row.i=which(align.list$Ialign$modified == idx_perm[i]) #row index of $Ialign containing the value idx_perm
        row.o1=row.i #By construction, property linking Ialign and Oalign: the row index of the first organ of an interval in $Oalign is the same that row index of that interval in $Ialign
        idx.o1=align.list$Oalign$modified[row.o1] #value of organ index of the 1st organ of interval idx_perm
        row.o2=which(align.list$Oalign$modified == idx.o1+1) #Find where in $Oalign is 'o2', the organ that "naturally" (ie without permutation) follows o1 
        row.i2=which(align.list$Ialign$modified == row.o2) #row.i2 is the row of the interval starting with o2 in $Ialign
        idx_i2=align.list$Ialign$modified[row.i2] #value of the interval starting with o2
        if (length(idx_i2)==0){#Particular case when o2 is the last organ, it does not start a next interval
          #change the value to avoid error hereafter and to allow the permutation
          idx_i2=idx_perm[i]}
        
        if ( (idx_i2 < idx_perm[i]) && (sum(in.seq$internodes[c(idx_i2:idx_perm[i])]) > i_threshold)  ){
          #Case1: skip permutation to avoid drift beyond i_threshold? Explanations of if conditions:
          #-Case occurs only if idx2<idx_perm[i] (if idx_i2 > idx_perm: either o2 is just after or even further (hence they will get closer after permutation), the permutation can be done)
          #-conditions whether after permutation, the final distance between o1 and o2 exceeds i_threshold
          if (verbose){cat(paste("Permutation skipped to avoid excessive drift: ($modified) organ n°",
                                 idx.o1,"too far away from organ n°", idx.o1+1,".\n"))}
        } 
        else { #Permutation executed
          test_proba=rbinom(1,1,proba) #generates randomly 0 or 1 depending on proba value
          if (test_proba==1){
            if (verbose){cat(paste("Permuting organs of interval n°",idx_perm[i],"in the $modified sequence.\n"))}
            #Alignment:
            before=as.character(align.list$Ialign$dtw) #to record if changes happen
            align.list=seq_simple_permut(align.list, idx = idx_perm[i], verbose=verbose)
            #Values: only change them if align.list has been modified
            if (  sum(as.character(align.list$Ialign$dtw) == before) !=  nrow(align.list$Ialign) ) {
              in.seq=value_simple_permut(in.seq, idx=idx_perm[i])
              count=count+1
            }
          }
        }
      }
      if (verbose){cat(paste0("Number of permutation performed = ", count, ".\n"))}
      #Outputs
      return(list(values=in.seq, I=align.list$Ialign, O=align.list$Oalign))
    }
  }
}

#########################################
##  functions to provide information   ##
#########################################
print_info=function(seg_errors, permutation, Noise_or_Measures){
  cat("reminder of main scenario parameters \n")
  print(paste("seg_errors =",seg_errors))
  print(paste("permutation =", permutation))
  print(paste("Noise or measures =", Noise_or_Measures))}