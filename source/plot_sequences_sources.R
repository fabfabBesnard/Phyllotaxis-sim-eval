###########################################################################
#### Copyright (C) 2021 - INRAe (Fabrice Besnard, RDP)
#### This script is a free software: you can redistribute it
#### and/or modify it under the terms of the GNU General Public
#### License (GNU GPL) as published by the Free Software Foundation, either
#### version 3 of the License, or (at your option) any later version.
#### Distributed without any warranty.
###########################################################################
#started 2020-09-20
# last edit: 2021-11-17
#Version v0

##Content
#Plotting one sequence 'seq_plot'
#plotting two sequences: 'multiseq_plot'
#plotting two alignements (e.g. ground truth and prediction) of two sequences: 'compare_plots'

##################
##   Libraries  ##
##################
pkgTest <- function(x){
  if (!require(x,character.only = TRUE))
  {
    install.packages(x,dep=TRUE)
    if(!require(x,character.only = TRUE)) stop("Package not found")
  }
}
suppressPackageStartupMessages(require(ggplot2))
suppressPackageStartupMessages(require(reshape2))
pkgTest("gridExtra")
suppressPackageStartupMessages(require(gridExtra))

#################################
## Plotting a single sequence  ##
#################################
seq_plot=function(seq, angle_ref=137){
  #Plot on top of each other angles and internode values (y-axis) function of interval order (xaxis)
  #seq must be a dataframe containing at least three fields: $intervals / $angles / $internodes
  #angle_ref [input]: value of the horizontal dashed line plotted in the divergence angle plot (default=137°)
  
  seq.long=melt(seq, id.vars = c("intervals"))
  #grid.arrange solution:
  p1=ggplot(subset(seq.long, variable=="angles"), aes(x=intervals, y=value))+
    geom_point()+geom_line()+
    geom_hline(yintercept=angle_ref, linetype="dashed")+
    ylim(0,360)+ylab("ANGLES")+xlab(NULL)+
    theme_minimal()+
    theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)))
  
  p2=ggplot(subset(seq.long, variable=="internodes"), aes(x=intervals, y=value))+
    geom_point()+geom_line()+
    ylab("INTERNODES")+xlab("interval order (base -> top)")+
    theme_minimal()+
    theme(axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0)))
  
  grid.arrange(p1, p2, nrow=2)
}

###########################################################################
# Plotting several unaligned sequences / or 2 sequences with realignment
###########################################################################
check_format_seq_df=function(seq.df, force.convert=TRUE, verbose=FALSE){
    #DESCRIPTION: check that a sequence dataframe is formatted with the following three fields
    #$intervals, $angles, $internodes
  #seq.df [input]: dataframe of phyllotaxis sequences for a plant
  #force.convert [input]: TRUE/FALSE, if TRUE try to convert the input into the expected format
  #verbose [input]: TRUE/FALSE, increase verbosity
  
  #Checking the format = dataframe
  if (!is.data.frame(seq.df)){
    if (force.convert){
      warning(deparse(substitute(seq.df)), "is not a dataframe. Forcing conversion")
      seq.df=as.data.frame(seq.df)}
    else {stop(deparse(substitute(seq.df)), "is not a dataframe. Aborting program")}
  }
  
  #Checking the field identity
  if (sum(colnames(seq.df) == c("intervals", "angles", "internodes")) != 3){
    #fields are different than the classical input from R-simulated data -> go for in-depth checks
    intervals.colidx=grep("intervals", colnames(seq.df), ignore.case = TRUE)
    intervals.colidx=ifelse(length(intervals.colidx), intervals.colidx, 0)
    angles.colidx=grep("angles", colnames(seq.df), ignore.case = TRUE)
    angles.colidx=ifelse(length(angles.colidx), angles.colidx, 0)
    internodes.colidx=grep("internodes", colnames(seq.df), ignore.case = TRUE)
    internodes.colidx=ifelse(length(internodes.colidx), internodes.colidx, 0)
    
    found.cols=c(intervals.colidx, angles.colidx, internodes.colidx)
    #print(found.cols)
    
    if (sum(found.cols == c(0,0,0))==3){
      #no expected fields were found, try to convert them all based on column positions
      if (ncol(seq.df)==3 & force.convert){ #only try the conversion if three input columns
        if (verbose) {
          warning("Unexpected colnames for ", deparse(substitute(seq.df)), ". They will be changed from ", paste0(colnames(seq.df), ", "), "to (intervals, angles, internodes)")}
        #change the names of the columns
        colnames(seq.df) = c("intervals", "angles", "internodes")
        #make sure that intervals are intervals:
        if (!is.integer(seq.df$intervals)){
          warning("Inappropriate format for intervals. True intervals will be computed")
          seq.df$intervals=1:nrow(seq.df)
        }
        }
      else { stop(paste("Expected field names 'intervals/angles/internodes' were not found in", deparse(substitute(seq.df)), "Aborting program")) }
    } else {
      if (sum(found.cols != c(0,0,0))==3){
        #all fields are found. They may be either case variations in field name or change(s) in their order inside the dataframe.
        #rearrange:
        if (verbose){ 
          print(paste("Unusual field names were found for", deparse(substitute(seq.df)), ":"))
          print(colnames(seq.df))}
        if (force.convert){
            warning("The three expected columns 'intervals/angles/internodes' were automatically guessed, please check that no errors mare made.")
          seq.df=rbind.data.frame(seq.df[,found.cols[1]], seq.df[,found.cols[2]], seq.df[,found.cols[3]])
          colnames(seq.df)=c("intervals", "angles", "internodes") }
        else { #do not try force conversion and abort
          stop("Aborting program") } } 
      else {
        #Only consider two cases: both ~angles and ~internodes have been found (and so you try to convert) or no such fields are found and so the program abort
        if (sum(found.cols[c(2,3)] != c(0,0))==2){
          if (verbose){ 
            print(paste("Unusual field names were found for", deparse(substitute(seq.df)), ":"))
            print(colnames(seq.df))}
          if (force.convert){
            if (verbose) {warning("Expected 'interval' column was not found in ", deparse(substitute(seq.df)), ": adding a proper interval column (angles and internodes values were OK)")}
            #the two fields angles and internodes have been found, just add intervals
            seq.df=cbind.data.frame(intervals=1:nrow(seq.df), seq.df[,found.cols[2]], seq.df[,found.cols[3]])
            colnames(seq.df)=c("intervals", "angles", "internodes") }
          else{
            stop(paste("Expected field names 'intervals/angles/internodes' were not found in", deparse(substitute(seq.df)), "Aborting program"))}
          }
        else {
          #in all other cases, just abort the program
          stop(paste("Expected field names 'intervals/angles/internodes' were not found in", deparse(substitute(seq.df)), "Aborting program"))}
      }
      
      }
  }
  return(seq.df)
}#end of function check_format_seq_df

multiseq_plot=function(mylist, align.df=NULL, prediction.eval= NULL, 
                       id.names=NULL, ref.first=TRUE, title=NULL, verbose=FALSE, 
                       angle_ref=137){
  #DESCRIPTION: plot multiple sequences on the same plot, like the "true" plant and one or several measure of it
  #input [mylist]: a list of phyllotaxis sequence df ($intervals, $angles, $internodes).
                      #In the case of two sequences representing a reference and a realigned sequence, give the reference first.
  #input [align.df]: a dataframe giving the interval alignment (output of function 'segmentation_errors'). 
                      #If given: the value of the modified sequence will plotted aligned to the reference sequence ;
                      #if given: the predicted sm-dtw code for alignment  will be written on top of each realigned value of the test sequence
                      #/!\ only works for 2 sequences given mylist
  #input [prediction.eval]: a dataframe containing the evaluation of prediction alignment for each indices of the test sequence
                      # (given by 'evaluate_align_prediction' function in eval_dtw_source.R)
                      # This will highlight in RED the interval where the prediction is evaluated as error and in ORANGE where the prediction is ambiguous
  #input [id.names]: a list of the name you want to give to the sequences in the same order as input mylist. 
  #By default, it takes the names of the elements of the list if they exist. If not, the first sequence is considered as the reference
  #angle_ref [input]: value of the horizontal dashed line plotted in the divergence angle plot (default=137°)
  
  ##Inputs
  if (is.null(id.names)){
    id.names=names(mylist)
    if (is.null(id.names)){
      if (ref.first){
        id.names=c("reference", paste0("seq", 2:length(mylist)))
      } else {
        id.names=c(paste0("seq", 1:length(mylist)))
      }
    }
  }
  if (!is.null(align.df)){
    if (verbose){
      cat("Option 'align.df' selected. In the input list of sequence values, reference must be first and the test sequence must be second. \n ") 
    }
    if(!is.data.frame(align.df) || 
       sum(c("reference","modified", "dtw") %in% colnames(align.df)) != 3){
      stop("Issue with align.df format. Program aborts")}
    if (length(mylist)>2){
      warning("Option align.df only works for 2 input sequences. The alignment will be ignored")
      align.df=NULL
    }
  }
  
  #Main body
  #####
  #1. Create a dataframe df containing all sequences, with the last field =$id (ref & test)
  #####
  #Initialize with the first element of the list:
  df=check_format_seq_df(mylist[[1]], force.convert=TRUE, verbose = verbose)
  df$id=id.names[1]
  if (is.null(align.df)){
    #loop over the rest of the list
    for (i in 2:length(mylist)){
      #tempdf=mylist[[i]]
      tempdf=check_format_seq_df(mylist[[i]], force.convert=TRUE, verbose = verbose)
      tempdf$id=id.names[i]
      df=rbind(df, tempdf)
    } }
  else { #2 sequences only, one is realigned against the other:
    df2=check_format_seq_df(mylist[[2]], force.convert=TRUE, verbose = verbose)
    #1. select modified values that are aligned to a reference interval
    aligned=df2[align.df$modified[which(!is.na(align.df$reference))],]
    #Replace the $modified interval by the reference intervals
    aligned$intervals=align.df$reference[which(!is.na(align.df$reference))]
    #add the name for this plot
    aligned$id=id.names[2]
    #add the dtw for this plot
    df$dtw="" #empty to initialize
    aligned$dtw=align.df$dtw[which(!is.na(align.df$reference))]
    if (!is.null(prediction.eval)){
      if (verbose){
        cat("An evaluation of the predicted realignment has been given. Errors will be highlighted in the plot. \n ")
      }
      #add the evaluation to the alignment df, just after the code:
      align.df$eval=prediction.eval$dtw.eval[align.df$modified]
      #also add this new field to the joined dataframe of ref/test
      df$eval="" #empty initialization
      aligned$eval=align.df$eval[which(!is.na(align.df$reference))]
    }
    df=rbind(df, aligned) #join data
    #2. record tails separately if any
    tail_idx_ref=which(is.na(align.df$reference)) # tails are "test" values for which no reference values (NA) are aligned
    if (length(tail_idx_ref)>0){
      #store a starting tail if any
      if (min(tail_idx_ref) < which(align.df$reference == 1) ) {
        tail_idx=tail_idx_ref[tail_idx_ref< which(align.df$reference == 1)]
        start_tail=df2[align.df$modified[tail_idx],]
        #change the indexes
        start_tail$intervals=sort(-seq(1:length(tail_idx))) #by convention, starting tails will end up at -1 to be plotted separated from reference values
        #add a name
        start_tail$id=paste0(id.names[2],"/starting_tail")
        #dtw
        start_tail$dtw=align.df$dtw[tail_idx]
        #eval
        if (!is.null(prediction.eval)){
          start_tail$eval=align.df$eval[tail_idx]
        }
        #bind
        df=rbind(df, start_tail)
        #suppress the starting tail indexes from record before analyzing trailing tail
        tail_idx_ref=tail_idx_ref[!(tail_idx_ref %in% tail_idx)]
      }
      
      #store a trailing tail if any
      if ( (length(tail_idx_ref)>0) &&
           (min(tail_idx_ref) > which(align.df$reference == max(align.df$reference, na.rm = TRUE)) ) ) {
        tail_idx= tail_idx_ref[tail_idx_ref > which(align.df$reference == max(align.df$reference, na.rm = TRUE))]
        trailing_tail=df2[align.df$modified[tail_idx],]
        #replace the intervals (modified) by a index realigned against the reference
        trailing_tail$intervals=seq((max(align.df$reference, na.rm = TRUE)+1), #trailing tails will start after the last reference interval, so they appear separated from reference values in the plot
                                    (max(align.df$reference, na.rm = TRUE)+length(tail_idx)), 1)
        #add a name
        trailing_tail$id=paste0(id.names[2],"/trailing_tail")
        #dtw
        trailing_tail$dtw=align.df$dtw[tail_idx]
        #eval
        if (!is.null(prediction.eval)){
          trailing_tail$eval=align.df$eval[tail_idx]
        }
        #bind
        df=rbind(df, trailing_tail)
      }
    }
  }
  #print(df)
  
  #####
  #2. Transform to a long df compatible with ggplot2
  #####
  if (is.null(align.df)){
    seq.long=melt(df, id.vars = c("intervals", "id"))
  }
  else {
    if(is.null(prediction.eval)){
      seq.long=melt(df, id.vars = c("intervals", "id", "dtw"))  
    } else {
      seq.long=melt(df, id.vars = c("intervals", "id", "dtw", "eval"))  
    }
  }
  #print(seq.long)
  
  #####
  ##3. PLOT
  #####
  # print(str(seq.long))
  # print(levels(as.factor(seq.long$id)))
  nb.id=length(levels(as.factor(seq.long$id)))
  if (!is.null(align.df)){
    #typical case of 2 sequences, one is ref (assumed to be first element in input mylist), the other is test
    #re-order levels to put reference first, followed by realigned sequences and its starting/trailing tails if any:
    seq.long$id=as.factor(seq.long$id)
    other.id=levels(seq.long$id)[levels(seq.long$id) != id.names[1]]
    seq.long$id=factor(seq.long$id, levels=c(id.names[1], other.id))
  }
  
  p1=ggplot(subset(seq.long, variable=="angles"), aes(x=intervals, y=value, color=id))+
    geom_point()+geom_line()+
    geom_hline(yintercept=angle_ref, linetype="dashed")+
    scale_color_manual(values = c("#F8766D", rep("#00BFC4", nb.id-1) ) ) +
    ylim(0,360)+ylab("ANGLES")+xlab(NULL)+
    theme_minimal()+
    theme(legend.position="bottom",
          axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)))
  p2=ggplot(subset(seq.long, variable=="internodes"), aes(x=intervals, y=value, color=id))+
    geom_point()+geom_line()+
    scale_color_manual(values = c("#F8766D", rep("#00BFC4", nb.id-1) ) ) +
    ylab("INTERNODES")+xlab("interval order (base -> top)")+
    theme_minimal()+
    theme(axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0)))
  
  #extract legend
  #https://stackoverflow.com/questions/13649473/add-a-common-legend-for-combined-ggplots
  #https://github.com/hadley/ggplot2/wiki/Share-a-legend-between-two-ggplot2-graphs
  g_legend<-function(a.gplot){
    tmp <- ggplot_gtable(ggplot_build(a.gplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    return(legend)}
  mylegend<-g_legend(p1)
  
  #extract a common title
  g_title<-function(a.gplot){
    tmp <- ggplot_gtable(ggplot_build(a.gplot))
    title_grob_idx <- which(tmp$layout$name == "title")
    title_grob <- tmp$grobs[[title_grob_idx]]
    return(title_grob)}
  
  #if align.df is not NULL, plot the dtw code
  if (!is.null(align.df)){
    #Trick to display chop "C" (NA values in the test sequences, so the label is not printed)
    #Replace values of the test sequence by the value of the reference sequence
    seq.long$dtw=as.factor(seq.long$dtw)
    if ("C" %in% levels(seq.long$dtw) ){
      chops.df=seq.long[seq.long$dtw == "C",]
      match.chop.ref=seq.long[(seq.long$id != id.names[2]) & (seq.long$intervals %in% chops.df$intervals), ]
      seq.long.print=seq.long
      seq.long.print[seq.long.print$dtw == "C",]$value=match.chop.ref$value
    }
    else { seq.long.print=seq.long}
    p1=p1+
      geom_text(data=subset(seq.long.print, variable=="angles"),
                aes(y=value, label=dtw, vjust=0), size=3, color="black")
    p2=p2+
      geom_text(data=subset(seq.long.print, variable=="internodes"), 
                aes(y=value, label=dtw, vjust=0), size=3, color="black")
  }
  
  #if prediction.eval is not NULL, highlight the possible errors
  if (!is.null(prediction.eval)){
    #evaluation is the same for angles or internodes: take one of them
    errors.highlight=seq.long[(seq.long$eval != "correct" & 
                                 seq.long$variable == "angles" & 
                                 seq.long$id != id.names[1]) ,] #case of ref.first=TRUE
    #suppress NA (coming from Chops in the reference intervals):
    errors.highlight=errors.highlight[!is.na(errors.highlight$eval),]
    #suppress duplicated coordinates coming from merge:
    errors.highlight=errors.highlight[!duplicated(errors.highlight$intervals),]
    p1=p1+
      geom_rect(data = errors.highlight, xmin = errors.highlight$intervals - 0.5, xmax = errors.highlight$intervals +
                  0.5, ymin = -Inf, ymax = Inf, colour=NA, fill = 'red', alpha = 0.1)
    internode_values=subset(seq.long.print, variable=="internodes")
    ymax=max(internode_values$value)
    p2=p2+
      geom_rect(data = errors.highlight, xmin = errors.highlight$intervals - 0.5, xmax = errors.highlight$intervals +
                  0.5, ymin = -Inf, ymax = Inf, colour=NA, fill = 'red', alpha = 0.1)+
      scale_y_continuous(limits = c(0,1.2*ymax))
    
  }
  
  #Final layout of the two plots (angles/internodes on top of each other)
  if (!is.null(title)){
    mytitle=g_title(p1+ggtitle(title)+
                      theme(plot.title = element_text(hjust = 0.5, size=12)))
    grid.arrange(mytitle,
                 arrangeGrob(p1 + theme(legend.position="none"),
                             p2 + theme(legend.position="none"),
                             nrow=2),
                 mylegend, 
                 nrow=3, heights=c(1, 10, 1))
  }
  else {
    grid.arrange(arrangeGrob(p1 + theme(legend.position="none"),
                             p2 + theme(legend.position="none"),
                             nrow=2),
                 mylegend, nrow=2,heights=c(10, 1))
  }
  #return(seq.long)
}

multiseq_plot_pdf=function(seq.ref, seq.test, 
                           true.align, id.names=NULL,
                           PlantID=NULL, pdf.name="Simulated_paired_sequences.pdf", verbose=FALSE){
  ## DESCRIPTION: outputs plots of simulated paired sequences in a pdf
  # seq.ref [input]: dataframe formated as $PlantID $angles $internodes, contains the values of the reference sequence
  # seq.test [input]: dataframe formated as $PlantID $angles $internodes, contains the values of the test sequence to be aligned against the reference sequence
  # true.align [input]: dataframe formated as $PlantID $reference $modified $dtw
  # #input [id.names]: a list of the name you want to give to the related sequences in the same order as input mylist, e.g: "reference" and "test"
  # PlantID [input]: optional, name(s) of a PlantID(s) (if several given as c("id1", "id2")) if only the plot for specific plants is wanted. Default=NULL, all PlantIDs will be processed
  # pdf.name [input]: name of the saved pdf
  # verbose [input]: increase verbosity
  #Checking/reformatting Inputs
  seq.ref$PlantID=as.factor(seq.ref$PlantID)
  seq.test$PlantID=as.factor(seq.test$PlantID)
  true.align$PlantID=as.factor(true.align$PlantID)
  if (!is.null(PlantID)){ #PlantID must be in the four input dataframes
    if (sum(PlantID %in% levels(seq.ref$PlantID)) != length(PlantID)  &
        sum(PlantID %in% levels(seq.test$PlantID)) != length(PlantID) &
        sum(PlantID %in% levels(true.align$PlantID)) != length(PlantID)){
      stop("Given PlantID is/are not found in all input dataframes. Check PlantID and/or input data.")
    }
  }
  
  #Main function body
  if (is.null(PlantID)){ plotted_ids=levels(seq.ref$PlantID) } else { plotted_ids=PlantID }
  if (length(plotted_ids==1)){pop_up_window=TRUE}
  plot_page=list() #initialize the list of plots
  i=1 #initialize a counter for the above list length
  if (verbose){
    cat("plots will be generated for: \n")
    cat(plotted_ids, "\n")  
  }
  for (id in plotted_ids){
    if (verbose) {print(paste("starting for", id))}
    #Subset per PlantID and rearrange them to use the function `multiseq_plot` 
    seq.ref.id=seq.ref[seq.ref$PlantID==id, ] #subset rows corresponding to the current PlantID
    seq.ref.id=check_format_seq_df(seq.ref.id, verbose = verbose) #add the 'interval' row and remove 'PlantID' row
    seq.test.id=seq.test[seq.test$PlantID==id, ] #same for seq.test
    seq.test.id=check_format_seq_df(seq.test.id, verbose = verbose)
    true.align.id=true.align[true.align$PlantID==id, ] #same subset for ground truth alignment
    
    #plot `multiseq_plot`
    p1=multiseq_plot(list(seq.ref.id, seq.test.id), 
                     align.df = true.align.id, 
                     id.names=id.names,
                     title=paste0("PlantID: ", id, "(simulated paired sequences)") )
    plot_page[[i]]=p1
    i=i+1 #update the next index of the list
    
  }#end of for loop
  
  multipage_plots <- marrangeGrob( plot_page, nrow=1, ncol=1)
  ggsave(file=pdf.name, multipage_plots, width = 210, height = 297, units = "mm")
}

#############################################################################
# Plotting ground truth alignments and prediction alignements (in one pdf page)
#############################################################################
compare_plots=function(seq.ref, seq.test, 
                       true.align, dtw.results, prediction.eval= NULL,
                       PlantID=NULL, 
                       id.names=NULL,
                       PDF=TRUE, pdf.name="Compare_Prediction_Plots.pdf", verbose=FALSE){
  ## DESCRIPTION: display simulation & predicted re-alignment on the same screen : on the top, the true alignment (simulated) between a test and a reference sequences ; on the bottom the predicted alignment between the same sequences by a program (e.g. dtw)  
  # seq.ref [input]: dataframe formated as $PlantID $angles $internodes, contains the values of the reference sequence
  # seq.test [input]: dataframe formated as $PlantID $angles $internodes, contains the values of the test sequence to be aligned against the reference sequence
  # true.align [input]: dataframe formated as $PlantID $reference $modified $dtw
  # dtw.results [input]: dataframe formated as $PlantID $intervals $reference $test $dtw $cost, output from function `convert_dtw_results` 
  # prediction.eval [input]: a dataframe containing the evaluation of prediction alignment for each indices of the test sequence
  # (given by 'evaluate_align_prediction' function in eval_dtw_source.R)
  # This will highlight in RED the interval where the prediction is evaluated as error and in ORANGE where the prediction is ambiguous
  # PlantID [input]: optional, name(s) of a PlantID(s) (if several given as c("id1", "id2")) if only the plot for specific plants is wanted. Default=NULL, all PlantIDs will be processed
  # PDF [input]: TRUE/FALSE: indicate whether a pdf of the plots should be printed (one PlantID per page)
  # pdf.name [input]: name of the saved pdf
  # verbose [input]: increase verbosity
  
  #Checking/reformatting Inputs
  seq.ref$PlantID=as.factor(seq.ref$PlantID)
  seq.test$PlantID=as.factor(seq.test$PlantID)
  dtw.results$PlantID=as.factor(dtw.results$PlantID)
  if (!is.null(PlantID)){ #PlantID must be in the four input dataframes
    if (sum(PlantID %in% levels(seq.ref$PlantID)) != length(PlantID)  &
        sum(PlantID %in% levels(seq.test$PlantID)) != length(PlantID) &
        sum(PlantID %in% levels(true.align$PlantID)) != length(PlantID) &
        sum(PlantID %in% levels(dtw.results$PlantID)) != length(PlantID) ){
      stop("Given PlantID is/are not found in all input dataframes. Check PlantID and/or input data.")
    }
  }
  
  #Main function body
  if (is.null(PlantID)){ plotted_ids=levels(seq.ref$PlantID) } else { plotted_ids=PlantID }
  if (length(plotted_ids==1)){pop_up_window=TRUE}
  sideside_plots=list() #initialize the list of plots
  i=1 #initialize a counter for the above list length
  if (verbose){
    cat("plots will be generated for: \n")
    cat(plotted_ids, "\n")  
  }
  for (id in plotted_ids){
    if (verbose) {print(paste("starting for", id))}
    #Subset per PlantID and rearrange them to use the function `multiseq_plot` 
    seq.ref.id=seq.ref[seq.ref$PlantID==id, ] #subset rows corresponding to the current PlantID
    seq.ref.id=check_format_seq_df(seq.ref.id, verbose = verbose) #add the 'interval' row and remove 'PlantID' row
    seq.test.id=seq.test[seq.test$PlantID==id, ] #same for seq.test
    seq.test.id=check_format_seq_df(seq.test.id, verbose = verbose)
    true.align.id=true.align[true.align$PlantID==id, ] #same subset for ground truth alignment
    dtw.results.id=dtw.results[dtw.results$PlantID==id, ] #same subset for dtw.results + reformat/select useful fields
    dtw.results.id=cbind.data.frame(reference=dtw.results.id$reference, modified=dtw.results.id$test, dtw=dtw.results.id$dtw)
    if (!is.null(prediction.eval)){
    prediction.eval.id=prediction.eval[prediction.eval$PlantID==id,] #same subset for the evaluation
    }
    
    #plot `multiseq_plot`
    p1=multiseq_plot(list(seq.ref.id, seq.test.id), align.df = true.align.id, 
                     id.names=id.names, 
                     title=paste0("PlantID: ", id, "/ True Simulation") )
    if (is.null(prediction.eval)){
      p2=multiseq_plot(list(seq.ref.id, seq.test.id), align.df = dtw.results.id, 
                       id.names=id.names, 
                       title=paste0("PlantID: ", id, "/ Prediction") )  
    } else {
      p2=multiseq_plot(list(seq.ref.id, seq.test.id), align.df = dtw.results.id, prediction.eval = prediction.eval.id,
                       id.names=id.names, 
                       title=paste0("PlantID: ", id, "/ Prediction") )  
    }
    sideside_plots[[i]]=p1
    sideside_plots[[i+1]]=p2
    i=i+2 #update the next index of the list
  }#end of for-loop
  
  if (PDF){#Saving as pdf
    multipage_plots <- marrangeGrob(sideside_plots, nrow=2, ncol=1)
    ggsave(file=pdf.name, multipage_plots, width = 210, height = 297, units = "mm")
  }
  #return(sideside_plots)
}