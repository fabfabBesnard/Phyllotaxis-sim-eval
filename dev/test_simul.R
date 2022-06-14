local.repo="~/Documents/RDP/MyProjects/ROMI/Data/Eval_AnglesAndInternodes/" #add the final '/'
setwd(paste0(local.repo, "tests"))
source(paste0(local.repo, "Phyllotaxis-sim-eval/source/sim_phyllo_sources.R"))
source(paste0(local.repo, "Phyllotaxis-sim-eval/source/plot_sequences_sources.R"))

#################################
#1. Make reference and plot it
N=19 #Number of intervals in the sequence (which contains N+1 organs)
alpha=137.5 #the average value of divergence angles (137°5 is the value for Fibonacci spiral phyllotaxis)
a_sd=18.5 #biological/natural gaussian noise on angle values
i_Gsd=0.8 #biological/natural gaussian noise on internode values
i_noise_pct=75 #scale the biological/natural noise on internode relative to the average internode length
seq=make_refseq(N, alpha, a_sd, i_Gsd, i_noise_pct)
align_init=make_align_list(N)


REF=make_refseq(N, alpha, a_sd, i_Gsd, i_noise_pct)
seq_plot(REF)

############################
#2. Make a measure of reference and plot it with REF
manual_anoise_sd=20
manual_inoise_sd=5
man_MEAS=make_measure(REF, manual_anoise_sd, manual_inoise_sd)
seq_plot(man_MEAS)

################
#3. compare the two plots
seq.list=list(REF, man_MEAS)
multiseq_plot(seq.list, title = "test")

##################
#4. Make a computer measure 
aut_MEAS=make_measure(REF, aut_anoise_sd, aut_inoise_sd)
multiseq_plot(list(REF, aut_MEAS))
seq.list=list(REF, man_MEAS, aut_MEAS)
multiseq_plot(seq.list)

###################
#5. Introduce sequence errors=SPLIT/MERGE only and plot sequences
#this sequence is modified by Merge & Split
#Set randomly an idx where the error will occur within refseq:
idx=round(runif(1,1,N+2), digits = 0) #modeled as uniform probability around all the sequence
#idx=31

manMEAS.error.M=value_merge(man_MEAS,idx)  #possible idx range is [1;n+2]: the idx taken by the false organ in the new sequence
multiseq_plot(list(REF, manMEAS.error.M))

manMEAS.error.S=value_split(man_MEAS,idx) #possible idx range is [1;n+1]: the idx of the true organ which is missed
multiseq_plot(list(REF, manMEAS.error.S))

manMEAS.chop=value_chops(man_MEAS,5,7)
multiseq_plot(list(REF, manMEAS.chop))

manMEAS.tail=value_chops(man_MEAS,1,1)
multiseq_plot(list(REF, manMEAS.tail))

###################
#6. Introduce segmentation errors
#test sequences
O1=data.frame(reference=c(1:10),
              modified=c(1:10),
              segmentation=c(rep("~",10)))
I1=data.frame(reference=c(1:9),
              modified=c(1:9),
              dtw=c(rep("~",9)))

ref=c(1,2,3,NA,4,5,NA,6,7,8)
modif=c(NA,NA,1,2,3,NA,4,5,6,7)
seg=c(rep("chop",2),"~","over", "~", "under","over", rep("~",3))
O2=data.frame(reference=ref,
              modified=modif,
              segmentation=seg)
I2=data.frame(reference=c(1,2,3,3,4,5,6,7),
              modified=c(NA,NA,1,2,3,4,5,6),
              dtw=c("C", "C", "M", "M", "SM", "SM", "~", "~"))

O3=data.frame(reference=ref,
              modified=c(NA,NA,1,2,3,NA,4,5,NA,NA),
              segmentation=c(rep("chop",2),"~","over", "~", "under","over", "~", rep("chop",2)))
I3=data.frame(reference=c(1:9),
              modified=c(NA,NA,1,2,3,3,4,NA,NA),
              dtw=c("C", "C", "M", "M", "S", "SM", "M", "C", "C"))

O4=data.frame(reference=c(1:7),
              modified=c(1:7),
              segmentation=c(rep("~",7)))
I4=data.frame(reference=c(1:6),
              modified=c(1:6),
              dtw=c(rep("~",6)))

O5=data.frame(reference=c(NA,NA,1:7),
              modified=c(1:9),
              segmentation=c("tail", "tail", rep("~",7)))
I5=data.frame(reference=c(NA,NA,1:6),
              modified=c(1:8),
              dtw=c("T","T",rep("~",6)))

O6=data.frame(reference=c(1:7,NA,NA),
              modified=c(1:9),
              segmentation=c(rep("~",7), "tail", "tail"))
I6=data.frame(reference=c(1:6, NA,NA),
              modified=c(1:8),
              dtw=c(rep("~",6), "T","T"))

O7=data.frame(reference=c(1:7),
              modified=c(NA,1:6),
              segmentation=c("chop", rep("~",6)))
I7=data.frame(reference=c(1:6),
              modified=c(NA, 1:5),
              dtw=c("C", rep("~",5)))

O8=data.frame(reference=c(1:7),
              modified=c(1:6,NA),
              segmentation=c(rep("~",6), "chop"))
I8=data.frame(reference=c(1:6),
              modified=c(1:5,NA),
              dtw=c(rep("~",5), "C"))

O9=data.frame(reference=c(1:7),
              modified=c(1:2, NA, 3:6),
              segmentation=c(rep("~",2), "under", rep("~",4)))
I9=data.frame(reference=c(1:6),
              modified=c(1,2,2,3:5),
              dtw=c("~", "S", "S", rep("~",3)))

#Test lists
test=list(Ialign=I1, Oalign=O1) #10 aligned organs
test2=list(Ialign=I2, Oalign=O2) #dtw=c("C", "C", "M", "M", "SM", "SM", "~", "~")

test3=list(Ialign=I3, Oalign=O3)

test4=list(Ialign=I4, Oalign=O4) # $ref == $modified
test5=list(Ialign=I5, Oalign=O5) # x 2 tails already starting $modified
test6=list(Ialign=I6, Oalign=O6) # x2 tails already at the end of $modified
test7=list(Ialign=I7, Oalign=O7) # 1st organ in ref is already chopped
test8=list(Ialign=I8, Oalign=O8) # last organ in ref is already chopped
test9=list(Ialign=I9, Oalign=O9) # test for a Split/Merge
test10=seq_remove(test9, 3)

#Make lists:
list23=make_align_list(23)

###################
##test the functions
#limit cases: test insertion/removal at the tails
View(test$Ialign)
View(test$Oalign)
seq_remove(test, 1)
seq_insert(test, 1)
seq_insert(test, 2)
seq_append(test,1)
seq_remove(test,10)

seq_insert(test,6)
seq_insert(test2,6)
View(test2$Oalign)

seq_remove(test,6)
newlist=seq_remove(test2,6)

##Test seq_append
View(test3$Oalign)
View(test3$Ialign)
seq_append(test3,1)
seq_append(test3,6)

View(test4$Oalign)
View(test4$Ialign)
seq_append(test4,1) #at start
seq_append(test4,8) #at the end

View(test5$Oalign)
View(test5$Ialign)
seq_append(test5,1)
seq_append(test5,10)

View(test6$Oalign)
View(test6$Ialign)
seq_append(test6,1)

View(test7$Oalign)
View(test7$Ialign)
seq_append(test7,1)

View(test8$Oalign)
View(test8$Ialign)
seq_append(test8,1)
seq_append(test8,7)
###################
#Test SM
View(test9$Oalign)
View(test9$Ialign)
#insert an organ in a gap of size 1: $modified organ n°3 is a under/over
seq_insert(test9,3)
View(test10$Oalign)
View(test10$Ialign)
#insert an organ in a gap of size 2: $modified organ n°3 is a under/over
seq_insert(test10,3)
seq_append(test8,7)

#test split on split
#remove an organ just before an existing gap of size 1:
seq_remove(test9,2)
#and just after -> test10
seq_remove(test9,3)
###################
#Test segmentation_errors
test #10 organs
seq=make_refseq(nrow(test$Ialign), alpha, a_sd, i_Gsd, i_noise_pct)
man_MEAS=make_measure(seq, manual_anoise_sd, manual_inoise_sd)
aut_MEAS=make_measure(seq, aut_anoise_sd, aut_inoise_sd)
multiseq_plot(list(seq, man_MEAS, aut_MEAS))

final=segmentation_errors(aut_MEAS, test, 
                    organ_loss = c(1,2,10), 
                    organ_gain = c(2,3))

#note: modify seq_insert to replace organ if there is a NA cell (instead of creating a new row)

segmentation_errors(aut_MEAS, test, 
                    organ_loss = c(1,2,8,10), 
                    organ_gain = c(4,5))
segmentation_errors(aut_MEAS, test, 
                    organ_loss = c(1,2,8,10), 
                    organ_gain = c(4,5),
                    verbose=TRUE)

AUT=segmentation_errors(aut1, listN1,
                        organ_loss=c(1,2,21,22,23,24),
                        organ_gain=c(6,6,12))

AUT=segmentation_errors(aut1, listN1,
                        organ_loss=c(6,6,12),
                        organ_gain=c(1,1,25,25), verbose = TRUE)

AUT=segmentation_errors(aut1, listN1,
                        organ_loss=c(6,6,12),
                        organ_gain=c(1,1,25), verbose=TRUE) 
AUT=segmentation_errors(aut1, listN1,
                        organ_gain=c(1,1,1), verbose = TRUE)
AUT=segmentation_errors(aut1, listN1,
                        organ_gain=c(25,25,25), verbose = TRUE)
###################
#insert several organs in the same interval:
segmentation_errors(aut_MEAS, test, 
                    organ_loss = c(1,2,8,10), 
                    organ_gain = c(4,4),
                    verbose=TRUE)
comp=segmentation_errors(aut_MEAS, test, 
                         organ_loss = c(1,2,8,10), 
                         organ_gain = c(4,5))

multiseq_plot(list(seq, man_MEAS, aut_MEAS, comp$values), 
              id.names = c("real_biol", "manual", "automated", "with_seg_errors"))


write.csv("test",AUT$values)

##Test input data file
data_test=read.delim2("dtw_input_plants.csv", header=TRUE)
###################
##test for realigned plots
####
# EXAMPLES derived from experiments
#ex: Col0_26_10_2018_B
N1=23
seq1=make_refseq(N1, alpha, a_sd, i_Gsd, i_noise_pct)
man1=make_measure(seq1, manual_anoise_sd, manual_inoise_sd)
aut1=make_measure(seq1, aut_anoise_sd, aut_inoise_sd)
#segmentation errors
listN1=make_align_list(N1)
AUT=segmentation_errors(aut1, listN1,
                        organ_loss=c(6,6,12),
                        organ_gain=c(1,1,15,25,25))

AUT=segmentation_errors(aut1, listN1,
                        organ_loss=c(6,6,12),
                        organ_gain=c(1,1,15,25,25))

AUT=segmentation_errors(aut1, listN1,
                        organ_loss=c(6,6,12),
                        organ_gain=c(15))

multiseq_plot(list(seq1, AUT$values), align.df = AUT$I, ref.first = TRUE)

multiseq_plot(list(man1, aut1), id.names = c("manual", "before_seg_errors"), ref.first = TRUE)
multiseq_plot(list(man1, AUT$values), id.names = c("manual", "mypipe"), align.df = AUT$I, ref.first = TRUE)

###################
#Improve display of multisep plot / title separated from plot1 and common for the entire figure
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

g_title<-function(a.gplot){
    tmp <- ggplot_gtable(ggplot_build(a.gplot))
    #title_grob_idx <- which(sapply(tmp$grobs, function(x) x$name) == "title")
    title_grob_idx <- which(tmp$layout$name == "title")
    title_grob <- tmp$grobs[[title_grob_idx]]
    return(title_grob)}

myplot=ggplot(seq, aes(x=intervals, y=angles, color=factor))+geom_point()
myplot.title=myplot+ggtitle(label="this is a plot")

tmp. <- ggplot_gtable(ggplot_build(myplot))
which(sapply(tmp.$grobs, function(x) x$name) == "title")
leg. <- which(sapply(tmp.$grobs, function(x) x$name) == "guide-box")
sapply(tmp.$grobs, function(x) x$name)
tmp. <- ggplot_gtable(ggplot_build(myplot.title))
which(sapply(tmp.$grobs, function(x) x$name) == "title")

mylegend=g_legend(myplot)

mytitle_grob=g_title(myplot)
mytitle_grob=g_title(myplot+ggtitle(label="this is a plot"))

grid.arrange(mytitle_grob,
             arrangeGrob(myplot + theme(legend.position="none"),
                         myplot + theme(legend.position="none"),
                         nrow=2),
             mylegend, nrow=3, heights=c(1, 8, 2))

multiseq_plot(list(seq1, AUT$values), align.df = AUT$I, title = "test")

###################
#Correct display bug in multiseq_plot for trailing_tails (2021-01-04)
N1=14
alpha=137.5
a_sd=18.5
i_Gsd=0.8
i_noise_pct=75
seq1=make_refseq(N1, alpha, a_sd, i_Gsd, i_noise_pct)
listN1=make_align_list(N1)
seq2=segmentation_errors(seq1,listN1,
                         organ_gain=c(1,1,16,16,16))
multiseq_plot(list(seq1, seq2$values), align.df = seq2$I, title="tails at both ends")

seq3=segmentation_errors(seq1,listN1,
                         organ_gain=c(1,1))
multiseq_plot(list(seq1, seq3$values), align.df = seq3$I, title="starting tail only")

seq4=segmentation_errors(seq1,listN1,
                         organ_gain=c(16,16))
multiseq_plot(list(seq1, seq4$values), align.df = seq4$I, title="trailing tail only")

seq5=segmentation_errors(seq1,listN1,
                         organ_gain=c(1,1),organ_loss = c(14,15))
multiseq_plot(list(seq1, seq5$values), align.df = seq5$I, title="starting tail + trailing chop")

seq6=segmentation_errors(seq1,listN1,
                         organ_loss = c(1,2,14,15))
multiseq_plot(list(seq1, seq6$values), align.df = seq6$I, title="chops at both ends")

seq7=segmentation_errors(seq1,listN1,
                         organ_gain=c(16,16), organ_loss = c(1,2))
multiseq_plot(list(seq1, seq7$values), align.df = seq7$I, title="starting chops + trailing tail")

seq8=segmentation_errors(seq1,listN1,
                         organ_gain=c(2,3,13,14), organ_loss = c(1,2,13,14))
multiseq_plot(list(seq1, seq8$values), align.df = seq8$I, title="complex events - no apparent effects on interval number")
multiseq_plot(list(seq1, seq8$values), title="starting chops + trailing tail")

###################
#Testing permutations (2021-03-02)
####
N1=20
alpha=137.5
a_sd=18.5
i_Gsd=0.8
i_noise_pct=75
seq1=make_refseq(N1, alpha, a_sd, i_Gsd, i_noise_pct)
listN1=make_align_list(N1)
res2=segmentation_errors(seq1,listN1,
                         organ_gain=c(1,1,16,16),
                         organ_loss = c(2,5,13))
#Note for later: better choose 'res' than 'seq' for list of $values, $I, $O...
multiseq_plot(list(seq1, res2$values), align.df = res2$I, title="test")

#Debug/update first 'check_align'
test1=list(res2$I, res2$O)
test2=list(res2$O, res2$I)
test3=test2
colnames(test3[[2]])=c("a", "b", "c")
check_align_list(res2, verbose = TRUE) #bad input
check_align_list(test1, verbose = TRUE) #input in the right order
check_align_list(test2, verbose = TRUE) #input in the reverse order, still works
check_align_list(test3, verbose = TRUE)

#Debug/Update 'dtw_code_compute'
a=c(6,2,9,3,14,1)

#Function seq_simple_permut
seq_simple_permut(test1, 22, verbose = TRUE) #error: out of range
seq_simple_permut(test1, 0,verbose = TRUE)  #error: out of range
a=seq_simple_permut(test1, 1, verbose = TRUE) #No modification because pre-existing segmentation errors
a=seq_simple_permut(test1, 3, verbose = TRUE) #No modification because pre-existing segmentation errors
a=seq_simple_permut(test1, 5, verbose = TRUE) #No modification because pre-existing segmentation errors
a=seq_simple_permut(test1, 12, verbose = TRUE) #No modification because pre-existing segmentation errors
a=seq_simple_permut(test1, 16, verbose = TRUE) #No modification because pre-existing segmentation errors
a=seq_simple_permut(test1, 4,verbose = TRUE) #Real example, it should work for both $I and $O

sum(as.character(a$Ialign$dtw) == as.character(test1[[1]]$dtw))
sum(a$Oalign$modified == test1[[2]]$modified)
  
#Function 'value_simple_permut'
#Only permutations (without segmentation errors)
seq1_perm=value_simple_permut(seq1, idx=4)
seq1_perm=value_simple_permut(seq1_perm, idx=11)
seq1_perm=value_simple_permut(seq1_perm, idx=14)
multiseq_plot(list(seq1, seq1_perm), title="just permutation(s)")
#testing border/limits
seq1_perm=value_simple_permut(seq1_perm, idx=1)
seq1_perm=value_simple_permut(seq1_perm, idx=20)
multiseq_plot(list(seq1, seq1_perm), title="just permutation(s)")
#With Segmentation errors
multiseq_plot(list(seq1, res2$values), align.df = res2$I, title="seg errors + Permutation") #before permutation
#the interval n°4 can be permuted (seg errors on both sides)
seq2_perm=value_simple_permut(res2$values, idx=4)
multiseq_plot(list(seq1, seq2_perm), align.df = res2$I, title="seg errors + Permutation")

#Testing wrapper for permutation
#seq1 (random), no segmentation errors
res=simple_measure_permutation(seq1, listN1, i_threshold=2, proba=1, verbose=TRUE)
seq_plot(res$values)
multiseq_plot(list(seq1, res$values), align.df = res$I, title="multiple measure permutations")
#seq1, no short internodes
seq_test=seq1
seq_test$internodes[which(seq_test$internodes <=2)]=10
res=simple_measure_permutation(seq_test, listN1, i_threshold=2, proba=1, verbose=TRUE)
multiseq_plot(list(seq_test, res$values), align.df = res$I, title="no permutations")
#Add more values below threshold
seq_test=seq1
seq_test$internodes[15]=0
res=simple_measure_permutation(seq_test, listN1, i_threshold=2, proba=1, verbose=TRUE)
multiseq_plot(list(seq_test, res$values), align.df = res$I, title="multiple measure permutations")
#testing cat permut and drift
seq_test$internodes[c(14:16)]=1
res=simple_measure_permutation(seq_test, listN1, i_threshold=2, proba=1, verbose=TRUE)
multiseq_plot(list(seq_test, res$values), align.df = res$I, title="multiple measure permutations")
#Test borders
seq_test=seq1
#Start
seq_test$internodes[1]=0
res=simple_measure_permutation(seq_test, listN1, i_threshold=2, proba=1, verbose=TRUE)
#End
seq_test$internodes[N1]=0
res=simple_measure_permutation(seq_test, listN1, i_threshold=2, proba=1, verbose=TRUE)

#Force multiple 0 values to test the permutations
n=4 #nber of imposed null internodes (there may be more if some internodes are already null)
null_internodes=sort(round(runif(4,min=1, max=nrow(seq1)), digits=0))
null_internodes=null_internodes[!duplicated(null_internodes)]
seq_test=seq1
seq_test$internodes[null_internodes]=0
seq_plot(seq_test)
res=simple_measure_permutation(seq_test, listN1, i_threshold=2, proba=1, verbose=TRUE)
multiseq_plot(list(seq_test, res$values), align.df = res$I, title="multiple measure permutations")

#Test consecutive permutations
seq_test=seq1
seq_test$internodes[c(4:8)]=0
res=simple_measure_permutation(seq_test, listN1, i_threshold=2, proba=1, verbose=TRUE)
multiseq_plot(list(seq_test, res$values), align.df = res$I, title="multiple measure permutations")

#Test proba / threshold
res=simple_measure_permutation(seq1, listN1, i_threshold=0, proba=0.5, verbose=TRUE)
multiseq_plot(list(seq1, res$values), align.df = res$I, title="multiple measure permutations")

#Real case scenario
N1=20
alpha=137.5
a_sd=18.5
i_Gsd=0.8
i_noise_pct=75
N1.seq=make_refseq(N1, alpha, a_sd, i_Gsd, i_noise_pct)
N1.list=make_align_list(N1)
N1.err.seg=segmentation_errors(N1.seq,N1.list,
                            organ_gain=c(1,1,3,4),
                            organ_loss = c(2,5,13,18,20,21))
N1.err.seg.perm=simple_measure_permutation(N1.err.seg$values, list(N1.err.seg$I, N1.err.seg$O),
                                           i_threshold = 2, proba = 1, verbose = TRUE)
multiseq_plot(list(seq1, N1.err.seg$values), align.df = N1.err.seg$I, title="step1: segmentation errors")
multiseq_plot(list(seq1, N1.err.seg.perm$values), align.df = N1.err.seg.perm$I, title="step2: adding permutations")

##################
#Pb of overlapping tails & chops (April 2021)
#Note: test with tails & chops
#Real case scenario
N1=19 #so 20 organs
alpha=137.5
a_sd=18.5
i_Gsd=0.8
i_noise_pct=75
N1.seq=make_refseq(N1, alpha, a_sd, i_Gsd, i_noise_pct)
N1.list=make_align_list(N1)

#test1 with reproducible error @seq.append -> ok debug in arg passing for 'dtw_code_compute'
GAIN=c(2,10,15,21,21,21)
LOSS=c(1,8,18)
N1.err.seg=segmentation_errors(N1.seq,N1.list,
                               organ_gain=GAIN,
                               organ_loss = LOSS, verbose = TRUE)
multiseq_plot(list(N1.seq, N1.err.seg$values), align.df = N1.err.seg$I)

#test2 @start
GAIN=c(4,10,15,21,21,21)
LOSS=c(1,2,3,4,8,18)
N1.err.seg=segmentation_errors(N1.seq,N1.list,
                               organ_gain=GAIN,
                               organ_loss = LOSS, verbose = TRUE)
multiseq_plot(list(N1.seq, N1.err.seg$values), align.df = N1.err.seg$I)

#test2 @end
GAIN=c(21,21,21,21)
LOSS=c(19,20)
N1.err.seg=segmentation_errors(N1.seq,N1.list,
                               organ_gain=GAIN,
                               organ_loss = LOSS, verbose = TRUE)
multiseq_plot(list(N1.seq, N1.err.seg$values), align.df = N1.err.seg$I)

##################
#Pb of overlapping tails & chops (May 2021)
#Note: test with tails & chops
#Real case scenario
#Real case scenario
N1=19 #so 20 organs
alpha=137.5
a_sd=18.5
i_Gsd=0.8
i_noise_pct=75
N1.seq=make_refseq(N1, alpha, a_sd, i_Gsd, i_noise_pct)
N1.list=make_align_list(N1)

#test1: progressive
#chop the end
GAIN=NULL
LOSS=c(18,19,20)
N1.err.seg=segmentation_errors(N1.seq,N1.list,
                               organ_gain=GAIN,
                               organ_loss = LOSS, verbose = TRUE)
multiseq_plot(list(N1.seq, N1.err.seg$values), align.df = N1.err.seg$I)

#overlapping Merge and chop: 1 chop just before -> OK
GAIN=c(17)
LOSS=c(18,19,20)
N1.err.seg=segmentation_errors(N1.seq,N1.list,
                               organ_gain=GAIN,
                               organ_loss = LOSS, verbose = TRUE)
multiseq_plot(list(N1.seq, N1.err.seg$values), align.df = N1.err.seg$I)

#overlapping Merge and chop: 1 chop at the 1st chop -> OK
GAIN=c(18)
LOSS=c(18,19,20)
N1.err.seg=segmentation_errors(N1.seq,N1.list,
                               organ_gain=GAIN,
                               organ_loss = LOSS, verbose = TRUE)
multiseq_plot(list(N1.seq, N1.err.seg$values), align.df = N1.err.seg$I)

#overlapping Merge and chop: 1 chop at the 2nd chop -> OK
GAIN=c(19)
LOSS=c(18,19,20)
N1.err.seg=segmentation_errors(N1.seq,N1.list,
                               organ_gain=GAIN,
                               organ_loss = LOSS, verbose = TRUE)
multiseq_plot(list(N1.seq, N1.err.seg$values), align.df = N1.err.seg$I)

#overlapping Merge and chop: 1 chop at the 2nd chop -> OK
GAIN=c(20)
LOSS=c(18,19,20)
N1.err.seg=segmentation_errors(N1.seq,N1.list,
                               organ_gain=GAIN,
                               organ_loss = LOSS, verbose = TRUE)
multiseq_plot(list(N1.seq, N1.err.seg$values), align.df = N1.err.seg$I)

#overlapping Merge and chop: 1 chop at the 2nd chop -> OK
GAIN=c(21)
LOSS=c(18,19,20)
N1.err.seg=segmentation_errors(N1.seq,N1.list,
                               organ_gain=GAIN,
                               organ_loss = LOSS, verbose = TRUE)
multiseq_plot(list(N1.seq, N1.err.seg$values), align.df = N1.err.seg$I)


GAIN=c(15)
LOSS=c(15)
N1.err.seg=segmentation_errors(N1.seq,N1.list,
                               organ_gain=GAIN,
                               organ_loss = LOSS, verbose = TRUE)
multiseq_plot(list(N1.seq, N1.err.seg$values), align.df = N1.err.seg$I)

###################
## Specific issue of overlapping Chops and and Tails at the end
## October 2021
################
N=19 #so 20 organs
alpha=137.5
a_sd=18.5
i_Gsd=0.8
i_noise_pct=75
init.seq=make_refseq(N, alpha, a_sd, i_Gsd, i_noise_pct)
init.align=make_align_list(N)

#Visualize the problem: starting tails are separated but trailing tails are fused and flag as C/T events
#####
GAIN=c(1,1,1, 21,21,21)
LOSS=c(1,2,3, 18,19,20)
test.seg=segmentation_errors(in.seq = init.seq, align.list = init.align,
                               organ_gain=GAIN,
                               organ_loss = LOSS, verbose = TRUE)

multiseq_plot(mylist=list(init.seq, test.seg$values), align.df = test.seg$I, 
              id.names=c("init", "test"), verbose=TRUE)

##Observe the print:
#Just before adding the Tails:
# [1] "step n° 9 ; type= loss idx =  18"
# ref.idx modif.idx type
# 1        1        -1 gain
# 2        1        -1 gain
# 3        1        -1 gain
# 4        1        -1 loss
# 5        2         0 loss
# 6        3         1 loss
# 7       18        16 loss
# 8       19        17 loss
# 9       20        18 loss
# 10      21        19 gain
# 11      21        19 gain
# 12      21        19 gain
# 
# $Ialign
# reference modified dtw
# 1         NA        1   T
# 2         NA        2   T
# 3         NA        3  TS
# 4          1        3  SS
# 5          2        3  SS
# 6          3        3   S
# 7          4        4   ~
#   8          5        5   ~
#   9          6        6   ~
#   10         7        7   ~
#   11         8        8   ~
#   12         9        9   ~
#   13        10       10   ~
#   14        11       11   ~
#   15        12       12   ~
#   16        13       13   ~
#   17        14       14   ~
#   18        15       15   ~
#   19        16       16   ~
#   20        17       NA   C
# 21        18       NA   C
# 22        19       NA   C
# 
# $Oalign
# reference modified segmentation
# 1         NA        1         tail
# 2         NA        2         tail
# 3         NA        3         tail
# 4          1       NA        under
# 5          2       NA        under
# 6          3       NA        under
# 7          4        4            ~
#   8          5        5            ~
#   9          6        6            ~
#   10         7        7            ~
#   11         8        8            ~
#   12         9        9            ~
#   13        10       10            ~
#   14        11       11            ~
#   15        12       12            ~
#   16        13       13            ~
#   17        14       14            ~
#   18        15       15            ~
#   19        16       16            ~
#   20        17       17            ~
#   21        18       NA        under
# 22        19       NA        under
# 23        20       NA         chop

#Then adding the 1st tail
# [1] "step n° 10 ; type= gain idx =  18"
# ref.idx modif.idx type
# 1        1        -2 gain
# 2        1        -2 gain
# 3        1        -2 gain
# 4        1        -2 loss
# 5        2        -1 loss
# 6        3         0 loss
# 7       18        15 loss
# 8       19        16 loss
# 9       20        17 loss
# 10      21        18 gain
# 11      21        18 gain
# 12      21        18 gain

# Reference sequence contains 20 organs.
# Currently modified sequence contains 17 organs.
# Appending organs at a sequence tail -> use 'seq_append' function instead. 
# Reference sequence contains 20 organs.
# Currently modified sequence contains 17 organs.
# After edition, the modified sequence now contains 18 organs.
# This simply adds values at the end of sequences (not a proper merge case). The function 'value_tails' was used instead.

# $Ialign
# reference modified dtw
# 1         NA        1   T
# 2         NA        2   T
# 3         NA        3  TS
# 4          1        3  SS
# 5          2        3  SS
# 6          3        3   S
# 7          4        4   ~
#   8          5        5   ~
#   9          6        6   ~
#   10         7        7   ~
#   11         8        8   ~
#   12         9        9   ~
#   13        10       10   ~
#   14        11       11   ~
#   15        12       12   ~
#   16        13       13   ~
#   17        14       14   ~
#   18        15       15   ~
#   19        16       16   ~
#   20        17       17  CT
# 21        18       NA   C
# 22        19       NA   C

# $Oalign
# reference modified segmentation
# 1         NA        1         tail
# 2         NA        2         tail
# 3         NA        3         tail
# 4          1       NA        under
# 5          2       NA        under
# 6          3       NA        under
# 7          4        4            ~
#   8          5        5            ~
#   9          6        6            ~
#   10         7        7            ~
#   11         8        8            ~
#   12         9        9            ~
#   13        10       10            ~
#   14        11       11            ~
#   15        12       12            ~
#   16        13       13            ~
#   17        14       14            ~
#   18        15       15            ~
#   19        16       16            ~
#   20        17       17            ~
#   21        18       18    chop/tail
# 22        19       NA        under
# 23        20       NA         chop
#####

#Modif the `seq_append` function, only work at the end of the sequences
#####
#1. create the ending chops
GAIN=NULL
LOSS=c(18,19,20)
test.seg=segmentation_errors(in.seq = init.seq, align.list = init.align,
                             organ_gain=GAIN,
                             organ_loss = LOSS, verbose = TRUE)
multiseq_plot(mylist=list(init.seq, test.seg$values), align.df = test.seg$I, 
              id.names=c("init", "test"), verbose=TRUE)

#2. add the trailing tails
test.list=list(test.seg$I, test.seg$O)
names(test.list)=c("I", "O")

#list.plus.tail=
seq_append(test.list, 18, verbose=TRUE)
#OK, expected result

#3. the other function for values (`value_tails`) does not need to be changed ? -> test it

#4. Test again the previous
GAIN=c(21, 21,21,21) #testing with 21, then c(21, 21), then c(21, 21, 21), etc...
LOSS=c(18,19,20)
test.seg=segmentation_errors(in.seq = init.seq, align.list = init.align,
                             organ_gain=GAIN,
                             organ_loss = LOSS, verbose = TRUE)
multiseq_plot(mylist=list(init.seq, test.seg$values), align.df = test.seg$I, 
              id.names=c("init", "test"), verbose=TRUE)
#OK, that's the right behavior

#5. further tests
GAIN=21
LOSS=c(17,18)
#all the tests below run ok and plots give the expected results
#        test1        test2     test3     test4       test5     test6
#GAIN  c(21,21,21,21)  20         21        21          21        21
#LOSS     NULL         20         20      c(19,20)   c(18,19)   c(17,18)
test.seg=segmentation_errors(in.seq = init.seq, align.list = init.align,
                             organ_gain=GAIN,
                             organ_loss = LOSS, verbose = TRUE)
multiseq_plot(mylist=list(init.seq, test.seg$values), align.df = test.seg$I, 
              id.names=c("init", "test"), verbose=TRUE)

###################
## Pb to initialise testnoise
## November 2021
################
N=19 #so 20 organs
alpha=137.5
a_sd=18.5
i_Gsd=0.8
i_noise_pct=75
init.seq=make_refseq(N, alpha, a_sd, i_Gsd, i_noise_pct)
init.align=make_align_list(N)

#####
GAIN=NULL
LOSS=NULL
test.seg=segmentation_errors(in.seq = init.seq, align.list = init.align,
                             organ_gain=GAIN,
                             organ_loss = LOSS, verbose = TRUE)

multiseq_plot(mylist=list(init.seq, test.seg$values), align.df = test.seg$I, 
              id.names=c("init", "test"), verbose=TRUE)

###################
## Change the way noise is encoded
## November 2021
################
N=19 #so 20 organs
alpha=137.5
a_sd=35
i_Gsd=0.8
i_noise_pct=75
init.seq=make_refseq(N, alpha, a_sd, i_Gsd, i_noise_pct)
init.align=make_align_list(N)
seq_plot(init.seq)

#Make a measure with noise
anoise=0.5
inoise=0.5
noise.seq=make_measure(init.seq, anoise_sd=anoise, inoise_sd = inoise, 
                       noise.scale = "sd", 
                       anoise.mean=0, inoise.mean=0, verbose = TRUE)
multiseq_plot(mylist=list(init.seq, noise.seq), id.names = c("ref", "noise"))

###################
## Implement natural permutations
## December 2021
################
N=19 #so 20 organs
alpha=137.5
a_sd=35
i_Gsd=0.8
i_noise_pct=75

test=c(0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0)
get_consecutive_idx(test, value=0)

init.seq=make_refseq(N, alpha, a_sd, i_Gsd, i_noise_pct, 
                     natural.permutation = TRUE,
                     verbose=TRUE)
seq_plot(init.seq)

init.seq=make_refseq(N, alpha, a_sd, i_Gsd, i_noise_pct, 
                     natural.permutation = TRUE, permutation.frequency=1,
                     verbose=TRUE)
seq_plot(init.seq)

###################
## Improving internodes realism
## December 2021
################
N=25 #so 26 organs
alpha=137.5
a_sd=22
permutation.frequency=1
i_Gsd=5
i_noise_pct=50

init.seq=make_refseq(N, alpha, a_sd, 
                     natural.permutation = TRUE, permutation.frequency=permutation.frequency,
                     i_Gsd, i_noise_pct,
                     i_beta=1, i_max=80, i_plateau=5,
                     verbose=FALSE)
seq_plot(init.seq)

###################
## Improving the measures of internodes (ob of null internodes)
## Jan 2022
################
N=25 #so 20 organs
alpha=137.5
a_sd=30
permutation.frequency=0.04
i_Gsd=1.5
i_beta=1.5
i_noise_pct=75
i_max=100
i_plateau=5

init.seq=make_refseq(N, alpha, a_sd,
                     natural.permutation = TRUE, permutation.frequency=permutation.frequency,
                     i_Gsd=i_Gsd, i_noise_pct=i_noise_pct,
                     i_beta=i_beta, i_max=i_max, i_plateau=i_plateau,
                     verbose=FALSE)
seq_plot(init.seq)

#Make a measure with noise
anoise=50
inoise=5
noise.seq=make_measure(init.seq, anoise_sd=anoise, inoise_sd = inoise, 
                       noise.scale = "absolute", 
                       anoise.mean=0, inoise.mean=0, verbose = TRUE)
multiseq_plot(mylist=list(init.seq, noise.seq), id.names = c("ref", "noise"))

