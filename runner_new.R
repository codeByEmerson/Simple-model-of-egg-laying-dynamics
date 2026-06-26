# Runs the model specified in model.R with various parameter choices
#   Note: Updated for Ahl = Ad in expDecay 

analyzeRun <- function(runData) {
    # generates a number of potentially informative plots
    pdf(width=4,height=4)
    plot(x=1:STOPTIME,y=-runData$Idat, type="l", col="red", 
        ylim=c(-max(runData$Idat)-1,max(runData$Adat)+1),lwd=3)
    points(x=1:STOPTIME,y=runData$Adat, type="l", col="green")
    points(x=1:STOPTIME,y=runData$Adat-runData$Idat, type="l", col="black")

    tp=120
    trimx<-1:min(STOPTIME,tp)
    trimy<-function(d) {return(d[1:min(STOPTIME,tp)])}

    plot(x=trimx,y=trimy(-runData$Idat), type="l", col="red", 
        ylim=c(-max(runData$Idat)-1,max(runData$Adat)+1),lwd=3)
    points(x=trimx,y=trimy(runData$Adat), type="l", col="green")
    points(x=trimx,y=trimy(runData$Adat-runData$Idat), type="l", col="black")

    tp=120
    trimx<-tp:min(STOPTIME,tp*10)
    trimy<-function(d) {return(d[tp:min(STOPTIME,tp*10)])}

    plot(x=trimx,y=trimy(-runData$Idat), type="l", col="red", 
        ylim=c(-max(runData$Idat)-1,max(runData$Adat)+1),lwd=3)
    points(x=trimx,y=trimy(runData$Adat), type="l", col="green")
    points(x=trimx,y=trimy(runData$Adat-runData$Idat), type="l", col="black")

    

    # we preassign the cycles vector to all zeroes to avoid accumulation slowdown
    #   in exchange, we have to toss all those zeroes out (including potentially
    #   legitimate ones, r.i.p.)
    cycleLens<-runData$cycles[runData$cycles!=0]

    print("CYCLES")
    print(length(cycleLens))

    print("ACTUAL MEAN AND VAR")
    print(mean(cycleLens))
    print(var(cycleLens))


    hist(cycleLens,breaks=30)
    plot(density(cycleLens))

    hist(log(cycleLens),breaks=30)
    plot(density(log(cycleLens)))

    # Predict the model outcomes for the specified parameters
    # Predictions are only accurate for the simple model !!!

    #mean <- (1)/Ap + ( (Ihl * log(As/Is) ) / log(1/2) )
    mean <- (1)/runData$param_Ap + ( (runData$param_Ihl * log(runData$param_As/(runData$param_Is+runData$param_As)) ) / log(1/2) ) - 1

    varX <- (1-runData$param_Ap)/(runData$param_Ap**2)

    print("PREDICTED MEAN AND VAR")
    print(mean)
    print(varX)
    
    source("stochmath.R")
    #Ap, As, Ad, Is, Ihl -> expected mean in the stochastic model
    # assumes Is_actual = (lambda(1)+1)*Is
    E_mean <- expmean_stoch(runData$param_Ap, runData$param_As, runData$param_Ad, 
                            runData$param_Is, runData$param_Ihl)
    print("PREDICTED MEAN STOCHMODEL")
    print(E_mean)

    # if we want the actual TIMES of the cycles, we can get this from the position
    #   of the cycleLengths in state$cycles
    # OR, potentially faster, by summing/replacing on the zeroless list
    q<-Reduce(sum, cycleLens, 0, accumulate=TRUE)
    q<-q[-1] # remove the initial zero added by 'reduce'ing
    # note that the values in q will be relative to TIME, not POSITION 
    #   (0-indexed, not 1-indexed)
    #   e.g: cycle @ position 194 occured at time 193, and will be reported as-such
    #print(q)

    source("slicerTool.R")
    
    # !!!WARNING WARNING WARNING TAKE HEED TAKE HEED TAKE HEED!!!
    # !!!WARNING WARNING WARNING TAKE HEED TAKE HEED TAKE HEED!!!
    # !!!WARNING WARNING WARNING TAKE HEED TAKE HEED TAKE HEED!!!
    # !!! YOU MUST SET THE tmin TO THE GENOTYPE-SPECIFIC VALUE !!!
    # tmin is calculated from actual experimental data per-genotype, represents 
    #   the cutoff between active and inactive phase egg-to-egg intervals
    
    # all tmin are given in SECONDS
    # make sure you divide the tmin by 60 to translate it into MINUTES
    
    # tmin for wt is 221
    # tmin for nlp3 is 493
    # tmin for flp11nlp7 is 200
    # tmin for tdc1 is 239
    # tmin for tph1 is 216
    
    # currently this is set for wt data: tmin=221/60
    slicer(q, STOPTIME, tMin=221/60, unvarMean=(mean-(1/runData$param_Ap)) )
}

source("model.R")

# how many minutes to simulate. 300k is the standard we adopted, although>100k
#   seems to yield ~stable results. smaller runs give more variable results
#   b/c they undersample the distribution of possible outcomes
STOPTIME=300000

###################
# TO RUN THE MODEL#
###################
# 1) Set tMin= your genotype tMin (listed above)
# 2) Set STOPTIME to a value of your choice (300k recommended)
# 3) Uncomment the output<-... and analyzeRun... lines for your desired genotype (only 1 genotype can be run at a time!!!)
# 4) Make sure STOPTIME in the output<- line matches the value declared above.
# 5) Run this file in-place as Rscript runner_new.R (it needs to share a directory
#       with model.R, stochmath.R, and slicerTool.R)

#================================================
#=============BASIC WT BELOW=====================
#================================================

#output<-runModel(STOPTIME=300000, Is=5.55, As=1, Ihl=3, Ad=1, Ap=0.0764, capA=TRUE, randIs=FALSE, randAs=FALSE)
#analyzeRun(output)

#================================================
#=============COMPLEX WT BELOW===================
#================================================

output<-runModel(STOPTIME=300000, Is=7, As=2, Ihl=4, Ad=1.75, Ap=0.125, capA=FALSE, randIs=TRUE, randAs=TRUE, expDecay=TRUE)
analyzeRun(output)

#===============================================
#==============FLP11NLP7 BELOW==================
#===============================================

#output<-runModel(STOPTIME=300000, Is=6.5, As=2, Ihl=3.5, Ad=1.75, Ap=0.125, capA=FALSE, randIs=TRUE, randAs=TRUE, expDecay=TRUE)
#analyzeRun(output)

#================================================
#=============TDC1 MUTANT BELOW==================
#================================================

#output<-runModel(STOPTIME=300000, Is=19, As=2, Ihl=10.5, Ad=1.75, Ap=0.20, capA=FALSE, randIs=TRUE, randAs=TRUE, expDecay=TRUE)
#analyzeRun(output)


#================================================
#=============TPH1 MUTANT BELOW==================
#================================================

#output<-runModel(STOPTIME=300000, Is=7, As=1.2, Ihl=4, Ad=4.5, Ap=0.0525, capA=FALSE, randIs=TRUE, randAs=TRUE, expDecay=TRUE, randN=6)
#analyzeRun(output)

#================================================
#=============NLP3 MUTANT BELOW==================
#================================================

#output<-runModel(STOPTIME=300000, Is=7, As=0.4, Ihl=4, Ad=1, Ap=0.06, capA=FALSE, randIs=TRUE, randAs=TRUE, expDecay=TRUE, randN=1)
#analyzeRun(output)

