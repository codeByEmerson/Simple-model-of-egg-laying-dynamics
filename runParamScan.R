# how many minutes to simulate. 300k is the standard we adopted, although>100k
#   seems to yield ~stable results. smaller runs give more variable results
#   b/c they undersample the distribution of possible outcomes
STOPTIME=25000

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

#use 25000

runAndCalcWT_mod2 <- function(STOPTIME=25000, uIs=7, uAs=2, uIhl=4, uAd=6, uAp=0.125, ucapA=FALSE, urandIs=TRUE, urandAs=TRUE, uexpDecay=TRUE, urandN=6, tMinSeconds=221) {

    real.wt<-read.csv("real_data/real_wt.csv")

    output<-runModel(STOPTIME, Is=uIs, As=uAs, Ihl=uIhl, 
                     Ad=uAd, Ap=uAp, capA=ucapA, randIs=urandIs, 
                     randAs=urandAs, expDecay=uexpDecay, randN=urandN)
    # tmin for wt is 221
    # tmin for nlp3 is 493
    # tmin for flp11nlp7 is 200
    # tmin for tdc1 is 239
    # tmin for tph1 is 216
    runDat<-analyzeRun(output, tMinSeconds)
    
    # for final output, we turn the list runDat into the data.frame dfput
    #  we throw out all the long-run data (e.g: step-by-step info), since it
    #  won't be used again after this
    dfput <- data.frame(
                t=output[['t']],
                I=output[['I']],
                A=output[['A']],
                cycleLen=output[['cycleLen']],
                cycleCount=output[['cycleCount']],
                param_Is=output[['param_Is']], 
                param_As=output[['param_As']], 
                param_Ihl=output[['param_Ihl']],
                param_Ad=output[['param_Ad']],
                param_Ap=output[['param_Ap']],
                param_capA=output[['param_capA']],
                param_randIs=output[['param_randIs']],
                param_randAs=output[['param_randAs']]
                )
    
    dfput[['meanInactiveLength']]<-runDat[['meanDat']]
    dfput[['histError']]<-calcHistoryError(runDat[['histDat']],real.wt)
    # 20.4 is wt mean
    dfput[['meanError']]<-rmsd(runDat[['meanDat']], 20.4)
    dfput[['varError']]<-rmsd(runDat[['varDat']], 158)
    
    dfput[['model30SlopeError']]<-rmsd(runDat[['model30slope']], -0.31) 
    dfput[['model30InterceptError']]<-rmsd(runDat[['model30intercept']], 16.5) 
    # div by max b/c slopeError is ~-1 to 1 and intercept error is ~-40 to 40
    dfput[['model30NetError']] <- (dfput[['model30SlopeError']]/max(dfput[['model30SlopeError']])) + (dfput[['model30InterceptError']]/max(dfput[['model30InterceptError']]))

    return(dfput)
}


runAndCalcFLP11NLP7_mod2 <- function(STOPTIME=25000, uIs=7, uAs=2, uIhl=3, uAd=6, uAp=0.125, ucapA=FALSE, urandIs=TRUE, urandAs=TRUE, uexpDecay=TRUE, urandN=6, tMinSeconds=200) {

    real.wt<-read.csv("real_data/real_flp11nlp7.csv")

    output<-runModel(STOPTIME, Is=uIs, As=uAs, Ihl=uIhl, 
                     Ad=uAd, Ap=uAp, capA=ucapA, randIs=urandIs, 
                     randAs=urandAs, expDecay=uexpDecay, randN=urandN)
    # tmin for wt is 221
    # tmin for nlp3 is 493
    # tmin for flp11nlp7 is 200
    # tmin for tdc1 is 239
    # tmin for tph1 is 216
    runDat<-analyzeRun(output, tMinSeconds)
    
    # for final output, we turn the list runDat into the data.frame dfput
    #  we throw out all the long-run data (e.g: step-by-step info), since it
    #  won't be used again after this
    dfput <- data.frame(
                t=output[['t']],
                I=output[['I']],
                A=output[['A']],
                cycleLen=output[['cycleLen']],
                cycleCount=output[['cycleCount']],
                param_Is=output[['param_Is']], 
                param_As=output[['param_As']], 
                param_Ihl=output[['param_Ihl']],
                param_Ad=output[['param_Ad']],
                param_Ap=output[['param_Ap']],
                param_capA=output[['param_capA']],
                param_randIs=output[['param_randIs']],
                param_randAs=output[['param_randAs']]
                )
    
    dfput[['meanInactiveLength']]<-runDat[['meanDat']]
    dfput[['histError']]<-calcHistoryError(runDat[['histDat']],real.wt)
    # 20.4 is wt mean
    dfput[['meanError']]<-rmsd(runDat[['meanDat']], 18.95)
    dfput[['varError']]<-rmsd(runDat[['varDat']], 123)
    
    dfput[['model30SlopeError']]<-rmsd(runDat[['model30slope']], -0.3) 
    dfput[['model30InterceptError']]<-rmsd(runDat[['model30intercept']], 16)
    # div by max b/c slopeError is ~-1 to 1 and intercept error is ~-40 to 40
    dfput[['model30NetError']] <- (dfput[['model30SlopeError']]/max(dfput[['model30SlopeError']])) + (dfput[['model30InterceptError']]/max(dfput[['model30InterceptError']]))

    return(dfput)
}


#uIs=use_Is, uAs=use_As, uIhl=use_Ihl, uAd=use_Ad, uAp=use_Ap, ucapA=use_capA, urandIs=use_randIs, urandAs=use_randAs
runAndCalcTDC1_mod2 <- function(STOPTIME=25000, uIs=32, uAs=2, uIhl=6.5, uAd=6, uAp=0.125, ucapA=FALSE, urandIs=TRUE, urandAs=TRUE, uexpDecay=TRUE, urandN=6, tMinSeconds=239) {

    real.wt<-read.csv("real_data/real_tdc1.csv")

    #STOPTIME, Is=uIs, As=uAs, Ihl=uIhl, 
    #Ad=uAd, Ap=uAp, capA=ucapA, randIs=urandIs, 
    #randAs=urandAs, expDecay=uexpDecay, randN=urandN)
    output<-runModel(STOPTIME, Is=uIs, As=uAs, Ihl=uIhl, 
                     Ad=uAd, Ap=uAp, capA=ucapA, randIs=urandIs, 
                     randAs=urandAs, expDecay=uexpDecay, randN=urandN)
    # tmin for wt is 221
    # tmin for nlp3 is 493
    # tmin for flp11nlp7 is 200
    # tmin for tdc1 is 239
    # tmin for tph1 is 216
    runDat<-analyzeRun(output, tMinSeconds)
    
    # for final output, we turn the list runDat into the data.frame dfput
    #  we throw out all the long-run data (e.g: step-by-step info), since it
    #  won't be used again after this
    dfput <- data.frame(
                t=output[['t']],
                I=output[['I']],
                A=output[['A']],
                cycleLen=output[['cycleLen']],
                cycleCount=output[['cycleCount']],
                param_Is=output[['param_Is']], 
                param_As=output[['param_As']], 
                param_Ihl=output[['param_Ihl']],
                param_Ad=output[['param_Ad']],
                param_Ap=output[['param_Ap']],
                param_capA=output[['param_capA']],
                param_randIs=output[['param_randIs']],
                param_randAs=output[['param_randAs']]
                )
    
    dfput[['meanInactiveLength']]<-runDat[['meanDat']]
    dfput[['histError']]<-calcHistoryError(runDat[['histDat']],real.wt)
    # 20.4 is wt mean
    dfput[['meanError']]<-rmsd(runDat[['meanDat']], 31.6)
    dfput[['varError']]<-rmsd(runDat[['varDat']], 222)
    
    dfput[['model30SlopeError']]<-rmsd(runDat[['model30slope']], -0.75)
    dfput[['model30InterceptError']]<-rmsd(runDat[['model30intercept']], 31) 
    # div by max b/c slopeError is ~-1 to 1 and intercept error is ~-40 to 40
    dfput[['model30NetError']] <- (dfput[['model30SlopeError']]/max(dfput[['model30SlopeError']])) + (dfput[['model30InterceptError']]/max(dfput[['model30InterceptError']]))
    
    

    return(dfput)
}

#output<-runModel(STOPTIME=300000, Is=7, As=0.5, Ihl=4, Ad=24, Ap=0.06, capA=FALSE, randIs=TRUE, randAs=TRUE, expDecay=TRUE, randN=1)
runAndCalcNLP3_mod2 <- function(STOPTIME=25000, uIs=7, uAs=0.5, uIhl=4, uAd=24, uAp=0.06, ucapA=FALSE, urandIs=TRUE, urandAs=TRUE, uexpDecay=TRUE, urandN=1, tMinSeconds=493) {

    real.wt<-read.csv("real_data/real_nlp3.csv")

    output<-runModel(STOPTIME, Is=uIs, As=uAs, Ihl=uIhl, 
                     Ad=uAd, Ap=uAp, capA=ucapA, randIs=urandIs, 
                     randAs=urandAs, expDecay=uexpDecay, randN=urandN)
    runDat<-analyzeRun(output, tMinSeconds)
    
    # for final output, we turn the list runDat into the data.frame dfput
    #  we throw out all the long-run data (e.g: step-by-step info), since it
    #  won't be used again after this
    dfput <- data.frame(
                t=output[['t']],
                I=output[['I']],
                A=output[['A']],
                cycleLen=output[['cycleLen']],
                cycleCount=output[['cycleCount']],
                param_Is=output[['param_Is']], 
                param_As=output[['param_As']], 
                param_Ihl=output[['param_Ihl']],
                param_Ad=output[['param_Ad']],
                param_Ap=output[['param_Ap']],
                param_capA=output[['param_capA']],
                param_randIs=output[['param_randIs']],
                param_randAs=output[['param_randAs']]
                )
    
    dfput[['meanInactiveLength']]<-runDat[['meanDat']]
    dfput[['histError']]<-calcHistoryError(runDat[['histDat']],real.wt)
    # 20.4 is wt mean
    dfput[['meanError']]<-rmsd(runDat[['meanDat']], 35.6)
    dfput[['varError']]<-rmsd(runDat[['varDat']], 353)
    
    dfput[['model30SlopeError']]<-rmsd(runDat[['model30slope']], -0.80)
    dfput[['model30InterceptError']]<-rmsd(runDat[['model30intercept']], 31)
    # div by max b/c slopeError is ~-1 to 1 and intercept error is ~-40 to 40
    dfput[['model30NetError']] <- (dfput[['model30SlopeError']]/max(dfput[['model30SlopeError']])) + (dfput[['model30InterceptError']]/max(dfput[['model30InterceptError']])) 


    return(dfput)
}

runAndCalcTPH1_mod2 <- function(STOPTIME=25000, uIs=7, uAs=1.1, uIhl=4, uAd=.15, uAp=0.0525, ucapA=FALSE, urandIs=TRUE, urandAs=TRUE, uexpDecay=TRUE, urandN=3, tMinSeconds=216) {

    print(urandN)
    real.wt<-read.csv("real_data/real_tph1.csv")

    output<-runModel(STOPTIME, Is=uIs, As=uAs, Ihl=uIhl, 
                     Ad=uAd, Ap=uAp, capA=ucapA, randIs=urandIs, 
                     randAs=urandAs, expDecay=uexpDecay, randN=urandN)
    runDat<-analyzeRun(output, tMinSeconds)
    
    # for final output, we turn the list runDat into the data.frame dfput
    #  we throw out all the long-run data (e.g: step-by-step info), since it
    #  won't be used again after this
    dfput <- data.frame(
                t=output[['t']],
                I=output[['I']],
                A=output[['A']],
                cycleLen=output[['cycleLen']],
                cycleCount=output[['cycleCount']],
                param_Is=output[['param_Is']], 
                param_As=output[['param_As']], 
                param_Ihl=output[['param_Ihl']],
                param_Ad=output[['param_Ad']],
                param_Ap=output[['param_Ap']],
                param_capA=output[['param_capA']],
                param_randIs=output[['param_randIs']],
                param_randAs=output[['param_randAs']]
                )
    
    dfput[['meanInactiveLength']]<-runDat[['meanDat']]
    dfput[['histError']]<-calcHistoryError(runDat[['histDat']],real.wt)
    # 20.4 is wt mean
    dfput[['meanError']]<-rmsd(runDat[['meanDat']], 30.2)
    dfput[['varError']]<-rmsd(runDat[['varDat']], 371)
    
    dfput[['model30SlopeError']]<-rmsd(runDat[['model30slope']], -0.30) 
    dfput[['model30InterceptError']]<-rmsd(runDat[['model30intercept']], 26)
    # div by max b/c slopeError is ~-1 to 1 and intercept error is ~-40 to 40
    dfput[['model30NetError']] <- (dfput[['model30SlopeError']]/max(dfput[['model30SlopeError']])) + (dfput[['model30InterceptError']]/max(dfput[['model30InterceptError']]))

    return(dfput)
}


