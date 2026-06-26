analyzeRun <- function(runData, tMinSeconds) {
    # we preassign the cycles vector to all zeroes to avoid accumulation slowdown
    #   in exchange, we have to toss all those zeroes out (including potentially
    #   legitimate ones, r.i.p.)
    cycleLens<-runData$cycles[runData$cycles!=0]

    
    mean <- (1)/runData$param_Ap + ( (runData$param_Ihl * log(runData$param_As/(runData$param_Is+runData$param_As)) ) / log(1/2) ) - 1
    
    
    # if we want the actual TIMES of the cycles, we can get this from the position
    #   of the cycleLengths in state$cycles
    # OR, potentially faster, by summing/replacing on the zeroless list
    q<-Reduce(sum, cycleLens, 0, accumulate=TRUE)
    q<-q[-1] # remove the initial zero added by 'reduce'ing
    # note that the values in q will be relative to TIME, not POSITION 
    #   (0-indexed, not 1-indexed)
    #   e.g: cycle @ position 194 occured at time 193, and will be reported as-such
    #print(q)

    source("paramscan_slicer.R")
    
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
    # 
    return( paramscan_slicer(q, STOPTIME, tMin=tMinSeconds/60, unvarMean=(mean-(1/runData$param_Ap)) ) )
}

rmsd <- function(model, real) {
    return( sqrt( (model-real)**2 ) )
}

calcHistoryError<- function(runDat, standard) {
    rmsd<-0
    i=1
    while (i<30) {
        # using RMSD for error measurement
        # 
        # divide standard by 60 b/c the real data in the standards
        #   are given in seconds, and model is in minutes
        rmsd<-rmsd + rmsd(runDat$y[i],standard$y[i]/60)
        i<-i+1
    }
    return(rmsd)
}

# use this to compare double's semi-safely
isSame<- function(db, val, tol=1e-6) {
    return(abs(db-val)<=tol)
}
