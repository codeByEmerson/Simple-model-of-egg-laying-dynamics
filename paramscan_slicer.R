# same code as slicerTool.R, just slightly modified to work with outputs from gridscanner

# slices up a time series of events to check for correlation between timepoints
# tMin = minimum time between events to include, use this for excluding active phases from consideration
paramscan_slicer <- function(eventTimes, endTime, tMin, unvarMean) {

    # have: list of eventTimes
    # want: list of event interval lengths
    eventLengths<-diff(eventTimes)
    eventLengthstrim<-eventLengths[eventLengths>tMin]
    
    print(paste("mean event lengths (nomaxcap), tmin=",tMin))
    print(mean(eventLengthstrim))
    meanEventLength<-mean(eventLengthstrim)
    varEventLength<-var(eventLengthstrim)
    
    
    # slicing time
    # slice at every possible timepoint that falls between two egg-laying events
    #   use this to calculate expected waiting time between events based on
    #   how long we have already waited. this is an extension of the simple case
    #   outlined below.
    #
    # * the simple case:
    #   pick a handful of points during a waiting period between two events
    #   if the two events are correlated, then the time waited (w) will correspond to
    #   the time until the next event (u) 
    #   e.g: if 1:1 match with a mean wait time of 20 minutes, then 20 + (w*-1) = u
    #   but if the two are uncorrelated, then 20 + (w*0) = u, aka expected time is 
    #   always just the mean.

    slices=1:(max(eventTimes)-1)

    i=1
    step=1
    tlast=0
    slicePairs<- vector("list",length(slices)-1)
    while(i<length(slices)){
        t=eventTimes[step]
        # slice between 2 time points
        if(tlast<=slices[i] && t>=slices[i]) {
            # exclude slices that fall within an active phase
            if(t-tlast>tMin){
                slicePairs[[ i ]]<-c(slices[i]-tlast,t-slices[i])
            } else {
                slicePairs [[ i ]] <-c(0,0)
            }
            i=i+1
        } else {
            tlast=t
            step=step+1
        }
    }
    
    x<-sapply(slicePairs, function(pair) {pair[1]} )
    y<-sapply(slicePairs, function(pair) {pair[2]} )

    # note that for actual usable data-sizes (millions of points) this will fail
    #   to give any insight into what is happening
    lmodel <- lm(y ~ x)
    
    #   when handling experimental data, I set mins=floor(x/60)   
    #     this bins thousands of second-frames into tens of minute-frames
    #     making the data easier to see/plot/handle
    #   BUT for the model this doesn't make sense
    #     (model timestep is 1 minute not 1 second) so we skip this reduction 
    data <- data.frame(x = x, y = y, mins=floor(x))
    data <- data[order(data$x), ]

    bp <- vector("list",max(data$mins)+1)
    i=1
    while(i<length(data$mins)) {
        bp[[ data$mins[[i]]+1 ]] <- c(bp[[ data$mins[[i]]+1 ]],data$y[[i]])
        i=i+1
    }
    
    bp[sapply(bp, is.null)] <- NA

    bp.medians<-sapply( bp, function(x) {median(x)} )
    vect.medians<-unlist(bp.medians)

    # package the output for nice plotting
    finaldat <- data.frame(x=1:length(vect.medians), 
                           y=vect.medians,
                           weight=unlist( sapply(bp, function(x) {length(x)}) ) 
                         ) 
    
    # if tMin>0, the t=0 point will always be set to zero, weird thing, shows up in egg-data too
    #  so we just prune it and move on
    if (tMin>0) {
        finaldat <- data.frame( x= finaldat$x[-1], # ditch point at 0,0
                                y= finaldat$y[-1],
                                weight=finaldat$weight[-1] )
    }
    
    
    # fits a single line to the mean inactive phase vs time passed data
    #  used to assess accuracy vs experimental data
    trimX<-finaldat$x[1:30]
    trimY<-finaldat$y[1:30]
    calcModel30slope <- NA
    calcModel30intercept <- NA
    result = tryCatch({
        model30 <- lm(trimY ~ trimX)
        calcModel30slope <- coef(model30)[2]     # second param is slope
        calcModel30intercept <- coef(model30)[1] # first param is intercept
    }, warning = function(w) {
        warning(w)
    }, error = function(e) {
        warning(e)
    }, finally = {
        # nothing
    })
    
    return(list(histDat=finaldat,meanDat=meanEventLength,varDat=varEventLength,
                model30slope=calcModel30slope,
                model30intercept=calcModel30intercept))
    
    # don't save the outputs
    # if desired, you can write all the outputs to a csv for further processing
    #   this is used for other data for fancier plotting, but is not
    #   necessary here.
    #write.csv(finaldat, "testing.csv")
}
