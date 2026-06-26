
# slices up a time series of events to check for correlation between timepoints
# tMin = minimum time between events to include, use this for excluding active phases from consideration
slicer <- function(eventTimes, endTime, tMin, unvarMean) {

    eventLengths<-diff(eventTimes)
    
    plot(density(log(eventLengths)), main="density(log(eventLengths))")
    
    eventLengthstrim<-eventLengths[eventLengths>tMin]
    
    eventLengthscut.den<-density(eventLengthstrim)
    
    eventLengthscut<-eventLengthstrim[eventLengthstrim<69]
    
    pdf(paste("slicergenny.pdf"), height=5,width=5)
    
    eventLengthscut.hist <- hist(eventLengthscut, breaks=seq((tMin),70,5),
                    main="Wild type", 
                    xlab="Time between egg-laying events (minutes)",
                    xlim=c(0,70),
                    xaxt="n",
                    col=rgb(0.2,.2,1))
    axis(1, at = c(0, tMin, 10, 20, 30, 40, 50, 60, 70),
            labels=c(0, round(tMin,2), 10, 20, 30, 40, 50, 60, 70), cex.axis=.5)
    lines(eventLengthscut.den$x, max(eventLengthscut.hist$counts)*eventLengthscut.den$y/max(eventLengthscut.den$y), lwd=5)
    text(55,15, paste("mean=",round(mean(eventLengthstrim), 3)))
    text(55,-2, paste("sd=",round(sd(eventLengthstrim), 3)))


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
    print(summary(lmodel))

    plot(x,y, main="model",
        xlab="Time since last event",
        ylab="Time until next event")
    abline(a=coef(lmodel)[1], b=coef(lmodel)[2], col="red")  # plots the line, a=intercept, b=slope 

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
    
    ylimmaxval=30
    
    # this plot shows just how "off" the linear fit can be
    plot(finaldat$x,finaldat$y, main="model, medians only, badlinear",
        xlab="Time since last event",
        ylab="Time until next event")
    abline(a=coef(lmodel)[1], b=coef(lmodel)[2], col="red")  # plots the line, a=intercept, b=slope 

    plot(finaldat$x,finaldat$y, main="model, medians only, badlinearsquat",
        xlab="Time since last event",
        ylab="Time until next event",
        ylim=c(0,ylimmaxval),
        xlim=c(0,30))
    abline(a=coef(lmodel)[1], b=coef(lmodel)[2], col="red")  # plots the line, a=intercept, b=slope 


    trimX<-finaldat$x[1:ceiling(17)]
    trimY<-finaldat$y[1:ceiling(17)]

    medModel <- lm(trimY ~ trimX)
    print("MEDMODEL")
    print(summary(medModel))
    
    secondHalfX<-finaldat$x[ceiling(17):30]
    secondHalfY<-finaldat$y[ceiling(17):30]
    
    secondHalfModel <- lm(secondHalfY ~ secondHalfX)
    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    print("~~~~~~~~~SECOND HALF MODEL~~~~~~")
    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    print(summary(secondHalfModel))

    plot(finaldat$x,finaldat$y, main="model, medians only, linefit BREAK AT 17 MANUAL",
        xlab="Time since last event",
        ylab="Time until next event",
        ylim=c(0,ylimmaxval),
        xlim=c(0,30))
    abline(a=coef(medModel)[1], b=coef(medModel)[2], col="red")  # plots the line, a=intercept, b=slope 
    abline(a=coef(secondHalfModel)[1], b=coef(secondHalfModel)[2], col="blue")  # plots the line, a=intercept, b=slope 


    #unvarMean=30
    trimX<-finaldat$x[1:ceiling(unvarMean)]
    trimY<-finaldat$y[1:ceiling(unvarMean)]

    medModel <- lm(trimY ~ trimX)
    print("MEDMODEL")
    print(summary(medModel))
    
    secondHalfX<-finaldat$x[ceiling(unvarMean):30]
    secondHalfY<-finaldat$y[ceiling(unvarMean):30]
    
    secondHalfModel <- lm(secondHalfY ~ secondHalfX)
    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    print("~~~~~~~~~SECOND HALF MODEL~~~~~~")
    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    print(summary(secondHalfModel))

    plot(finaldat$x,finaldat$y, main="model, medians only, linefit unvarmean",
        xlab="Time since last event",
        ylab="Time until next event",
        ylim=c(0,ylimmaxval),
        xlim=c(0,30))
    abline(a=coef(medModel)[1], b=coef(medModel)[2], col="red")  # plots the line, a=intercept, b=slope 
    abline(a=coef(secondHalfModel)[1], b=coef(secondHalfModel)[2], col="blue")  # plots the line, a=intercept, b=slope 

    # slap the slope on there for reference
    text(20, 20, paste("slope:",coef(medModel)[2]), col="red")
    # slap the slope on there for reference
    text(20, 15, paste("slope:",coef(secondHalfModel)[2]), col="blue")
    
    plot(finaldat$x,finaldat$y, main="model, medians only",
        xlab="Time since last event",
        ylab="Time until next event",
        ylim=c(0,20),
        xlim=c(0,30))
    abline(a=coef(medModel)[1], b=coef(medModel)[2], col="red")  # plots the line, a=intercept, b=slope 
    abline(a=coef(secondHalfModel)[1], b=coef(secondHalfModel)[2], col="blue")  # plots the line, a=intercept, b=slope 

    # slap the slope on there for reference
    text(20, 20, paste("slope:",coef(medModel)[2]), col="red")
    # slap the slope on there for reference
    text(20, 15, paste("slope:",coef(secondHalfModel)[2]), col="blue")

    plot(finaldat$x,finaldat$y, main="model, medians only, first half only",
        xlab="Time since last event",
        ylab="Time until next event",
        ylim=c(0,50),
        xlim=c(0,50))
    abline(a=coef(medModel)[1], b=coef(medModel)[2], col="red")  # plots the line, a=intercept, b=slope 
    
    # slap the slope on there for reference
    text(20, 20, paste("slope:",coef(medModel)[2]), col="red")
    
    
    plot(finaldat$x,finaldat$y, main="model, medians only, first half only",
        xlab="Time since last event",
        ylab="Time until next event",
        ylim=c(0,30),
        xlim=c(0,30))
    abline(a=coef(medModel)[1], b=coef(medModel)[2], col="red")  # plots the line, a=intercept, b=slope 

    # slap the slope on there for reference
    text(20, 20, paste("slope:",coef(medModel)[2]), col="red")
    
    unvarMean=30
    trimX<-finaldat$x[1:ceiling(unvarMean)]
    trimY<-finaldat$y[1:ceiling(unvarMean)]
    
    medModel <- lm(trimY ~ trimX)
    print("MEDMODEL-unvarmean=30")
    print(summary(medModel))

    plot(finaldat$x,finaldat$y, main="model, medians only, linefit unvarmean=30",
        xlab="Time since last event",
        ylab="Time until next event",
        ylim=c(0,20),
        xlim=c(0,30))
    abline(a=coef(medModel)[1], b=coef(medModel)[2], col="orange")  # plots the line, a=intercept, b=slope 

    # slap the slope on there for reference
    text(20, 20, paste("slope:",coef(medModel)[2]), col="orange")
    


    plot(finaldat$x,finaldat$y, main="model, medians only, linefit unvarmean=30",
        xlab="Time since last event",
        ylab="Time until next event",
        ylim=c(0,30),
        xlim=c(0,30))
    abline(a=coef(medModel)[1], b=coef(medModel)[2], col="orange")  # plots the line, a=intercept, b=slope 
    text(20, 20, paste("slope:",coef(medModel)[2]), col="orange")
    
    
    # don't save the outputs
    # if desired, you can write all the outputs to a csv for further processing
    #   this is used for the experimental data for fancier plotting, but is not
    #   necessary here.
    #write.csv(finaldat, "output.csv")
}
