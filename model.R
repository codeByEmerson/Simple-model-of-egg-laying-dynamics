# Simple model of egg-laying dynamics
#
# Goal: 
#   Demonstrate that an exponentially decaying inhibitory signal +
#     stochastic short-lived activating signal can combine to yield
#     similar on-off dynamics as exhibited by the actual egg-laying
#     circuit.

# half life is in minutes
decay <- function(conc, halflife, dt) {
    current <- conc * (1/2) ** (dt/halflife)
    return(current)
}

step = function(state, Is, As, Ihl, Ad, Ap, capA, randIs, randAs, expDecay, threshAdjust, randN) {
    if ((state$A-state$I) >= 0+threshAdjust) {
        # whether or not to randomly vary the amount of Is released
        if (randIs) {
            # rpois(1,lambda=1)+1 simulates egg-laying, min=1, max~=6
            # divide by 2 to recenter it properly w.r.t. Is
            # e.g: with no var: return Is*1s
            # with var: rpois(1,lambda=1)+1 mean=2
            #   so divide by 2 s.t. still 1-centered
            eggs<-rpois(1,lambda=1)+1
            state$I<-state$I+(Is*eggs/2)
        }
        else {
            state$I<-state$I+Is 
        }
        state$Istep[[state$t]] <- 1
        state$cycles[[state$t]]<-state$cycleLen
        state$cycleLen<-0
        state$cycleCount<-state$cycleCount+1
    }
    if (Ap >= runif(1)) {
        Aup <- As
        if (randAs) {
            Aup <- As*randN*runif(1)**(randN-1)
            if (randN==1) {
                Aup <- As*runif(1)
            }
        }
        
        # Note: capping A only seems to matter for large values of A, or when
        #   the decay rate of A is small (etc.)
        #   probably something like: As*Ap >= Ad
        # if we are capping the maximum value of A at As:
        if (capA) {
            if (expDecay) {
                state$A<-Aup 
                # boost it to bypass insta-decay
                state$A<-decay(state$A, Ad, -1) # this is -1 b/c we are increasing by 
                # decay step-size s.t. A+=As for this time step, then on the next step it will decay 
                #  e.g: A+=As + Ad - Ad
            }
            else { # linear decay
                state$A<-min(state$A+Aup+Ad,Aup+Ad) # eliminates A-stacking completely
            }
        }
        else { # else just pump it up
            if (expDecay) {
                state$A<-state$A+Aup 
                # boost it to bypass insta-decay
                state$A<-decay(state$A, Ad, -1) # 
            }
            else {
                state$A<-state$A+Aup+Ad
            }
        }
    }
    
    state$I<-decay(state$I, Ihl, 1)
    
    if (expDecay) {
        # use exponential decay, since this is complex model
        #  no longer using (As/Ad)/2 for halflife
        #  instead just use Ad as Ahl
        state$A<-decay(state$A, Ad, 1)
    }
    # use linear method, since this is simple model
    else { 
        state$A<-max(0,state$A-Ad)
    }
        
    state$Idat[[state$t]]<-state$I
    state$Adat[[state$t]]<-state$A
    
    state$cycleLen<-state$cycleLen+1
        
    state$t<-state$t+1
    return(state)
}

# runs the model until the stop-time is reached
# we do this instead of recursive calls to avoid stackoverflow loops
resolveModelTime <- function(state, Is, As, Ihl, Ad, Ap, capA, randIs, randAs, expDecay, threshAdjust, randN) {
    lastprintT<-0
    while (state$t <= STOPTIME) {
        if (state$t-lastprintT > STOPTIME/3) {
            #print(state$t)
            lastprintT<-state$t
        }
        state <- step(state, Is, As, Ihl, Ad, Ap, capA, randIs, randAs, expDecay, threshAdjust, randN)
    }
    return(state)
}

runModel <- function(STOPTIME, Is, As, Ihl, Ad, Ap, capA, randIs, randAs, expDecay, threshAdjust, randN) {
    if (missing(expDecay)) {
        expDecay <- FALSE # default value
    }
    
    if (missing(randN)) {
        randN <- 6 # default value
    }
    
    if (missing(threshAdjust)) {
        threshAdjust <- 0 # default value
    }

    # initialize the list that stores all model information as it runs
    #   t=current model running time
    #   I=current value of I at time==t
    #   A=current value of A at time==t
    #   Idat stores value of I at time==index
    #   Istep is a binary vector, 1=I was increased at time==index, else 0
    #   Adat stores value of A at time==index
    #   cycleLen=length of current cycle
    #   cycles is a vector, when a cycle ends, the cycleLen is stored at index==time, else 0
    #   cycleCount is the total number of cycles (equiv to length(cycles[cycles!=0]) )
    state=list(t=1,I=0,A=0,
                Idat=double(STOPTIME), # pre-allocate to avoid massive copy slowdown
                Istep=double(STOPTIME),
                Adat=double(STOPTIME),
                cycleLen=0,
                cycles=double(STOPTIME),
                cycleCount=0,
                param_Is=Is, # store the model parameters in the state so we don't
                param_As=As, #    have to keep passing them around for analysis
                param_Ihl=Ihl, #  later on
                param_Ad=Ad,
                param_Ap=Ap,
                param_capA=capA,
                param_randIs=randIs, # these are additions used by model 2 !
                param_randAs=randAs,
                param_expDecay=expDecay,
                param_threshAdjust=threshAdjust,
                param_randN=randN)
                
    output<-resolveModelTime(state, Is, As, Ihl, Ad, Ap, capA, randIs, randAs, expDecay, threshAdjust, randN)
    
    return(output)
}

