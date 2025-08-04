# Toy model of egg-laying dynamics
#
# Goal: 
#   Demonstrate that an exponentially decaying inhibitory signal +
#     stochastic short-lived activating signal can combine to yield
#     similar on-off dynamics as exhibited by the actual egg-laying
#     circuit.

# Interesting breakpoints
#   If Ad <= Ap*As, then A+I will accumulate over time
#   If Id is linear: Id/Is ~= rate of fire, !VERY SMALL! distortion from Ad
#       max distortion is ~= As/Id 
#           (e.g: +10 in a 1 per second regime = net 10-tick boost)
#           (this excludes double-fire leading to bigger boost rarely)
#   If Id is exp: ROF 100% depends on As, As~=10% converges towards linear


# dA/dt = At - 5 + (15 if rand(P(A)) else 0)
# dI/dt = decay(It, hl=33) + (100 if A-I>=0 else 0)
#
# As P(A) (probability of A being increased) goes -> 0, variance -> infinity
# As P(A) -> 1, variance -> 0
#
# If A is capped (only increase when =0), then maximum possible effect on linear model is:
#    As / (Is/Id) (step size in A / (step size in I / linear decay in I))
#    e.g: for As=10, Is=100, Id=1, max effect = 10%
#    
#



# half life is in seconds
decay <- function(conc, halflife, dt) {
    current <- conc * (1/2) ** (dt/halflife)
    return(current)
}

# removes 0's from a vector
killduds <- function(l) {
    return(unlist(lapply(l, function(x) {x[x!=0]})))
}


Is=100
As=5#10

Ihl=33

Id=1
Ad=5

Ap=.05 # odds that A 'fires'

# .01 197 vs 214 
# .01 222 vs 214
#.0011086 vs 1114

STOPTIME=60000
STOPCYCLES=300

step = function(state) {
    if ((state$A-state$I) >= 0) {
        state$I<-state$I+Is
        state$cycles[[state$t]]<-state$cycleLen
        state$cycleLen<-0
        state$cycleCount<-state$cycleCount+1
        #print("HELLO IT PINGED")
    }
    if (Ap >= runif(1)) {
        # remove the possibility of multiple-fires stacking, forces decay-intervals
        #if (state$A==0) {
            #state$A<-state$A+As+Ad # add decay-step extra to counteract immediate decay
        #}
        state$A<-min(state$A+As+Ad,As+Ad) # eliminates A-stacking completely
    }
    
    state$I<-decay(state$I, Ihl, 1)
    state$A<-max(0,state$A-Ad)
        
    state$Idat[[state$t]]<-state$I
    state$Adat[[state$t]]<-state$A
    
    state$cycleLen<-state$cycleLen+1
        
    state$t<-state$t+1
    return(state)
}

# runs the model until the stop-time is reached
# we do this instead of step calling itself to avoid dumb stackoverflow loops
resolveModelTime <- function(state) {
    lastprintT<-0
    while (state$t <= STOPTIME) {
        if (state$t-lastprintT > STOPTIME/20) {
            print(state$t)
            lastprintT<-state$t
        }
        state <- step(state)
    }
    return(state)
}

# runs the model until the stop-cycle count is reached
# we do this instead of step calling itself to avoid dumb stackoverflow loops
resolveModelCycles <- function(state) {
    lastprintCycle<-0
    while (state$cycleCount <= STOPCYCLES) {
        if (state$cycleCount-lastprintCycle > STOPCYCLES/20) {
            print(state$cycleCount)
            lastprintCycle<-state$cycleCount
        }
        state <- step(state)
    }
    return(state)
}


state=list(t=1,I=0,A=0,
            Idat=double(STOPTIME),
            Adat=double(STOPTIME),
            cycleLen=0,
            cycles=double(STOPTIME),
            cycleCount=0)

output<-resolveModelTime(state)

pdf(width=10,height=5)
plot(x=1:STOPTIME,y=-output$Idat, type="l", col="red", 
    ylim=c(-max(output$Idat)-10,max(output$Adat)+10),lwd=3)
points(x=1:STOPTIME,y=output$Adat, type="l", col="green")
points(x=1:STOPTIME,y=output$Adat-output$Idat, type="l", col="black")

trimx<-1:min(STOPTIME,2000)
trimy<-function(d) {return(d[1:min(STOPTIME,2000)])}

plot(x=trimx,y=trimy(-output$Idat), type="l", col="red", 
    ylim=c(-max(output$Idat)-10,max(output$Adat)+10),lwd=3)
points(x=trimx,y=trimy(output$Adat), type="l", col="green")
points(x=trimx,y=trimy(output$Adat-output$Idat), type="l", col="black")


# we preassign the cycles vector to all zeroes to avoid accumulation slowdown
#   in exchange, we have to toss all those zeroes out (including potentially
#   legitimate ones, r.i.p.)
cycleLens<-killduds(output$cycles)

print("ACTUAL MEAN AND VAR")
print(mean(cycleLens))
print(var(cycleLens))


hist(cycleLens,breaks=30)
plot(density(cycleLens))

hist(log(cycleLens),breaks=30)
plot(density(log(cycleLens)))

# prediction time

#mean <- (1)/Ap + ( (Ihl * log(As/Is) ) / log(1/2) )
mean <- (1)/Ap + ( (Ihl * log(As/(Is+As)) ) / log(1/2) )

varX <- (1-Ap)/(Ap**2)

print("PREDICTED MEAN AND VAR")
print(mean)
print(varX)


# predict the PMF
PMF<- function(x, p, const) {
    return(p * ((1-p)**(x-1)) )
}

const <- ( (Ihl * log(As/(Is+As)) ) / log(1/2) )

minPMF<-1
maxPMF<-max(cycleLens)

plot(minPMF:maxPMF+const, PMF(minPMF:maxPMF, Ap, const), type="p", main="PMF of cycles for model")


dat<-sample(x = c(minPMF:maxPMF)+const, 100000, replace = T, prob = PMF(minPMF:maxPMF, Ap, const)) 
datsmol<-sample(x = c(minPMF:maxPMF)+const, 2000, replace = T, prob = PMF(minPMF:maxPMF, Ap, const)) 
datsmoller<-sample(x = c(minPMF:maxPMF)+const, 200, replace = T, prob = PMF(minPMF:maxPMF, Ap, const)) 

plot(density(dat), main="real data in blue, math prediction in black (100000 cycles)")
points(density(cycleLens), type="l", col="blue")



plot(density(datsmol), main="real data in blue, math prediction in black (2000 cycles)")
points(density(cycleLens), type="l", col="blue")

plot(density(datsmoller), main="real data in blue, math prediction in black (200 cycles)")
points(density(cycleLens), type="l", col="blue")



