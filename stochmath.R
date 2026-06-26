# formulas for mean and variance of the "fancier" model (stochastic levels of Is and As)

expmean_stoch <- function(Ap, As, Ad, Is, Ihl) {
    # expected value of At is hard to "solve" for, so we calculate it manually to depth-3
    # 2**3 possible paths to current time, but note that value is draw-order-dependent
    #   because you cannot decay below zero, so earlier events are decayed "extra"
    # e.g: f1 0 f3 ~= f1 f2 0 b/c f3 will experience 2Ad vs f2 experiencing 1Ad
    E_at <- As + max( (Ap**0 * (1-Ap)**3 * (0*As - 0*Ad)), 0) + # 0 0 0   
                 max( (Ap**1 * (1-Ap)**2 * (1*As - 1*Ad)), 0) + # f 0 0    
                 max( (Ap**1 * (1-Ap)**2 * (1*As - 2*Ad)), 0) + # 0 f 0
                 max( (Ap**1 * (1-Ap)**2 * (1*As - 3*Ad)), 0) + # 0 0 f
                 max( (Ap**3 * (1-Ap)**0 * (3*As - 3*Ad)), 0) + # f f f
                 max( (Ap**2 * (1-Ap)**1 * (2*As - 2*Ad)), 0) + # f f 0
                 max( (Ap**2 * (1-Ap)**1 * (2*As - 3*Ad)), 0) + # 0 f f
               # last case is f 0 f, have to split the sum since the oldest event could amount to nothing
                 (Ap**2 * (1-Ap)**1 * ( max((1*As - 1*Ad),0) + max((1*As-2*Ad),0) ) ) # f 0 f              
    print("###########@@@@@@@@@@@@@@@@")  
    print("E_at:")
    print(E_at)
    print("###########@@@@@@@@@@@@@@@@")
    # expected value of tmin (this language is weird/wrong, since the actual minimal 
    # time is very very small via repeated max-As vs 1-minimal Is, which is not what we are 
    # interested in/not what this is. this is the expected minimal wait time given that As fires
    # at the earliest possible opportunity once It <= As
    # E_tmin <- Ihl * (log(E_at / (E_at + Is)) / log(1/2))
    
    # expected value of Is is just 
    #  sum of (odds of i number of eggs * Is for i number of eggs)
    # note that min number of eggs is +1'd, so pulling a zero -> 1 egg, 1 -> 2 eggs, etc.
    psum <- 0
    i = 0
    while (i<100) {
        valIs <- ((i+1)*Is) / 2
        oddsIs <- dpois(i, lambda=1) 
        psum <- psum + (oddsIs * ( log(E_at / (E_at + valIs)) / log(1/2) ) )
        i<-i+1
    }
    E_tmin <- Ihl * psum
    print("E_tmin:")
    print(E_tmin)
    print("1/Ap:")
    print(1/Ap)
    # expected value of t when you cycle is just the minimum value plus expected wait time for As
    E_t <- 1/Ap + E_tmin
    # this is not quite correct, actually it should be something like:
    # 1/Ap * odds of single As being sufficient + 1/Ap**2 * odds of double As being necessary+sufficient + ...
    return(E_t)
}
