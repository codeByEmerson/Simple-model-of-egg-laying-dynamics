
Is=100
As=10#10

Ihl=10

Id=1
Ad=5

Ap=0.15 # odds that A 'fires'


mean <- 1/Ap + ( (Ihl * log(As/Is) ) / log(1/2) )

E_x_2 <- ( (1/Ap)**2 + ( (Ihl * log(As/Is) ) / log(1/2) ))**2

#VARIANCE IN VALUE OF A IS (As**2)*Ap - ((As*Ap)**2)
#THIS IS NOT VARIANCE IN EXPECTED WAIT TIME TO VALUE OF A=As!!!

varAp <- (Ap**2 * (1/ (1-(1-Ap**2))**2 ) ) - (1/Ap)**2

E_x2 <- ( (1/Ap**2) + ( (Ihl * log(As/Is) ) / log(1/2) ))**2

varX <- (1/Ap**2)

print(mean)
print(varX)
