# clear memory
rm(list=ls(all=TRUE))

# set working directory to your folder
setwd("C:\\Users\\Phunmiii_o\\Desktop\\Micheal Project\\Data\\Nigeria\\overnutrition\\csv")
# Load malnutrition data
data <- read.csv("Normal.csv", header=T)


# preallocate a file to save results in
results=matrix(0,35,10)

  # extract proportion malnourished
  Pmal<- as.vector (data[,3])
  # extract mortality data
  qx<-data[,2]
  
  s=length(qx) # nr of ages
  I=diag(s) # identity matrix
  # define the Markov chain
  #P=rbind(cbind(U,rep(0,35)),c(1-colSums(U),1))
  numtypes=2
  
  for(i in 1:numtypes){
    es=rep(1,s) # vector of ones
    # build transition matrix (assigns elements to the off diagonals)
    px=1-qx
    U=matrix(0, nrow=s, ncol=s)
    U[row(U)-1 == col(U) ] <- px[-35]
    U[s,s]=px[35]
    
    # define the Markov chain
    P=rbind(cbind(U,rep(0,s)),c(1-colSums(U),1))
    
    if(i==1){
      # reward for every being malnourished
      m1=rep(1,s)
      m2=rep(1,s)
      m3=rep(1,s)
    } else {
      # reward that not everyone is malnurished
      m1=Pmal
      m2=Pmal
      m3=Pmal   
    }
    
    # Define reward matrices with rewards associated with the transitions
    R1=matrix(0,nrow=s+1,ncol=s+1)
    R2=R3=R1
    # to place entries on the subdiagonal match row to column and then subtract 1 from row
    R1[row(R1)-1 == col(R1)] <-m1
    R1[s+1,s+1]=0
    
    R2[row(R2)-1 == col(R2)] <-m2
    R2[s+1,s+1]=0
    
    R3[row(R3)-1 == col(R3)] <-m3
    R3[s+1,s+1]=0
    ##alternatively, assume fixed rewards:
    #R2fixed=R1^2
    #R3fixed=R1^3
    
    Z=cbind(diag(s),rep(0,s)) # truncation matrix
    e=rep(1,s+1) #vector of ones (36x1)
    
    #R11=Z%*%R1%*%t(Z)
    #R22= Z%*%R2%*%t(Z)
    # moments of the accumulated rewards
    rho1 = solve(I-t(U))%*%(Z%*%((t(P*R1))%*%e))
    #rhoo2 = t(solve(I-(U)))%*%(Z%*%((t(P*R2))%*%e)+2*(t(U*R11))%*%rho1)
    #rhoo3 = t(solve(I-(U)))%*%(Z%*%((t(P*R3))%*%e)+3*(t(U*R22))%*%rho1 +3*(t(U*R11))%*%rhoo2)
    rho2 = solve(I-t(U))%*%(Z%*%((t(P*R2))%*%e)+2*Z%*%(t(P*R1))%*%t(Z)%*%rho1)
    rho3 = solve(I-t(U))%*%(Z%*%((t(P*R3))%*%e)+3*Z%*%(t(P*R2))%*%t(Z)%*%rho1 + 3*Z%*%(t(P*R1))%*%t(Z)%*%rho2);
    
    # Statistics of Longevity
    
    var=rho2-rho1*rho1
    std=sqrt(var)
    cv=std/rho1
    
    vv=var
    pick=which(vv==0)
    vv[pick]=1
    vv=vv^(-3/2)
    vv[pick]=0
    vv=as.vector(vv)
    dd=diag(vv,nrow=length(vv))
    skew=dd%*%(rho3-3*rho2*rho1+ 2*rho1*rho1*rho1)
    
    #Another way to obtain the Coefficient of variation
    #inv.rho1 = 1/(rho1)
    #Diag.rho1=matrix(0, nrow(inv.rho1), nrow(inv.rho1))
    #Diag.rho1[row(Diag.rho1) == col(Diag.rho1) ] <- inv.rho1
    #cvv = Diag.rho1%*%std
    
    # save Longevity stats (type 1) in the first 5 columns, and Healthy LE in the next five
    
    if(i==1){
      results[,1:5]<-cbind(rho1,var,std,cv,skew)
    } else {
      results[,6:10]<-cbind(rho1,var,std,cv,skew)
    }
    
  }
  #write.csv(results[,8], file = "stdfixed.csv")
  #write.csv(results[,10], file = "skewfixed.csv")
  # create age transition matrices
  
  #Create a 35 by 35 matrix with 1's on the subdiagonal
  ones=matrix(1,s-1,1)
  D=matrix(0,s,s)
  D[row(D)-1==col(D)]<-ones
  
library(matrixcalc)
#to obtain the vec operator of D
  vecY=vec(D)
      g=nrow(vecY)
#forming a diagonal matrix with vecY on the diagonals and zeros elswhere
  DvecY=matrix(0, nrow= g, ncol=g)
  DvecY[row(DvecY) == col(DvecY) ] <- vecY
  
#kronecker product of Identity matrix and Vectors of one
   bbD=kronecker(I,es)
#exponential of survival probability vector
  Uexp = exp(-qx)
##forming a diagonal matrix with Uexp on the diagonals and zeros elswhere
  DP=matrix(0, nrow= s, ncol=s)
  DP[row(DP) == col(DP) ] <- Uexp
#obtaining sensitivity of the transition matrix to the mortality schedule
dvecU = -(DvecY%*%bbD%*%DP)


Zer35by1 = as.matrix(rep(0,s))
Zer1by35 = t(Zer35by1)
#creating an identity matrix of order 1
I_1=diag(1)
C1 = kronecker(rbind(I,Zer1by35), rbind(I,Zer1by35))
C2 = kronecker(rbind(I,Zer1by35), rbind(Zer35by1,I_1))
#obtaining the 
kr_Iby1 = kronecker(I,t(es))
dvecP = (C1 - (C2%*%kr_Iby1))%*%dvecU
dvecR = kronecker(t(Z),e)

#forming a diagonal matrix with vecR on the diagonals and zeros elswhere
vecR= vec(R1) #vecR1 =vecR2 =vecR3 since R1=R2=R3
Diag.vecR1=matrix(0, nrow(vecR), nrow(vecR))
Diag.vecR1[row(Diag.vecR1) == col(Diag.vecR1) ] <- vecR

#forming a diagonal matrix with the vec of Markov Chain ,vecP, on the diagonals and zeros elswhere
vecP= vec(P) 
Diag.vecP=matrix(0, nrow(vecP), nrow(vecP))
Diag.vecP[row(Diag.vecP) == col(Diag.vecP) ] <- vecP

#Obtaining V 
V1 = kronecker(t(e),Z)%*%commutation.matrix(36,36)%*%((Diag.vecR1%*%dvecP)+(Diag.vecP%*%dvecR))
V2=V3=V1 #Since R1=R2=R3
#obtaining X
X1 = kronecker(t(rho1),I)%*%commutation.matrix(35,35)%*%(dvecU)
X2 = kronecker(t(rho2),I)%*%commutation.matrix(35,35)%*%(dvecU)
X3 = kronecker(t(rho3),I)%*%commutation.matrix(35,35)%*%(dvecU)

#Sensitivity of first moment of LMP
drho1= t(solve(I-(U)))%*%(V1 + X1) #Also the sensitivity of mean LMP to mortality rate

#Obtaining Rtilde
Rtilde.1 = Z%*%R1%*%t(Z) 
Rtilde.2 = Z%*%R2%*%t(Z) 
Rtilde.3 = Z%*%R3%*%t(Z)

#Obtaining vec of Rtilde
vec.Rtilde1 =vec(Rtilde.1)
vec.Rtilde2 =vec.Rtilde3 =vec.Rtilde1 
#forming a diagonal matrix with the vec of Rtilde on the diagonals and zeros elswhere
Diag.vecRtilde=matrix(0, nrow(vec.Rtilde1), nrow(vec.Rtilde1))
Diag.vecRtilde[row(Diag.vecRtilde) == col(Diag.vecRtilde) ] <- vec.Rtilde1 #Since R1 =R2=R, Diag(Rtilde1) that of 2 and 3

#Obtaining the derivative of vec of R_tilde
dvecR.tilde.1 = kronecker(Z,Z)%*%dvecR
dvecR.tilde.2 = dvecR.tilde.3 = dvecR.tilde.1

#Obataining the vec of The matrix of durvival probabilities U
vecU = vec(U)

#forming a diagonal matrix of vec U with one the diagonals and zeros elswhere
Diag.vecU=matrix(0, nrow(vecU), nrow(vecU))
Diag.vecU[row(Diag.vecU) == col(Diag.vecU) ] <- vecU 

#Obataining W1,1, W1,2
W1.1= kronecker(t(rho1),I)%*%commutation.matrix(35,35)%*%((Diag.vecRtilde%*%dvecU)+(Diag.vecU%*%dvecR.tilde.1))+ t((U%*%Rtilde.1))%*%drho1
W1.2= kronecker(t(rho1),I)%*%commutation.matrix(35,35)%*%((Diag.vecRtilde%*%dvecU)+(Diag.vecU%*%dvecR.tilde.1))+ t((U%*%Rtilde.2))%*%drho1

#Obtaining the sensitivity of the second moment of LMP
drho2 = t(solve(I-(U)))%*%(V2 + 2*W1.1 + X2)

#Obataining W2,1
W2.1= kronecker(t(rho2),I)%*%commutation.matrix(35,35)%*%((Diag.vecRtilde%*%dvecU)+(Diag.vecU%*%dvecR.tilde.1))+ t((U%*%Rtilde.1))%*%drho2
#Obtaining the sensitivities of the third moment of LMP
drho3 = t(solve(I-(U)))%*%(V3 + 3*W1.2 + 3*W2.1 + X3)


#Obtain the diagonal matrix with the first moment on the diagonals and zero elsewhere
#inv.rho1 = 1/(rho1)
Diag.rho1=matrix(0, nrow(rho1), nrow(rho1))
Diag.rho1[row(Diag.rho1) == col(Diag.rho1) ] <- rho1

#Obtain the diagonal matrix with the CV on the diagonals and zero elsewhere
Diag.cv=matrix(0, nrow(cv), nrow(cv))
Diag.cv[row(Diag.cv) == col(Diag.cv) ] <- cv

#Obtain the diagonal matrix with inverse of Standard deviation on the diagonals and zero elsewhere
Inv.std = 1/std
Diag.std.inv=matrix(0, nrow(Inv.std), nrow(Inv.std))
Diag.std.inv[row(Diag.std.inv) == col(Diag.std.inv) ] <- Inv.std

#Obtain the diagonal matrix with inverse of Sthe first moment on the diagonals and zero elsewhere

inv.rho1 = 1/(rho1)
Diag.rho1.inv=matrix(0, nrow(inv.rho1), nrow(inv.rho1))
Diag.rho1.inv[row(Diag.rho1.inv) == col(Diag.rho1.inv) ] <- inv.rho1

#Obtain the diagonal matrix with inverse of the square of first moment on the diagonals and zero elsewhere
inv.rho1.sqr = 1/(rho1**2)
Diag.rho1sq.inv=matrix(0, nrow(inv.rho1.sqr), nrow(inv.rho1.sqr))
Diag.rho1sq.inv[row(Diag.rho1sq.inv) == col(Diag.rho1sq.inv) ] <- inv.rho1.sqr

vecI = vec(I)
Diag.vecI=matrix(0, nrow(vecI), nrow(vecI))
Diag.vecI[row(Diag.vecI) == col(Diag.vecI) ] <- vecI

#Obtaining the sensitivity of the statistics of LMP

dvar.rho = drho2 - 2*Diag.rho1%*%drho1
dSD.rho = 1/2*Diag.std.inv%*%dvar.rho
dCV.rho = (Diag.rho1.inv%*%dSD.rho)- kronecker(t(cv),Diag.rho1.inv)%*%Diag.vecI%*%kronecker(I,es)%*%drho1
#hhhh = (Diag.rho1.inv%*%dSD.rho)- kronecker(std,I)%*%Diag.rho1sq.inv%*%Diag.vecI%*%kronecker(I,es)%*%drho1

#For sensitivity to proportion of individuals that are malnourished
#R1prop= e%*%t(m1)%*%Z
#R2prop=R1prop
#vprop = m2 - (m1*m1)

#Elasticity of LMP statistics to mortality rate

#Obtain the diagonal matrix with the inverse of Variance on the diagonals and zero elsewhere
inv.var = 1/var
Diag.var.inv=matrix(0, nrow(var), nrow(var))
Diag.var.inv[row(Diag.var.inv) == col(Diag.var.inv) ] <- inv.var

#Obtain the diagonal matrix with the inverse of CV on the diagonals and zero elsewhere
inv.cv = 1/cv
Diag.cv.inv=matrix(0, nrow(cv), nrow(cv))
Diag.cv.inv[row(Diag.cv.inv) == col(Diag.cv.inv) ] <- inv.cv


#forming a diagonal matrix with mortality rate on the diagonals and zeros elswhere
D.mortality=matrix(0, nrow=s, ncol=s)
D.mortality[row(D.mortality) == col(D.mortality) ] <- qx

Elasticity.rho1 = Diag.rho1.inv%*%drho1%*% D.mortality          
Elasticity.std = Diag.std.inv%*%dSD.rho%*% D.mortality
Elasticity.var = Diag.var.inv%*%dvar.rho%*% D.mortality
Elasticity.CV = Diag.cv.inv%*%dCV.rho%*% D.mortality

#Plot the sensitivity of stattistics at age 15, 24, 35, 49 for each

#AT ate 15
plot(data[,1],drho1[1,],type="o",ylim=c(-5,6),col="green",xlab ="Age", ylab="Sensitivity of Mean LMP")
plot(data[,1],dvar.rho[1,],type="o",ylim=c(-5,66),col="green",xlab ="Age", ylab="Sensitivity of Variance")

#AT age 24
plot(data[,1],drho1[10,],type="o",ylim=c(-5,7.5),col="green",xlab ="Age", ylab="Sensitivity of Mean LMP at age 24")
plot(data[,1],dvar.rho[1,],type="o",ylim=c(-2.5,63),col="green",xlab ="Age", ylab="Sensitivity of Variance at age 24")

#AT age 35
plot(data[,1],drho1[21,],type="o",ylim=c(-3.5,10),col="green",xlab ="Age", ylab="Sensitivity of Mean LMP at age 35")
plot(data[,1],dvar.rho[21,],type="o",ylim=c(0,49),col="green",xlab ="Age", ylab="Sensitivity of Variance at age 35")

#AT age 49
plot(data[,1],drho1[21,],type="o",ylim=c(-3.5,10),col="green",xlab ="Age", ylab="Sensitivity of Mean LMP at age 35")
plot(data[,1],dvar.rho[21,],type="o",ylim=c(0,49),col="green",xlab ="Age", ylab="Sensitivity of Variance at age 35")

write.csv(drho1, file = "drho1U.csv")
write.csv(dvar.rho, file = "dvarU.csv")
write.csv(dSD.rho, file = "dSDU.csv")
write.csv(dCV.rho, file = "dCVU.csv")

write.csv(Elasticity.rho1, file = "Erho1U.csv")
write.csv(Elasticity.std, file = "EstdU.csv")
write.csv(Elasticity.var, file = "EvarU.csv")
write.csv(Elasticity.CV, file = "ECVU.csv")

