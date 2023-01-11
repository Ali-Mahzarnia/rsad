boootcca=function(xlist=xlist, 
                  penalty=penalty, 
                  nsim=1000, alpha=0.05, samplesize=dim(xlist[[1]])[1] )
{
  xlistorig=xlist;
  xorig=xlist[[1]];
  yorig=xlist[[2]]
  zorig=xlist[[3]];
  # uboots=matrix(NA,dim(xorig)[2], nsim)
  # vboots=matrix(NA,dim(zorig)[2], nsim)
  corboots=matrix(NA, nsim,1)

  for (ii in 1:nsim) {
    
    bootind=sample(samplesize, replace = TRUE)
    x=xorig[bootind,];y=yorig[bootind,];z=zorig[bootind,];
    
    #remove indeceis with sd =0 o.w the multicca doesnt work with standardize=T
    bootindu=0
    for (i in 1:dim(x)[2]) if(sd(x[,i])==0 ) {bootindu=rbind(bootindu,i);}
    if (length(bootindu)>1){
      bootindu=bootindu[2:dim(bootindu)[1]]
      x=x[,-bootindu] }
    
    bootindw=0
    for (i in 1:dim(y)[2]) if(sd(y[,i])==0 ) {bootindw=rbind(bootindw,i);}
    if (length(bootindw)>1){
      bootindw=bootindw[2:dim(bootindw)[1]]
      y=y[,-bootindw] }
    
    bootindv=0
    for (i in 1:dim(z)[2]) if(sd(z[,i])==0 ) {bootindv=rbind(bootindv,i);}
    if (length(bootindv)>1){
      bootindv=bootindv[2:dim(bootindv)[1]]
      z=z[,-bootindv] }
    
    xlist_plug = list (x, y, z)
    
    outemp <- MultiCCA(xlist_plug, type=c("standard", "standard", "standard"),
                       penalty=penalty, ncomponents=1,  standardize = TRUE)
    corboots[ii]=outemp$cors
    cat('number', ii, 'cor',  corboots[ii] ,'\n')
    
    # uouttemp=uboots[,ii]
    # if(sum(bootindu)>0){
    # uouttemp[bootindu]=0
    # uouttemp[-bootindu]=outemp$u}
    # else {  uouttemp =outemp$u }
    # uboots[,ii]=uouttemp
    # 
    # vouttemp=vboots[,ii]
    # if(sum(bootindv)>0){
    # vouttemp[bootindv]=0
    # vouttemp[-bootindv]=outemp$v}
    # else {  vouttemp= outemp$v}
    # vboots[,ii]=vouttemp
    # 
  }
  
  lowerconfcor= quantile(corboots,  alpha/2)
  higherconfcor=quantile(corboots,  1-alpha/2)
  zstat=mean(corboots)/sd(corboots)
  pvalu=2*pnorm(abs(zstat), lower.tail = FALSE)
returnlist=list( 'corboots'=corboots, 'lowerconfcor'=lowerconfcor,
                'higherconfcor'=higherconfcor,
                'zstat'=zstat, 'pvalu'=pvalu )
return(returnlist)
  
}






set.seed(123)
boot=boootcca(xlist=xlist, 
              penalty=penalty, 
              nsim=1000, alpha=0.05, samplesize=dim(xlist[[1]])[1] )


boot$lowerconfcor
boot$higherconfcor
# boot$zstat
# boot$pvalu

# 
# apply(boot$vboots, 1, quantile, probs = c(0.025))
# apply(boot$vboots, 1, quantile, probs = c(0.975))
# statsv=apply(boot$vboots, 1, mean)/apply(boot$vboots, 1, sd)
# statsv
# 2*pnorm(abs(statsv), lower.tail = FALSE)



#### permutation confindence intervals
# quantile(perm.out$corperms,  0.025)
# quantile(perm.out$corperms,  0.975)
