phylostruct<-function(samp,tree,env=NULL,metric=c("psv","psr","pse","psc","sppregs"),null.model=c("frequency","richness","both"),runs=10,alpha=0.05,fam="binomial"){

  metric<-match.arg(metric)
  null.model<-match.arg(null.model)
  if(metric=="sppregs")
  {
  nulls<-t(replicate(runs,sppregs(randomizeSample(samp,null.model=null.model),env,tree,fam=fam)$correlations))
  obs<-sppregs(samp,env,tree,fam=fam)$correlations
  mean.null<-apply(nulls,2,mean)
  quantiles.null<-t(apply(nulls,2,quantile,probs=c(alpha/2,1-(alpha/2))))
  
  return(list(metric=metric,null.model=null.model,runs=runs,obs=obs,mean.null=mean.null
                ,quantiles.null=quantiles.null,phylo.structure=NULL,nulls=nulls))


  } else {

    nulls<-switch(metric,
                       psv = replicate(runs,mean(psv(randomizeSample(samp,null.model=null.model),tree,compute.var=FALSE)[,1])),
                       psr = replicate(runs,mean(psr(randomizeSample(samp,null.model=null.model),tree,compute.var=FALSE)[,1])),
                       pse = replicate(runs,mean(pse(randomizeSample(samp,null.model=null.model),tree)[,1])),
                       psc = replicate(runs,mean(psc(randomizeSample(samp,null.model=null.model),tree)[,1])))
    quantiles.null<-quantile(nulls,probs=c(alpha/2,1-(alpha/2)))
    mean.null<-mean(nulls)
    mean.obs<-switch(metric,
                       psv = mean(psv(samp,tree,compute.var=FALSE)[,1]),
                       psr = mean(psr(samp,tree,compute.var=FALSE)[,1]),
                       pse = mean(pse(samp,tree)[,1]),
                       psc = mean(psc(samp,tree)[,1]))

    if(mean.obs<=quantiles.null[1])
    {phylo.structure="underdispersed"
    } else {if(mean.obs>=quantiles.null[2]){
    phylo.structure="overdispersed"} else {phylo.structure="random"}
    }
    
    return(list(metric=metric,null.model=null.model,runs=runs,mean.obs=mean.obs,mean.null=mean.null
                ,quantiles.null=quantiles.null,phylo.structure=phylo.structure,null.means=nulls))
  }
}
