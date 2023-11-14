# Copula packages
library(copula)
library(VineCopula)
library(rvinecopulib)

# Graphical packages
library(ggplot2)
library(reshape2)
library(ggraph)

# Auxiliary package
library(combinat)

#####################################################
### Fit and visualize vine structure
### Input: Pobs.h=pseudo-observations
### Output: vine copula tree (first level)
#####################################################

Plot.Vine=function(Pobs.h){ 
Vine.h=vine(Pobs.h)
    plot(Vine.h$copula,tree=1,var_names="use",edge_labels="family")
}


#####################################################
### Visualize heat maps for model-based
### Spearman's rho and Kendall's tau
### Input: Pobs.h=pseudo-observations 
###        N=sample size (must be large)
###        type="spearman" or "kendall"
### Output: heat map
#####################################################

# Heatmap for model-based Spearman's rho
Heat.Map=function(Pobs.h,N=10^3,type="spearman"){
  vinecopula=vine(Pobs.h)$copula
  Data.sim=rvinecop(N,vinecopula)

  cormat <- round(cor(Data.sim,method=type),2)
  cormat[lower.tri(cormat)]<- NA
  melted_cormat <- reshape2::melt(cormat, na.rm=TRUE)
  
  ggheatmap=ggplot(data = melted_cormat, aes(Var2, Var1, fill = value))+
    geom_tile(color = "white")+
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                         midpoint = 0, limit = c(-1,1), space = "Lab", 
                         name="Correlation") +
    theme_minimal()+ 
    theme(axis.text.x = element_text(angle = 45, vjust = 1,size = 9, hjust = 1))+
    coord_fixed()
  
  ggheatmap.fin=ggheatmap + 
    geom_text(aes(Var2, Var1, label = value), color = "black", size = 4) +
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      panel.grid.major = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank(),
      axis.ticks = element_blank(),
      legend.justification = c(1, 0),
      legend.position = c(0.5, 0.7),
      legend.direction = "horizontal")+
    guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                                 title.position = "top", title.hjust = 0.5))
  
  return(ggheatmap.fin)
}


# Barplot of model-based Spearman's rho (or Kendall's tau) 
# between the first variable (aka Imbalance) 
# against all other variables

Barplot.cor=function(Pobs.h,N=10^3,type="spearman"){
  vinecopula=vine(Pobs.h)$copula
  Data.sim=rvinecop(N,vinecopula)
  
  cormat <- round(cor(Data.sim,method=type),2)
  
  Data.bar <- data.frame(Variable=colnames(Data.sim)[-1],  
                         Correlation=cormat[1,2:ncol(cormat)])
  
  Data.bar$Variable=factor(Data.bar$Variable,levels = Data.bar$Variable)
  
  barplot.aux=ggplot(Data.bar, aes(x=Variable, y=Correlation,fill=Correlation)) + 
    geom_bar(stat = "identity",col="black",size=0.25)+
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                         midpoint = 0, limit = c(-1,1), space = "Lab")+
    theme_minimal()+
    theme(axis.text.x = element_text(angle = 45, vjust = 1,size = 9, hjust = 1))+
    theme(legend.position="none")+
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          panel.border=element_rect(linetype=1,fill=NA),
          axis.ticks = element_blank()) + 
    guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                                 title.position = "top", title.hjust = 0.5))+
    coord_cartesian(ylim=c(-0.4,0.4))
    return(barplot.aux)  
}

#####################################################
### Calculate all induced pairwise association 
### measures (useful for extracting CI's) from B replications
### Input: Pobs.h=pseudo-observations 
###        N.B=sample size (must be large)
###        B=bootstrap replications (for CI's)
### Output: heat map
#####################################################

B.fun=function(Pobs.h,N.B=10^3){
  Vine.h=vine(Pobs.h)
  B.sample=rvinecop(nrow(Pobs.h),Vine.h$copula)
  B.vine=vine(B.sample)
  B.sample.aux=rvinecop(N.B,B.vine$copula)
  B.out.sp =as.vector(cor(B.sample.aux,method="spearman"))
  B.out.ken=as.vector(cor(B.sample.aux,method="kendall"))
  B.out.UTDC=as.vector(fitLambda(B.sample.aux, 
                                 method = "Schmid.Schmidt",lower.tail = FALSE)) 
  B.out.LTDC=as.vector(fitLambda(B.sample.aux,
                                 method = "Schmid.Schmidt",lower.tail = TRUE))
  return(cbind(B.out.sp,B.out.ken,B.out.UTDC,B.out.LTDC))
}

Corr.full=function(Pobs.h,N.B=10^3,B=50){
return(replicate(B,B.fun(Pobs.h,N.B)))
}


#####################################################
### Feature selection (Best Subset Selection) according to Kendall's d.f.
### Input: Pobs.h=pseudo-observations 
###        N=sample size (must be large)
###        type="spearman" or "kendall"
### Output: the subset (of a given cardinality) of variables 
###         that is most correlated (according to Kendall's d.f) 
###         with the first variable
#####################################################

BSS.cor=function(Pobs.h,N=10^4,type="spearman"){
  Vine.h=vine(Pobs.h)
  Copula.h=Vine.h$copula
  Sim.data=rvinecop(N,Copula.h)
  DIM=ncol(Sim.data)
  Matrix1=matrix(rep(1,N*DIM),ncol=DIM)
  
  Cor1.others=cor(Sim.data,method=type)[1,-1]
  Sign=as.numeric(na.omit(sign(Cor1.others)))
  Sim.data.rotated=1-Sim.data

  Sim.data[,1+which(Sign<0)]=Sim.data.rotated[,1+which(Sign<0)]
  
  my_vec=2:DIM
  my_combi = unlist(lapply(1:length(my_vec),    # Get all combinations
                           combinat::combn, 
                           x = my_vec,
                           simplify = FALSE), 
                           recursive = FALSE)
  
  cor.rho.bss=1:(length(my_combi))

      for(i in 1:(length(my_combi))){
        simdata.aux=Matrix1
        simdata.aux[,my_combi[[i]]]=Sim.data[,my_combi[[i]]]
        K.X=C.n(simdata.aux,simdata.aux,smoothing="none")
        cor.rho.bss[i]=cor(Sim.data[,1],K.X,method=type)
      }
  
  my_combi.dim=1:length(my_combi)
      for(i in 1:length(my_combi)){
        my_combi.dim[i]=length(my_combi[[i]])  
      }
  
  cor.rho.max.set=list()
  cor.rho.max.dim=vector()

  lev=0
      for(j in 1:(DIM-1)){
        cor.rho.max.dim[j]=max(cor.rho.bss[my_combi.dim==j])
        cor.rho.max.set[[j]]=my_combi[[lev+which.max(cor.rho.bss[my_combi.dim==j])]]
        lev=lev+sum(my_combi.dim==j)
      }
  
  BSS=list(Sign.cor=Sign,
               Best.subset.cor=cor.rho.max.dim,
               Best.variables.cor=cor.rho.max.set)#,
  return(BSS)
}

#####################################################
### Feature selection (Forward Selection) according to Kendall's d.f.
### Input: Pobs.h=pseudo-observations 
###        N=sample size (must be large)
###        type="spearman" or "kendall"
### Output: the subset (of a given cardinality) of variables 
###         that is most correlated (according to Kendall's d.f) 
###         with the first variable
#####################################################

FS.cor=function(Pobs.h,N=10^4,type="spearman"){
  Vine.h=vine(Pobs.h)
  Copula.h=Vine.h$copula
  Sim.data=rvinecop(N,Copula.h)
  DIM=ncol(Sim.data)
  Matrix1=matrix(rep(1,N*DIM),ncol=DIM)
  
  Cor1.others=cor(Sim.data,method=type)[1,-1]
  Sign=as.numeric(na.omit(sign(Cor1.others)))
  Sim.data.rotated=1-Sim.data

  Sim.data[,1+which(Sign<0)]=Sim.data.rotated[,1+which(Sign<0)]

  Tau.Kn=rep(-999,DIM-1)
  Var.Tau.Kn=rep(-999,DIM-1)
  
  Index=2:DIM           
  Index.search=vector()
  
  for(j in 1:(DIM-1)){
    
    cor.multi.vector=1:length(Index)
    
    for(i in 1:(length(Index))){
      simdata.aux=Matrix1
      simdata.aux[,c(Index.search,Index[i])]=Sim.data[,c(Index.search,Index[i])]
      K.X=C.n(simdata.aux,simdata.aux,smoothing="none")
      cor.multi.vector[i]=cor(Sim.data[,1],K.X,method=type)
    }
    
    loc.max=Index[which.max(abs(cor.multi.vector))]
    
    Var.Tau.Kn[j]=loc.max
    Tau.Kn[j]=cor.multi.vector[which.max(abs(cor.multi.vector))]
    
    Index=Index[Index!=loc.max]
    Index.search=c(Index.search,loc.max)
  }
  FSS=list(Best.variables.cor=Var.Tau.Kn,
           Sign.cor=Sign[Var.Tau.Kn-1],
           Best.subset.cor=Tau.Kn)
  return(FSS)
}


