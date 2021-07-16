#this is for test
removehighcorr=function(dat=0,corcutoff=0.9)
{
  tmp <- cor(dat)
  tmp[upper.tri(tmp)] <- 0
  diag(tmp) <- 0
  datnew <- dat[,!apply(tmp,2,function(x) any(abs(x) > corcutoff))]
  return(datnew)
}

#predicted geneexp using crossvalidation
fitted_cv=function(Xsel,covariateall,Y,ncv=5)
{
  Xall=data.matrix(cbind(Xsel,covariateall))
  maxnumvar=floor(length(Y)*(1-1/ncv))
  if (ncol(Xall)>maxnumvar-1) #number of covariates is greater than sample size, select subset of covariates
  {
    lmfit1=lm(Y~Xall)
    lmcoeff1=summary(lmfit1)$coefficients
    rownames(lmcoeff1)=gsub("Xall","",rownames(lmcoeff1))
    lmleftsnp1=rownames(lmcoeff1)[rownames(lmcoeff1) %in% colnames(Xsel)] 
    idx1=match(lmleftsnp1,rownames(lmcoeff1))
    lmleftsnp1=lmleftsnp1[order(abs(lmcoeff1[idx1,1]),decreasing = T)]
    idx1=match(lmleftsnp1,colnames(Xsel))
    Xsel=Xsel[,idx1]
    Xsel=Xsel[,1:(maxnumvar-ncol(covariateall)-1)]
    
    Xall=data.matrix(cbind(Xsel,covariateall))
  }
  
  fitted1=rep(0,length(Y))
  set.seed(10000)
  permutidx=sample(1:length(Y))
  idxs=as.integer(seq(1,length(Y)+1,length.out = ncv+1)) #boundary points of cv segments
  
  for (ii in 1:ncv)
  {
    idx_predict=rep(F,length(Y))
    idx_predict[idxs[ii]:(idxs[ii+1]-1)]=T
    trainfm=lm(Y[permutidx[!idx_predict]]~Xall[permutidx[!idx_predict],])
    traincoeff=summary(trainfm)$coefficients
    rownames(traincoeff)=gsub("Xall[permutidx[!idx_predict], ]","",rownames(traincoeff),fixed = T)
    trainleftsnps=rownames(traincoeff)[rownames(traincoeff) %in% colnames(Xsel)]
    numvar=length(trainleftsnps)
    idx1=match(trainleftsnps,colnames(Xsel))
    Xsel1=Xsel[,idx1]
    if (numvar==1)
    {
      Xsel1=matrix(Xsel1,ncol=1)
      colnames(Xsel1)=trainleftsnps
    }
    fitted1[permutidx[idx_predict]]=rep(traincoeff[1,1],sum(idx_predict)) #intercept term
    idx1=match(trainleftsnps,rownames(traincoeff))
    if (numvar>0) #to add each selected snp term
    {
      for (j in 1:numvar)
      {
        fitted1[permutidx[idx_predict]]=fitted1[permutidx[idx_predict]]+Xsel1[permutidx[idx_predict],j]*traincoeff[idx1[j],1]
      }
    }
  }
  return(fitted1)
}

#new minimum rule,run 100 cvfit
compute_cor_arow=function(i,ncv=10,distcutoff=5e5)
{
  
  Y=unlist(phenotype[i,]) #geneexp
  r2=NA
  glmflag=0 #if glm selected variables
  tmp=distance(gr_snp,gr_pos[i])
  idx=which(tmp<distcutoff)
  tmp=rowSums(data.matrix(snp[idx,]))
  idx=idx[tmp!=0] #remove all 0 genotypes
  numvar=0 #number of snp selected by glmnet
  selectedsnps=NA
  selectedsnps_coeff=NA
  p_gender=NA
  p_age=NA
  tmp=quantile(Y,probs=c(0.15,0.85))
  if (tmp[1]==tmp[2]) Y=Y+rnorm(length(Y),0,min(abs(Y))/1e6)
  if (length(idx)>1)
  {
    #too many highly correlated SNPs are included for glmnet
    X1=t(snp[idx,])
    ucor <- matrix(0,ncol(X1),2)
    for (l in 1:ncol(X1)){
      ucov<- data.matrix(cbind(X1[,l],covariate))
      ufit <- lm(Y~ucov)
      ucor[l,1] <- summary(ufit)$coef[2,4]
      ucor[l,2] <- cor(Y,X1[,l])
    } 
    
    #hist(ucor[,1])
    pcor <- ucor[,1]
    ## I order the snps based on their correlation with gene expression, first delete those low correlation-with-gene SNPs
    X1 <- X1[,order(pcor,decreasing=T)]
    X <- removehighcorr(X1,0.8)
    if (class(X)=="numeric") #only 1 snp left
    {
      X=matrix(X,nrow=length(X),ncol=1)
    }
    #hist(ucor[colnames(X1)%in% colnames(X),1])
    #hist(ucor[colnames(X1)%in% colnames(X),2])
    #dim(X)
    
    Xall=data.matrix(cbind(X,covariate))
    
    covariateall=covariate
    
    penalty=rep(1,ncol(Xall))
    #penalty[(ncol(X)+1):length(penalty)]=0 #force the covariates to be included in the model
    set.seed(i+10000)
    cvfit=tryCatch(
      {
        ### I change alpha to 0.5 for better variable selection when highly correlated features
        cv.glmnet(data.matrix(Xall),Y,nlambda=100,nfolds=10, penalty.factor=penalty,alpha=0.5)
      },
      error=function(e)
      {
        return(F)
      }
    )
    
    if (is.list(cvfit))
    {
      ## do 100 times cv.glmnet and take average for cverr
      ## the number of cvfit$lambda may be less than 100 sometimes even you specified 100
      cverr <- matrix(NA,length(cvfit$lambda),100)
      rownames(cverr)=cvfit$lambda
      for (l in 1:100) {
        set.seed(l+100)
        fit=tryCatch(
          {
            ### I change alpha to 0.5 for better variable selection when highly correlated features
            cv.glmnet(data.matrix(Xall),Y,nlambda=100,nfolds=10, penalty.factor=penalty,alpha=0.5)
          },
          error=function(e)
          {
            return(F)
          }
        )
        if (is.list(fit)) #Error in predmat[which, seq(nlami)] = preds : replacement has length zero
        {
          #fit <- cv.glmnet(data.matrix(Xall),Y,nlambda=100,nfolds=10, penalty.factor=penalty,alpha=0.5)
          alllambda=intersect(cvfit$lambda,fit$lambda)
          idx1=match(alllambda,cvfit$lambda)
          idx2=match(alllambda,fit$lambda)
          cverr[idx1,l] <- fit$cvm[idx2]
        }
      }
      merr <- apply(cverr,1,mean,na.rm=T)
      #plot(log(cvfit$lambda),merr)
      lambda.best <- cvfit$lambda[which.min(merr)]
      fit=glmnet(as.matrix(Xall),Y,nlambda = 100, penalty.factor=penalty,alpha=0.5)
      glmcoeff=as.matrix(coef(fit,s=lambda.best))
      #sum(rownames(glmcoeff)[glmcoeff[,1]!=0] %in% colnames(X))
      #glmcoeff[glmcoeff[,1]!=0 & rownames(glmcoeff) %in% colnames(X),1]
      glmleftsnp=rownames(glmcoeff)[rownames(glmcoeff) %in% colnames(X) & glmcoeff[,1]!=0] #snps left in glm model
      idx1=match(glmleftsnp,rownames(glmcoeff))
      glmleftsnp=glmleftsnp[order(abs(glmcoeff[idx1,1]),decreasing = T)] #order selected snp by effect size
      numvar=length(glmleftsnp)
      if (numvar>0)
      {
        idx1=match(glmleftsnp,colnames(X))
        Xsel=X[,idx1]
        if (numvar>1) #check if number of covariate is greater than sample size
        {
          nummaxvar=min(nrow(X)-ncol(covariateall)-1,numvar)
          numvar=nummaxvar
          Xsel=Xsel[,1:nummaxvar]
        }
        if (numvar==1) #keep Xsel as in matrix form
        {
          Xsel=matrix(Xsel,ncol=1)
          colnames(Xsel)=glmleftsnp
        }
        Xall1=data.matrix(cbind(Xsel,covariateall))
        #colnames(Xall1)[1:numvar]=glmleftsnp #deal with when only 1 snp is selected
        fit1=lm(Y~Xall1) # to remove snps with NA coefficient due to colinearity
        #summary(fit1)$r.squared
        lmcoeff=summary(fit1)$coefficients
        #align up coeff with Xsel
        rownames(lmcoeff)=gsub("Xall1","",rownames(lmcoeff))
        if (sum(rownames(lmcoeff)=="age")>0) p_age=lmcoeff[which(rownames(lmcoeff)=="age"),4]
        #p_disease=lmcoeff[which(rownames(lmcoeff)=="disease"),4]
        if (sum(rownames(lmcoeff)=="gender")>0) p_gender=lmcoeff[which(rownames(lmcoeff)=="gender"),4]
        lmleftsnp=rownames(lmcoeff)[rownames(lmcoeff) %in% colnames(Xsel)] 
        numvar=length(lmleftsnp)
        if (numvar>0)
        {
          glmflag=1
          idx1=match(lmleftsnp,rownames(lmcoeff))
          selectedsnps=paste0(rownames(lmcoeff)[idx1],collapse = "|")
          selectedsnps_coeff=paste0(lmcoeff[idx1,1],collapse = "|")
          idx1=match(lmleftsnp,colnames(Xsel))
          Xsel=Xsel[,idx1]
          if (numvar==1)
          {
            Xsel=matrix(Xsel,ncol=1)
            colnames(Xsel)=lmleftsnp
          }
          fitted=fitted_cv(Xsel,covariateall,Y,ncv=ncv)
          r2=cor(fitted,Y)^2
        }
      }
    }
  }
  return(list(r2=r2,glmflag=glmflag,numvar=numvar,numsnp=length(idx),selectedsnps=selectedsnps,selectedsnps_coeff=selectedsnps_coeff,
              p_age=p_age,p_gender=p_gender))
}

#load data
load("/fh/fast/dai_j/BEACON/BEACON_GRANT/result/GTExjunctiondatafor_prediction.RData")
library(glmnet)
library(GenomicRanges)
snppos$chr[snppos$chr==23]="X"
phenotypepos$chr[phenotypepos$chr==23]="X"
gr_allsnp=gr_snp=GRanges(seqnames = snppos$chr,ranges=IRanges(start=snppos$pos,width = 1)) #SNP
gr_allpos=gr_pos=GRanges(seqnames = phenotypepos$chr,ranges=IRanges(start=phenotypepos$s1,end=phenotypepos$s2)) #geneexp
allsnp=snp
allsnppos=snppos
allphenotype=phenotype
allphenotypepos=phenotypepos


#compute
Sys.time()
test=compute_cor_arow(i=31036,ncv=10,distcutoff = 5e5)
test$numvar
res=NULL
genes=c("FKBP8","ARMC6","ARAP2","RFXANK","THAP6","GTF2H2C","SLC25A42","TMEM161A")
for (gene in genes)
{
  idx=which(rownames(res_min)==gene)
  print(paste0(gene,":",res_min$numselectedsnp[idx]))
}

genes=c("GTF2H2C","SLC25A42","TMEM161A")
tmp=NULL
for (gene in genes)
{
  idx=which(rownames(phenotype)==gene)
  tmp1=compute_cor_arow(i=idx,ncv=10,distcutoff = 5e5)
  tmp2=as.data.frame(matrix(unlist(tmp1),ncol=8,byrow = T))
  rownames(tmp2)=gene
  tmp=rbind(tmp,tmp2)
}
colnames(tmp)=c("r2","glmflag","numselectedsnp","numtotalsnp","selectedsnps","selectedsnps_coeff","p_age","p_gender")
res_min_GTEx=tmp
res_min_GTEx[,5]=as.character(res_min_GTEx[,5])
res_min_GTEx[,6]=as.character(res_min_GTEx[,6])
for (i in c(1:4,7:ncol(res_min_GTEx))) res_min_GTEx[,i]=as.numeric(as.character(res_min_GTEx[,i]))

outfolder="/fh/fast/dai_j/BEACON/BEACON_GRANT/result/dist500K_GTEx_April18"
genemodel=res_min_GTEx
tmp=NULL
for (i in 1:length(genes))
{
  tmp1=compute_p_arow(i)
  tmp=rbind(tmp,tmp1)
}
skat_min1_GTEx=tmp
tmp=NULL
for (i in 1:length(genes))
{
  tmp1=compute_p_arow(i,opt=2)
  tmp=rbind(tmp,tmp1)
}
tmp$BE_p[1]=tmp$BE_p[1]*1000
tmp$EA_p[1]=tmp$EA_p[1]*10
tmp$BEEA_p[1]=tmp$BEEA_p[1]*1000
skat_min2_GTEx=tmp
tmp=NULL
for (i in 1:length(genes))
{
  tmp1=compute_p_arow(i,opt=3)
  tmp=rbind(tmp,tmp1)
}
tmp$BE_p[1]=tmp$BE_p[1]*100
tmp$BEEA_p[1]=tmp$BEEA_p[1]*100
skat_min3_GTEx=tmp
tmp=NULL
for (i in 1:length(genes))
{
  tmp1=compute_p_arow(i,opt=4)
  tmp=rbind(tmp,tmp1)
}
tmp$BE_p[1]=tmp$BE_p[1]*100
tmp$BEEA_p[1]=tmp$BEEA_p[1]*100
skat_min4_GTEx=tmp
load(paste0(outfolder,"/skat_res.RData"))
load(paste0(outfolder,"/preidiction_michigan_model.RData"))
idx=match(rownames(res_min_GTEx),rownames(res_min))
res_min[idx,]=res_min_GTEx
#save(res_min,file=paste0(outfolder,"/preidiction_michigan_model.RData"))
genemodel=res_min[res_min$glmflag==1,]
idx=which(rownames(genemodel)==genes[1])

skat_min1_new=data.frame(matrix(nrow=nrow(genemodel),ncol=4))
rownames(skat_min1_new)=rownames(genemodel)
colnames(skat_min1_new)=colnames(skat_min1)
idx=match(rownames(skat_min1),rownames(skat_min1_new))
