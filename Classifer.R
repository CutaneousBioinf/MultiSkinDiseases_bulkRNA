MA.classification.loo <- function(tempinput.pheno, tempinput.exp, tempinput.genes){
  tempclass <- unique(tempinput.pheno[,2])
  tempsamples <- tempinput.pheno[,1]
  
  NB.result <- matrix(numeric(), length(tempsamples), length(tempclass))
  rownames(NB.result) <- tempsamples
  colnames(NB.result) <- tempclass
  
  classprior <- rep(1/length(tempclass),length(tempclass))
  names(classprior) <- tempclass
  
  ### LOO
  for (i in 1:dim(tempinput.pheno)[1]){
    print(i)
    temp.p <- tempinput.pheno[-i,]
    temp.e <- tempinput.exp[tempinput.genes,-i]
    s <- tempsamples[i]
    
    ####################################################################
    ### genes x class parameters: mean / variance
    class.para <- list()
    for (c in names(classprior)){
      tempresult <- temp.e[, temp.p[,2]==c]
      class.para[[c]] <- t(apply(tempresult,1,function(x){c(mean(x),sd(x)^2)}))
    }
    
    ####################################################################
    ### pick genes that do not have variance =0 in any of the class
    tempnotgenes <- c()
    for (c in names(class.para)){
      tempnotgenes <- c(tempnotgenes, rownames(class.para[[c]])[class.para[[c]][,2]==0])
    }
    tempnotgenes <- unique(tempnotgenes)
    tempgenes <- setdiff(tempinput.genes, tempnotgenes)
    
    ####################################################################
    ### samples x class log likelihood
    class.logL <- matrix(numeric(), 1, length(tempclass))
    colnames(class.logL) <- tempclass
    
    for (c in names(classprior)){
        tempdnorm <- dnorm(tempinput.exp[tempgenes,s], mean=class.para[[c]][tempgenes,1], sd=sqrt(class.para[[c]][tempgenes,2]))
        class.logL[1,c] <- sum(log(tempdnorm))
    }
    
    ####################################################################
    ### result samples x class Log PriorLiklihood
    class.LogPL <- matrix(numeric(), 1, length(tempclass))
    colnames(class.LogPL) <- tempclass
    
    class.LogPL[1,] <- (log(classprior)+class.logL[1,])
    
    
    ####################################################################
    ### result samples x class posterior probability
    class.pp <- matrix(numeric(), 1, length(tempclass))
    colnames(class.pp) <- tempclass
    
    ### log(p1 + p2) = log(1 + p1/p2) + log(p2) = log(p2/p2 + p1/p2) + log(p2) = log(p2+p1) – log(p2) + log(p2)
    class.pp[1,] <- exp(class.LogPL[1,])/sum(exp(class.LogPL[1,]))
    
    ##################################
    NB.result[s,] <- class.pp[1,]
  }
  NB.result
}


### sites.anno stores all different tissue types annotations
s="skin"
tempp <- table(site.anno[[s]]$Disease)
tempdiseases <- names(tempp)[tempp>=5]

### site.DEGList_scale_voom_invnorm stores the normalized expression data for different tissues
temp.data <- site.DGEList_scale_voom_invnorm[[s]][,site.anno[[s]]$Disease %in% tempdiseases]
temp.pheno <- cbind(colnames(temp.data),as.character(site.anno[[s]]$Disease[site.anno[[s]]$Disease %in% tempdiseases]))


### pick the top 15 genes induced by the different cytokines: cytokineinduced 

tempcyt <-  lapply(names(cytokineinduced), function(x) { data.frame(cytokineinduced[[x]]) %>%mutate(cytokine = x, gene=rownames(cytokineinduced[[x]])) }) 
temp.cyt.upgenes <- tibble(as.data.frame(do.call(rbind,tempcyt))) %>% filter(lfcMLE>=1) %>% arrange(pvalue) %>% group_by(cytokine) %>% slice(1:15)

### using top 15 genes for each cytokine
res.naivebayes <- MA.classification.loo(temp.pheno, temp.data, intersect(rownames(temp.data),unique(temp.cyt.upgenes$gene)))

res.naivebayes.rank <- t(apply(res.naivebayes,1,function(x){rank(-x)}))
res.naivebayes.rank.truepheno <- character(dim(temp.pheno)[1])

for (i in 1:dim(temp.pheno)[1]){
res.naivebayes.rank.truepheno[i] <- res.naivebayes.rank[i,temp.pheno[i,2]]
}

