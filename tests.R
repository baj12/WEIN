#create a random matrix with 100 samples and 100 genes

# 168 is the number of all combinations (cond * sex * age * cmv * mutation)
# 6 is the number of replicats
library(DESeq2)
replicates = 600
resultingDDS = list()
for(replicates in c(600)){
  nSamples = 168 * replicates
  # nSamples = 4 * replicates
  nGenes = 100
  targetExpression = 100
  cnts <- matrix(rnbinom(n=nSamples * nGenes, mu=targetExpression, size=nSamples), ncol=nSamples)
  cnts[5:50,] <- matrix(rnorm(n=nSamples * 46, mean=targetExpression, sd=1), ncol=nSamples)
  
  rownames(cnts) = paste0("gene",1:100)
  
  
  # create some special cases
  # no change , no expression
  cnts[1,] = rnbinom(n=nSamples, mu=2, size=nSamples)
  rownames(cnts)[1] = "noChange"
  
  # no change , high expression
  cnts[2,] = 
    rnbinom(n=nSamples, mu=100, size=nSamples)
  rownames(cnts)[2] = "noChange2"
  
  
  # we want one control / stim
  # sex male/female
  # age: 20 30 40 50 60 70 80
  # cmv pos/neg
  # mutation A T G
  cond <- data.frame(control = factor(rep(c("control", "stim1"), each=nSamples/2)))
  cond$cmv = factor(x = "pos",levels = c("pos", "neg"))
  cond$sex = factor(x = "male",levels = c("male", "female"))
  cond$age = factor(x = "20",levels = c("20", "30", "40", "50", "60", "70", "80"))
  cond$mutation = factor(x = "A",levels = c("A", "T", "G"))
  
  condFactors = c(1, 100)
  sexFactors = c(1, 2)
  ageFactors = c(1, 2, 4, 6, 8, 10, 12)
  cmvFactors = c(1, 2)
  mutFactors = c(1, 10 ,20)
  
  names(condFactors) = levels(cond$control)
  names(cmvFactors) = levels(cond$cmv)
  names(sexFactors) = levels(cond$sex)
  names(ageFactors) = levels(cond$age)
  names(mutFactors) = levels(cond$mutation)
  
  
  idx = 1
  mat = matrix(nrow=0,ncol = 5)
  # mat = matrix(nrow=0,ncol = 2)
  for(co in levels(cond$control)){
    for(se in levels(cond$sex)){
    for(ag in levels(cond$age)){
    for(cm in levels(cond$cmv)){
      for(mu in levels(cond$mutation)){
      mat = rbind(mat, (matrix(rep(c(co, se, ag, cm, mu), each=replicates),nrow=replicates)))
      # mat = rbind(mat, (matrix(rep(c(co, cm), each=replicates),nrow=replicates)))
      fact = sexFactors[se] + condFactors[co] + ageFactors[ag] + cmvFactors[cm] + mutFactors[mu]
      # fact = condFactors[co] +cmvFactors[cm] 
      cols = idx:(idx+replicates-1)
      cnts[5:50,cols] = cnts[5:50,idx] + fact * targetExpression
      idx = idx + replicates
          }
        }
      }
    }
  }
  cnts = round(cnts)
  df = data.frame(mat,stringsAsFactors = T)
  colnames(df) = c("condition", "sex", "age", "cmv", "mutation")
  # colnames(df) = c("condition", "cmv")
  
  # object construction
  dds <- DESeqDataSetFromMatrix(cnts, DataFrame(df), ~ condition + cmv + age + sex + mutation +
                                  condition:age + condition:sex + condition:mutation +  condition:cmv +
                                  cmv:age + cmv:sex + cmv:mutation +
                                  age:sex + age:mutation +
                                  sex:mutation)
  #standard analysis
  dds <- DESeq(dds, parallel = T,fitType='local')
  # res <- results(dds)
  # results(dds, contrast=c("condition","control","stim1"))
  resultingDDS[[replicates]] = dds
}

length(resultingDDS)
{
  dds = resultingDDS[[600]]
  res <- results(dds)
  coefTab = as.data.frame(coef(dds[5:50]))
  # apply(coefTab,2,sd)
  2^apply(coefTab,2,mean)
}


# base plot
gene1Counts = plotCounts(dds, gene=geneName, 
                         intgroup=intgrpNoInt, 
                         returnData=T,
                         transform = T)
gene1Counts$name = geneName
gene2Counts = plotCounts(dds, gene="noChange", 
                         intgroup=intgrpNoInt, 
                         returnData=T,
                         transform = T)
gene2Counts$name = "noChange"
gene3Counts = plotCounts(dds, gene="noChange2", 
                         intgroup=intgrpNoInt, 
                         returnData=T,
                         transform = T)
gene3Counts$name = "noChange2"
gene1Counts = rbind(gene1Counts, gene2Counts, gene3Counts)
# normalize (should use rlog???)
# gene1Counts[rownames(gene1Counts), "count"] = log2(gene1Counts[rownames(gene1Counts), "count"]+1)

png(filename = "~/Downloads/3genes.png", width=25, height = 14,units = "cm", res = 400)
p = stripchart(formula(paste(countsName, "~ name")), 
               data=gene1Counts, vertical=TRUE,las=2, method="jitter",
               cex.axis=0.7, main="3 different genes", 
               # ylim=ylim, 
               ylab = "count", pch=20)
dev.off()


gene1Counts[rownames(gene1Counts), "count"] = log2(gene1Counts[rownames(gene1Counts), "count"]+1)
png(filename = "~/Downloads/3genesLog.png", width=25, height = 14,units = "cm", res = 400)
stripchart(formula(paste(countsName, "~ name")), 
           data=gene1Counts, vertical=TRUE,las=2, method="jitter",
           cex.axis=0.7, main="3 different genes", 
           # ylim=ylim, 
           ylab = "log2(count)", pch=20)
dev.off()
fact = sexFactors[se] + condFactors[co] + ageFactors[ag] + cmvFactors[cm] + mutFactors[mu]
# intgrpNoInt = c("condition", "sex", "age" ,"cmv", "mutation" )
intgrpNoInt = c("condition", "cmv")

gene1Counts = plotCounts(dds, gene=geneName, 
                         intgroup=intgrpNoInt, 
                         returnData=T,
                         transform = T)
gene1Counts$col = gene1Counts$condition
gene2Counts = plotCounts(dds, gene=geneName, 
                         intgroup=intgrpNoInt, 
                         returnData=T,
                         transform = T)
gene2Counts$col = gene2Counts$condition
gene2Counts$condition = "all"

p = ggplot(data=rbind(gene2Counts,gene1Counts),
           aes(x=condition, y=count))  +
  geom_jitter(postition=position_jitter(0.1),cex=0.5, aes(color=col))+
  geom_boxplot(inhereit.aes = F, aes(y = count),notch = T, alpha=0)
p
ggsave(filename = "~/Downloads/splitControl.png", plot = p,device = "png",width = 25, height = 14,units = "cm", dpi = 400)
gene1Counts = gene1Counts[sample(nrow(gene1Counts)),]


p = ggplot(data=gene1Counts,
           aes(x=condition, y=count, colour = condition))  +
  geom_jitter(width = 0.1, cex=0.1)+
  geom_boxplot(inherit.aes = T, aes(y = count),notch = T,width  = 0.3, alpha=0 )
ggsave(filename = "~/Downloads/Control.box.png", plot = p,device = "png",width = 25, height = 14,units = "cm", dpi = 400)


p= ggplot(data=gene1Counts,
          aes(x=condition, y=count, colour = condition))  +
  ggdist::stat_slab(
    ## custom bandwidth
    adjust = 0.5,
    ## move geom to the right
    justification = -.2
    ,
    # ## remove slab interval
    .width = 0,
    point_colour = NA
  ) +
  geom_jitter(width = 0.1, cex=0.1)+
  geom_boxplot(inherit.aes = T, aes(y = count),notch = T,width  = 0.3, alpha=0 )
p
ggsave(filename = "~/Downloads/Control.box.hist.png", plot = p,device = "png",width = 25, height = 14,units = "cm", dpi = 400)




gene1Counts$comb = factor(paste(gene1Counts$condition, gene1Counts$mutation, gene1Counts$cmv,sep =":"))
# gene1Counts
p=  ggplot(data=gene1Counts,
           aes(x=comb, y=count, colour = comb))  +
  ggdist::stat_slab(
    ## custom bandwidth
    adjust = 0.5,
    ## move geom to the right
    justification = -.2
    ,
    # ## remove slab interval
    .width = 0,
    point_colour = NA
  ) +
  geom_jitter(width = 0.1, cex=0.1)+
  geom_boxplot(inherit.aes = T, aes(y = count),notch = T,width  = 0.5, alpha=0 )+ 
  theme(axis.text.x = element_text(angle = 60, hjust=1))


p
ggsave(filename = "~/Downloads/mut.cmv.png", plot = p,device = "png",width = 25, height = 14,units = "cm", dpi = 400)

gene1Counts$comb = factor(paste(gene1Counts$condition, gene1Counts$mutation, gene1Counts$cmv, gene1Counts$sex,sep =":"))

p=  ggplot(data=gene1Counts,
           aes(x=comb, y=count, colour = comb))  +
  ggdist::stat_slab(
    ## custom bandwidth
    adjust = 0.5,
    ## move geom to the right
    justification = -.2
    ,
    # ## remove slab interval
    .width = 0,
    point_colour = NA
  ) +
  geom_jitter(width = 0.1, cex=0.1)+
  geom_boxplot(inherit.aes = T, aes(y = count),notch = T,width  = 0.5, alpha=0 )+ 
  theme(axis.text.x = element_text(angle = 60, hjust=1))

p
ggsave(filename = "~/Downloads/mut.cmv.sex.png", plot = p,device = "png",width = 25, height = 14,units = "cm", dpi = 400)

gene1Counts$comb = factor(paste(gene1Counts$condition, gene1Counts$mutation, gene1Counts$cmv, gene1Counts$sex, gene1Counts$age,sep =":"))

p=  ggplot(data=gene1Counts,
           aes(x=comb, y=count, colour = comb))  +
  ggdist::stat_slab(
    ## custom bandwidth
    adjust = 0.5,
    ## move geom to the right
    justification = -.2
    ,
    # ## remove slab interval
    .width = 0,
    point_colour = NA
  ) +
  geom_jitter(width = 0.1, cex=0.1)+
  geom_boxplot(inherit.aes = T, aes(y = count),notch = T,width  = 0.5, alpha=0 )+ 
  theme(axis.text.x = element_text(angle = 60, hjust=1), legend.position = "none")
p
ggsave(filename = "~/Downloads/mut.cmv.sex.age.png", plot = p,device = "png",width = 25, height = 14,units = "cm", dpi = 400)

gene2Counts = gene1Counts[gene1Counts$condition=="control",]
p=  ggplot(data=gene2Counts,
           aes(x=comb, y=count, colour = comb))  +
  ggdist::stat_slab(
    ## custom bandwidth
    adjust = 0.5,
    ## move geom to the right
    justification = -.2
    ,
    # ## remove slab interval
    .width = 0,
    point_colour = NA
  ) +
  geom_jitter(width = 0.1, cex=0.1)+
  geom_boxplot(inherit.aes = T, aes(y = count),notch = T,width  = 0.5, alpha=0 )+ 
  theme(axis.text.x = element_text(angle = 60, hjust=1), legend.position = "none")
p
ggsave(filename = "~/Downloads/contOnly.mut.cmv.sex.age.png", plot = p,device = "png",width = 25, height = 14,units = "cm", dpi = 400)

gene2Counts = gene2Counts[gene2Counts$mutation=="A",]
p=  ggplot(data=gene2Counts,
           aes(x=comb, y=count, colour = comb))  +
  ggdist::stat_slab(
    ## custom bandwidth
    adjust = 0.5,
    ## move geom to the right
    justification = -.2
    ,
    # ## remove slab interval
    .width = 0,
    point_colour = NA
  ) +
  geom_jitter(width = 0.1, cex=0.1)+
  geom_boxplot(inherit.aes = T, aes(y = count),notch = T,width  = 0.5, alpha=0 )+ 
  theme(axis.text.x = element_text(angle = 60, hjust=1), legend.position = "none")
p
ggsave(filename = "~/Downloads/contOnly.mut.A.cmv.sex.age.png", plot = p,device = "png",width = 25, height = 14,units = "cm", dpi = 400)

gene2Counts = gene2Counts[gene2Counts$cmv=="pos",]
p=  ggplot(data=gene2Counts,
           aes(x=comb, y=count, colour = comb))  +
  ggdist::stat_slab(
    ## custom bandwidth
    adjust = 0.5,
    ## move geom to the right
    justification = -.2
    ,
    # ## remove slab interval
    .width = 0,
    point_colour = NA
  ) +
  geom_jitter(width = 0.1, cex=0.1)+
  geom_boxplot(inherit.aes = T, aes(y = count),notch = T,width  = 0.5, alpha=0 )+ 
  theme(axis.text.x = element_text(angle = 60, hjust=1), legend.position = "none")
p
ggsave(filename = "~/Downloads/contOnly.mut.a.cmv.pos.sex.age.png", plot = p,device = "png",width = 25, height = 14,units = "cm", dpi = 400)


gene2Counts = gene2Counts[gene2Counts$sex=="male",]
p=  ggplot(data=gene2Counts,
           aes(x=comb, y=count, colour = comb))  +
  ggdist::stat_slab(
    ## custom bandwidth
    adjust = 0.5,
    ## move geom to the right
    justification = -.2
    ,
    # ## remove slab interval
    .width = 0,
    point_colour = NA
  ) +
  geom_jitter(width = 0.1, cex=0.1)+
  geom_boxplot(inherit.aes = T, aes(y = count),notch = T,width  = 0.5, alpha=0 )+ 
  theme(axis.text.x = element_text(angle = 60, hjust=1), legend.position = "none")
p
ggsave(filename = "~/Downloads/contOnly.mut.a.cmv.pos.sex.m.age.png", plot = p,device = "png",width = 25, height = 14,units = "cm", dpi = 400)

# ~ condition + sex + age + cmv + mutation
geneName = "gene10"
plotCoefficients(dds, geneName, legend=TRUE)

coefTab = as.data.frame(coef(dds[5:50]))

targetExpression * (sexFactors[1] + condFactors[1] + ageFactors[1] + mean(cmvFactors) + median(mutFactors))
mean(cmvFactors) + mean(mutFactors)

# 
# # moderated log2 fold changes
# resultsNames(dds)
# resLFC <- lfcShrink(dds, coef=2, type="apeglm")
# 
# # an alternate analysis: likelihood ratio test
# ddsLRT <- DESeq(dds, test="LRT", reduced= ~ 1)
# resLRT <- results(ddsLRT)
# 
# geneName = "gene10"
# plotCoefficients(dds, geneName, legend=TRUE)

