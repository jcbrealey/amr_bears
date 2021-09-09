library(ggplot2)
library(reshape2)
library(dplyr)
library(RColorBrewer)
library(vegan)
library(pheatmap)
library(cowplot)
library(ccrepe)
set.seed(4242)

####### MANUSCRIPT: Wild animal oral microbiomes reflect the history of  #######
#######             human antibiotics use                                #######
####### Corresponding authors:  Jaelle Brealey (jaelle.brealey@ntnu.no)  #######
#######              Katerina Guschanski (katerina.guschanski@ed.ac.uk)  #######
####### BioRxiv DOI: https://doi.org/10.1101/2020.12.22.423960           ####### 

### R analysis and figure generation ###

##### functions #####
# taxonomy pull functions
get_tax_level_f <- function(level, lin) {
  if (length(grep(level, lin)) == 0) {
    return("NA")
  } else {
    return(gsub(".{1,2}:","",grep(level,lin, value=T)))
  }
}

get_taxonomy_f <- function(taxIDs) {
  tax.df <- data.frame(taxon=as.character(), sk=as.character(), p=as.character(), c=as.character(), o=as.character(), f=as.character(), g=as.character(), s=as.character())
  for (t in taxIDs) {
    lin <- strsplit(readLines(curl::curl(paste0("https://taxonomy.jgi-psf.org/simple/id/sc/",t),"r")), split = ";")[[1]]
    lin.df <- data.frame(taxon=t,
                         sk=get_tax_level_f("sk:",lin),
                         p=get_tax_level_f("p:",lin),
                         c=get_tax_level_f("c:",lin),
                         o=get_tax_level_f("o:",lin),
                         f=get_tax_level_f("f:",lin),
                         g=get_tax_level_f("g:",lin),
                         s=get_tax_level_f("s:",lin))
    tax.df <- rbind(tax.df, lin.df)
  }
  return(tax.df)
}

# abundance filtering function
abundance_filter_f <- function(count.mat, cutoff) {
  #count.mat must be in format of samples = rows and taxa = columns
  #row.names must be sample IDs
  #cutoff must be proportion (e.g. 0.001 = 0.1% relative abundance)
  count.filt <- count.mat
  prop.mat <- as.data.frame(prop.table(as.matrix(count.filt), margin = 1))
  prop.mat$Sample <- row.names(prop.mat)
  prop.mat.m <- melt(prop.mat, by="Sample", value.name = "Proportion", variable.name = "Taxon")
  samples <- as.character(prop.mat$Sample)
  for (s in samples) {
    exclude.taxa <- as.vector(subset(prop.mat.m, Sample==s & Proportion < cutoff)$Taxon)
    count.filt[s,exclude.taxa] = 0
  }
  return(count.filt)
}

# AMR functions for processing CARD blast results
# Add ARO ID to index
add_ARO_f <- function(p,d) {
  a <- aro.index[which(aro.index$Protein.Accession == as.character(p) & aro.index$DNA.Accession == as.character(d)), "ARO.Accession"]
  if (length(a) == 1) {
    return(as.character(a))
  } else {
    return("NA")
  }
}

# Functions for getting one hit per query
hitID_2_ARO_f <- function(a, d) {
  info <- index[which(index$DNA.Accession == as.character(d) & index$ARO.Accession == as.character(a)), c("Protein.Accession","AMR.Gene.Family","Drug.Class","Resistance.Mechanism")]
  if (nrow(info) != 1) {
    info <- data.frame(Protein.Accession="NA",AMR.Gene.Family="NA",Drug.Class="NA",Resistance.Mechanism="NA")
  }
  p <- as.character(info$Protein.Accession)
  g <- as.character(info$AMR.Gene.Family)
  d <- as.character(info$Drug.Class)
  r <- as.character(info$Resistance.Mechanism)
  return(paste(c(p,g,d,r), collapse = "|"))
}

amr_2_ARO_f <- function(amr, s = " ") {
  amr$V12 <- strsplit(as.character(amr$V11), s)
  amr$V11 <- sapply(amr$V12, function(x){x[[1]]})
  amr$V12 <- sapply(amr$V12, function(x){x[[2]]})
  names(amr) <- c("query.name","hit.name","evalue","bitscore","aln.len","per.match","n.mismatch","n.gapopen","hit.len","hit.start.pos","hit.end.pos","sample_species")
  amr$DNA.acc <- sapply(amr$hit.name, function(x){strsplit(as.character(x), "\\|")[[1]][2]})
  amr$ARO.ID <- sapply(amr$hit.name, function(x){strsplit(as.character(x), "\\|")[[1]][5]})
  amr$info <- mapply(hitID_2_ARO_f, amr$ARO.ID, amr$DNA.acc)
  amr$info2 <- strsplit(as.character(amr$info), "\\|")
  amr$Protein.Accession <- sapply(amr$info2, function(x){x[[1]]})
  amr$AMR.Gene.Family <- sapply(amr$info2, function(x){x[[2]]})
  amr$Drug.Class <- sapply(amr$info2, function(x){x[[3]]})
  amr$Resistance.Mechanism <- sapply(amr$info2, function(x){x[[4]]})
  amr <- amr[,c(1:14,17:20)]
  return(amr)
}

bestHit_f <- function(amr.df) {
  amr.df$Drug.Class <- as.character(amr.df$Drug.Class)
  amr.df$keep <- NA
  for (q in unique(amr.df$query.name)) {
    amr.df[which(amr.df$query.name == q),"keep"]  <- amr.df[which(amr.df$query.name == q),"bitscore"] == max(amr.df[which(amr.df$query.name == q),"bitscore"])
    if (nrow(subset(amr.df, query.name == q & keep == "TRUE")) > 1) {
      # multiple hits with same score
      amr.df[which(amr.df$query.name == q & amr.df$keep == "TRUE"),"keep"] <- "multiple"
      temp.q.sub <- subset(amr.df, query.name == q & keep == "multiple")
      if (nrow(temp.q.sub[!duplicated( temp.q.sub[,c("AMR.Gene.Family","Drug.Class","Resistance.Mechanism")]),]) == 1) {
        #only one combination of gene/drug/mech, keep first hit
        amr.df[which(amr.df$query.name == q)[1],"keep"] <- "multiple_T"
        amr.df[which(amr.df$query.name == q)[c(2:length(which(amr.df$query.name == q)))],"keep"] <- "multiple_F"
      } else if (nrow(temp.q.sub[!duplicated( temp.q.sub[,c("AMR.Gene.Family","Drug.Class","Resistance.Mechanism")]),]) > 1) {
        #seems like often it's the drug class that has differs - maybe can combine info somehow?
        if (nrow(temp.q.sub[!duplicated( temp.q.sub[,c("AMR.Gene.Family","Resistance.Mechanism")]),]) == 1) {
          amr.df[which(amr.df$query.name == q)[1],"keep"] <- "multiple_T"
          amr.df[which(amr.df$query.name == q)[c(2:length(which(amr.df$query.name == q)))],"keep"] <- "multiple_F"
          amr.df[which(amr.df$query.name == q)[1],"Drug.Class"] <- paste(unique(amr.df[which(amr.df$query.name == q),"Drug.Class"]), sep = "", collapse = ";")
        }
      }
    }
  }
  return(amr.df)
}

# Function for getting sample by amr.var (e.g. AMR.Gene.Family) and normalising by bacterial read ct
amr_bysample_norm_f <- function(amr.sub, meta, by.var, ct.var, prop = TRUE) {
  bysample <- as.data.frame(table(amr.sub[,by.var], amr.sub$sample_species))
  bysample <- dcast(bysample, formula = Var1 ~ Var2, fun.aggregate = sum, value.var = "Freq")
  names(bysample)[1] <- by.var
  bysample <- bysample[which(rowSums(bysample[,c(2:ncol(bysample))]) > 0),]
  row.names(bysample) <- bysample[,by.var]
  bysample <- as.data.frame(t(bysample[,2:ncol(bysample)]))
  bysample$SampleID <- row.names(bysample)
  if (!prop) {
    # return counts not normalised by read number
    bysample <- merge(bysample, meta, by = "SampleID", all = TRUE)
    return(bysample)
  } else {
    bysample <- merge(meta[,c("SampleID",ct.var)], bysample, by = "SampleID", all = FALSE)
    n <- ncol(bysample)
    bysample.prop <- bysample
    bysample.prop[,3:n] <- bysample[,3:n] / bysample[,ct.var]
    bysample.prop <- merge(bysample.prop[,c(1,3:n)], meta, by = "SampleID", all = TRUE)
    # assume if first column is all NA, rest will be too
    bysample.prop[which(is.na(bysample.prop[,2])),2:(n-1)] <- 0
    return(bysample.prop)
  }
}

# confounding factor stats investigation - only works if outcome.var is a continuous variable
uni_var_stats_f <- function(df.div, outcome.var, meta.vars) {
  stats.df <- data.frame(meta.var=as.character(), test=as.character(),
                         p.value=as.numeric(), stat.val=as.numeric())
  for (pred.var in meta.vars) {
    if (paste(class(df.div[,pred.var]), collapse = ",") %in% c("numeric","integer")) {
      p <- cor.test(df.div[,outcome.var], df.div[,pred.var], method = "spearman", exact = FALSE)$p.value
      e <- cor.test(df.div[,outcome.var], df.div[,pred.var], method = "spearman", exact = FALSE)$estimate[[1]]
      t <- "cor.spearman"
    } else {
      p <- kruskal.test(df.div[,outcome.var] ~ df.div[,pred.var])$p.value
      e <- kruskal.test(df.div[,outcome.var] ~ df.div[,pred.var])$statistic[[1]]
      t <- "rank.kruskal"
    }
    stats.df <- rbind.data.frame(stats.df, cbind.data.frame(meta.var=pred.var, test=t, p.value=p, stat.val=e))
  }
  return(stats.df)
}

# ARG diversity and subsampling
arg_num_f <- function(data, indices) {
  if(class(data) == "character") {
    # long format
    a.div <- length(unique(data[indices]))
  } else {
    # short format
    d <- data[indices,]
    a.div <- specnumber(colSums(d, na.rm = T))
  }
  return(a.div)
}

# subsample for sequencing effort
arg_num_subsample_f <- function(data, strata.var, ARG.var, depth, perm, r) {
  tests <- data.frame(Test=as.numeric(), Variable=as.character())
  for (x in levels(data[,strata.var])) {
    for (i in seq(1,perm,1)){
      t <- sample(data[which(data[,strata.var] == x),ARG.var], depth, replace = r)
      t2 <- arg_num_f(t)
      tests <- rbind.data.frame(tests, data.frame(Test=t2, Variable=x))
    }
  }
  return(tests)
}

# subsample for sample size
arg_num_subsample_SAMPLES_f <- function(data, strata.var, ARG.cols, depth, perm, r) {
  tests <- data.frame(Test=as.numeric(), Variable=as.character())
  for (x in levels(data[,strata.var])) {
    for (i in seq(1,perm,1)){
      t <- sample(data[which(data[,strata.var] == x),"SampleID"], depth, replace = r)
      t2 <- arg_num_f(data[which(data[,"SampleID"] %in% t), ARG.cols])
      tests <- rbind.data.frame(tests, data.frame(Test=t2, Variable=x))
    }
  }
  return(tests)
}


##### meta data ##### 
amr.meta <- read.csv("data/DC_AMR_sample_metadata_210909.csv", sep = ",", header = TRUE)
row.names(amr.meta) <- amr.meta$SampleID

bears.sequenced <- as.character(amr.meta[which(amr.meta$Sample.type == "Bear"),"SampleID"])
summary(amr.meta[bears.sequenced,"Reads.n.nonhost"])

# set up variables
amr.meta$Spec.Period <- factor(amr.meta$Spec.Period, 
                               levels=c("pre1950","1950.70","1970.85","1985.00","2000.20","blanks"), 
                               ordered = TRUE)
amr.meta$Sample.type <- factor(amr.meta$Sample.type,
                               levels=c("Bear","Ext.blank","LP.blank","Swab"))
amr.meta$Ext.batch <- as.factor(amr.meta$Ext.batch)
amr.meta$LP.batch <- as.factor(amr.meta$LP.batch)
amr.meta$IN.batch <- as.factor(amr.meta$IN.batch)

## check qPCR values
ggplot(amr.meta, aes(Sample.type, log(LP.qPCR.copies.ul)))+
  geom_boxplot(outlier.colour = NA)+
  geom_point(aes(colour = Sample.type), size = 3, position = position_jitterdodge())+
  labs(x = "Sample type", y = "LibPrep qPCR concentration")+
  theme_classic() + theme(axis.text = element_text(size = 12, colour = "black"), axis.title = element_text(size = 13),
                          legend.position = "none")

##### SourceTracker results on unfiltered Kraken2 species & genus level assignments ##### 

#             unknown       soil   labcontam    gut       skin    plaque  calculus
#                grey     purple    pink       blue      green    yellow    orange
colors.st <- c("#999999","#5E4FA2","#993299","#3288BD","#ABDDA4","#FFFFBF","#FDAE61")

st.all.m <- melt(amr.meta[,c("SampleID","st.unknown","st.soil","st.labcontam","st.gut","st.skin","st.plaque","st.calculus")],
                 by="SampleID")

# order samples by host and proportion of oral taxa (then skin, then lab, the soil)
st.all.m$SampleID.2 <- factor(st.all.m$SampleID, 
                              levels = unique(st.all.m[order(amr.meta$Sample.type, amr.meta$st.sum.oral, amr.meta$st.skin, amr.meta$st.labcontam, amr.meta$st.soil),"SampleID"]), ordered = TRUE)

tiff("figures/Supp_Figure_S1_SourceTracker.tiff", units = "in", width = 12, height = 5, res = 600, compression = "zip")
ggplot(st.all.m, aes(SampleID.2,value,fill=variable))+
  geom_col(position="stack")+
  labs(x="Sample", y="Proportion of source", fill="Source")+
  scale_fill_manual(values=colors.st[1:7],
                    labels=c("Unknown","Soil","Laboratory","Human gut",
                             "Human skin","Human plaque","Human calculus"))+
  geom_vline(xintercept = 79.5, color = "black")+
  geom_vline(xintercept = 79.5+29, color = "black")+
  geom_vline(xintercept = 79.5+29+8, color = "black")+
  geom_hline(yintercept = 0.05, color = "black", linetype = 3, size = 1)+
  annotate(geom="text", x = 40, y = 1.02, label = "Bear calculus")+
  annotate(geom="text", x = 94, y = 1.02, label = "Extraction blanks")+
  annotate(geom="text", x = 113, y = 1.02, label = "LP")+
  annotate(geom="text", x = 117.5, y = 1.02, label = "S")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.title = element_text(size = 13), axis.text = element_text(color = "black", size = 12), axis.text.x = element_text(angle = 90, size = 9),
        legend.title = element_text(size = 13), legend.text = element_text(size = 12))
dev.off()

##### Bracken oral bacteria processing ##### 
oral.taxa.ids <- as.character(t(read.delim("data/DC_AMR_oral_bacteria_taxonIDs_200306.txt", sep = " ", header = FALSE)))
oral.taxonomy <- read.csv("data/DC_AMR_oral_bacteria_taxonomy_200306.csv", header = TRUE)

k.mcat <- read.delim("data/DC_AMR_kraken2bracken_otu_table_190827.txt", sep = '\t', row.names = 1)

# Bracken pre-filtering totals
amr.meta$Bracken.readct.total <- sapply(amr.meta$SampleID, function(x){ 
  sum(k.mcat[,as.character(x)])
})
amr.meta$Bracken.richness.unfilt <- sapply(as.character(amr.meta$SampleID), function(x) {
  specnumber(k.mcat[,x])
})

# abundance filtering @ 0.05% relative abundance
k.filt <- abundance_filter_f(as.data.frame(t(k.mcat)), 0.0005)
k.filt <- k.filt[,colSums(k.filt) > 0]

amr.meta$Bracken.richness.abundfilt <- sapply(as.character(amr.meta$SampleID), function(x) {
  specnumber(k.filt[x,])
})

amr.meta$Bracken.readct.total.abundfilt <- sapply(as.character(amr.meta$SampleID), function(x){ 
  sum(k.filt[x,]) 
})

# flag bear samples not passing QC
# flag samples with low Bracken read ct totals 
max(subset(amr.meta, Sample.type %in% c("Ext.blank","LP.blank"))$Bracken.readct.total)
#12390
samples.lowseqdepth <- as.character(subset(amr.meta, !Sample.type %in% c("Ext.blank","LP.blank") & 
                                             Bracken.readct.total < max(subset(amr.meta, Sample.type %in% c("Ext.blank","LP.blank"))$Bracken.readct.total))$SampleID)
#none

# flag samples with low oral SourceTracker proportions
samples.loworal <- as.character(subset(amr.meta, Sample.type == "Bear" & st.sum.oral < 0.05)$SampleID)
#22 samples

## ordination showing swabs/blanks different from bear calculus
k.filt.clr <- apply(k.filt, MARGIN = c(1,2), FUN = function(x) x+1)
k.filt.clr <- as.data.frame(t(apply(k.filt.clr, 1, function(x){log(x) - mean(log(x))})))
k.filt.euc <- vegdist(k.filt.clr, method = "euclidean")
k.filt.euc.nmds <- metaMDS(k.filt.euc, k = 3, try = 20, trymax = 150, maxit = 2000)

k.filt.euc.nmds$stress
# 0.1108
k.filt.euc.nmds.scores <- as.data.frame(scores(k.filt.euc.nmds))
k.filt.euc.nmds.scores$SampleID <- row.names(k.filt.euc.nmds.scores)
k.filt.euc.nmds.scores <- merge(k.filt.euc.nmds.scores, amr.meta, by = "SampleID")
k.filt.euc.nmds.scores$Sample.type2 <- factor(sapply(as.character(k.filt.euc.nmds.scores$SampleID), function(x){
  if(x %in% samples.loworal) {
    "Bear_loworal"
  } else {
    as.character(k.filt.euc.nmds.scores[which(k.filt.euc.nmds.scores$SampleID == x), "Sample.type"])
  }
}), levels = c("Bear","Bear_loworal","Swab","Ext.blank","LP.blank"))

tiff("figures/Supp_Figure_S1_NMDS_BearsBlanks.tiff", units = "in", width = 7, height = 7, res = 600, compression = "zip")
ggplot(k.filt.euc.nmds.scores,aes(NMDS1, NMDS2, shape = Sample.type2, color = st.sum.oral))+
  geom_vline(xintercept = 0, linetype = 3, size = 1, color = "grey42") + geom_hline(yintercept = 0, linetype = 3, size = 1, color = "grey42")+
  geom_point(size = 4)+
  scale_shape_manual(values = c(16,15,17,3,4))+
  scale_color_gradient2(mid = "#96bfeb", high = "#fffd57", midpoint = 0.05)+
  scale_x_continuous(limits = c(-110,110)) + scale_y_continuous(limits = c(-110,110))+
  labs(shape = "Sample type", color = "Oral source proportion")+
  coord_fixed(ratio = 1)+
  theme_bw() + theme(panel.grid = element_blank(), legend.position = "right",
                     axis.text = element_text(size = 12, colour = "black"), axis.title = element_text(size = 13), panel.border = element_rect(size = 1.2))
dev.off()

# Stats for supp figure
adonis(k.filt.euc ~ k.filt.euc.nmds.scores$Sample.type + k.filt.euc.nmds.scores$st.sum.oral + k.filt.euc.nmds.scores$Bracken.readct.total.abundfilt, method = "euclidean", permutations = 1000)
"                                                      Df SumsOfSqs MeanSqs F.Model      R2   Pr(>F)    
k.filt.euc.nmds.scores$Sample.type                      3     88946 29648.7  5.8434 0.12810 0.000999 ***
k.filt.euc.nmds.scores$st.sum.oral                      1     26730 26729.5  5.2681 0.03850 0.000999 ***
k.filt.euc.nmds.scores$Bracken.readct.total.abundfilt   1     10391 10390.5  2.0478 0.01496 0.033966 *  
Residuals                                             112    568274  5073.9         0.81844             
Total                                                 117    694341                 1.00000
"

# subset to oral bacteria and remove flagged bear samples
amr.meta.filt <- amr.meta[which(!amr.meta$SampleID %in% c(samples.loworal, samples.lowseqdepth)),]
k.oral <- k.filt[as.character(amr.meta.filt$SampleID),oral.taxa.ids]

# sum species to genus level
k.oral.genus <- as.data.frame(sapply(as.character(unique(oral.taxonomy$g)), function(x){
  taxa.ids <- as.character(oral.taxonomy[which(oral.taxonomy$g == x),"taxon"])
  if (length(taxa.ids) > 1) {
    rowSums(k.oral[,taxa.ids])
  } else {
    k.oral[,taxa.ids]
  }
}))
k.oral.genera.n <- ncol(k.oral.genus)

k.oral.genus$SampleID <- row.names(k.oral.genus)
k.oral.genus$Bracken.oral.ct <- sapply(as.character(k.oral.genus$SampleID), function(x){
  rowSums(k.oral[x,])
})
k.oral.genus$Unknown.genus <- k.oral.genus$Bracken.oral.ct - rowSums(k.oral.genus[,1:k.oral.genera.n])
k.oral.genus.prop <- as.data.frame(prop.table(as.matrix(k.oral.genus[,c(1:k.oral.genera.n,58)]), 1))

amr.meta.filt$Reads.prop.oral.bacteria <- amr.meta.filt$Reads.n.oral.bacteria / amr.meta.filt$Reads.n.nonhost

## plots showing differences between blanks and bears
# absolute read counts (P6unmap)
ggplot(amr.meta.filt, aes(Sample.type, log10(Reads.n.nonhost)))+
  geom_boxplot(outlier.colour = NA)+
  geom_point(aes(colour = Sample.type), size = 3, position = position_jitterdodge())+
  labs(x = "Sample type", y = "No. of nonhost reads (log)")+
  theme_classic() + theme(axis.text = element_text(size = 12, colour = "black"), axis.title = element_text(size = 13),
                          legend.position = "none")

# oral bacteria (prop of P6unmap)
ggplot(amr.meta.filt, aes(Sample.type, Reads.prop.oral.bacteria))+
  geom_boxplot(outlier.colour = NA)+
  geom_point(aes(colour = Sample.type), size = 3, position = position_jitterdodge())+
  labs(x = "Sample type", y = "Proportion of oral bacteria")+
  theme_classic() + theme(axis.text = element_text(size = 12, colour = "black"), axis.title = element_text(size = 13),
                          legend.position = "none")

amr.meta.filt$Sample.type <- factor(amr.meta.filt$Sample.type, levels = c("Bear","Swab","Ext.blank","LP.blank"))
tiff("figures/Supp_Figure_S1_OralBacteriaReads.tiff", units = "in", width = 6, height = 4, res = 600, compression = "zip")
ggplot(amr.meta.filt, aes(Sample.type, Reads.n.oral.bacteria))+
  geom_boxplot(outlier.colour = NA)+
  geom_point(aes(colour = Sample.type), size = 3, position = position_jitterdodge())+
  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  scale_x_discrete(labels = c("Bear calculus","Museum swab","Extraction blank","Library blank"))+
  labs(y = "Number of oral bacteria reads", x = "Sample source")+
  theme_classic() + theme(axis.text = element_text(size = 12, colour = "black"), axis.title = element_text(size = 13),
                          legend.position = "none")
dev.off()

summary(amr.meta.filt[which(amr.meta.filt$Sample.type == "Bear"),"Reads.n.oral.bacteria"])
"   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
   6168  221375  701630 1135500 1527951 8480762 "
amr.meta.filt[which(amr.meta.filt$Sample.type == "Swab"),"Reads.n.oral.bacteria"]
"5486, 82838"
summary(amr.meta.filt[which(amr.meta.filt$Sample.type == "Ext.blank"),"Reads.n.oral.bacteria"])
"   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
   0.00    7.00   18.00   40.28   49.00  245.00 "
summary(amr.meta.filt[which(amr.meta.filt$Sample.type == "LP.blank"),"Reads.n.oral.bacteria"])
"   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  0.000   0.000   0.000   2.125   3.750   8.000 "

## read count by spec period
ggplot(amr.meta.filt, aes(Spec.Period, Reads.n.oral.bacteria))+
  geom_boxplot(outlier.colour = NA)+
  geom_point(aes(colour = Spec.Period), size = 3, position = position_jitterdodge())+
  labs(y = "Number of oral bacteria reads", x = "Time period")+
  theme_classic() + theme(axis.text = element_text(size = 12, colour = "black"), axis.title = element_text(size = 13),
                          legend.position = "none")

##### meta data for bears (no blanks) #####
ua.amr.meta <- amr.meta.filt[which(amr.meta.filt$Sample.type == "Bear"),]

ua.amr.meta$Spec.Period <- droplevels(ua.amr.meta$Spec.Period)
ua.amr.meta$Spec.county <- factor(ua.amr.meta$Spec.county, 
                                  levels = c("Norrbotten","Vasterbotten","Jamtland","Vasternorrland","Dalarna","Gavleborg","Uppsala"))
ua.amr.meta$Sampling.batch <- factor(ua.amr.meta$Sampling.batch,
                                     levels = c("2017.08","2018.08","2018.11"), ordered = T)

sample.id.byyear <- as.character(ua.amr.meta[order(ua.amr.meta$Spec.coll.year),"SampleID"])[c(1:5,55:57,6:54)]

spec.period.cols <- c("#430053","#3b528b","#379b97","#5dc863","#fde61e")
names(spec.period.cols) <- levels(ua.amr.meta$Spec.Period)

# placeholder coordinates for known county but unknown locality
ua.amr.meta[which(!ua.amr.meta$Spec.known.locality & ua.amr.meta$Spec.county == "Dalarna"),
            c("Spec.longitude","Spec.latitude")] <- as.data.frame(matrix(c(60.6,60.6,60.6,14.5,15.0,15.5), ncol = 2))
ua.amr.meta[which(!ua.amr.meta$Spec.known.locality & ua.amr.meta$Spec.county == "Norrbotten"),
            c("Spec.longitude","Spec.latitude")] <- as.data.frame(matrix(c(67.0,67.0,67.0,67.0,21.0,21.5,22,22.5), ncol = 2))

# Basic information on bears for results & methods
summary(ua.amr.meta$Spec.coll.year)
table(ua.amr.meta$Spec.coll.year.type)
table(ua.amr.meta$Spec.sex)
table(ua.amr.meta$Spec.age.years)
summary(as.numeric(ua.amr.meta[which(!ua.amr.meta$Spec.age.years %in% c("Adult","notAdult","unknown")),"Spec.age.years"]))

##### AMR analysis - processing #####
# AMR processing of blast CARD results

# CARD / ARO files
index <- read.delim("CARD_ARO_index/AMR_CARD_mod190211_aro_categories_index.csv")
index <- index[!duplicated(index[,c(2:5)]),]
dim(index) # 2396 x 5
length(unique(index$DNA.Accession)) #2176

aro.index <- read.delim("CARD_ARO_index/AMR_CARD_mod190211_aro_index.csv")
dim(aro.index) #2464
length(unique(aro.index$ARO.Accession)) #2460
subset(as.data.frame(table(aro.index$ARO.Accession)),Freq > 1)
aro.index[which(aro.index$ARO.Accession %in% c("ARO:3000506","ARO:3003072","ARO:3003893","ARO:3003900")),]
# ARO:3003900 entries identical (excluding model IDs)
# ARO:3003893 entries identical (excluding model IDs)
# ARO:3000506 entries have different DNA accessions
# ARO:3003072 entries have different Protein and DNA accessions
aro.index <- aro.index[!duplicated(aro.index[,c("ARO.Accession","Model.Name","ARO.Name","Protein.Accession","DNA.Accession")]),]

index$ARO.Accession <- mapply(add_ARO_f, index$Protein.Accession, index$DNA.Accession)
# 8 "NAs":
temp.nas <- index[which(index$ARO.Accession == "NA"),c("Protein.Accession","DNA.Accession")]
aro.index[which(aro.index$Protein.Accession %in% temp.nas$Protein.Accession),c("ARO.Accession","ARO.Name","Protein.Accession","DNA.Accession")]
# all cases where multiple ARO accessions with one DNA/Protein accession combination - add to index manually using ARO gene family from index and ARO Name from aro.index
temp.nas.i <- row.names(index[which(index$ARO.Accession == "NA"),])
index[temp.nas.i,"ARO.Accession"] <- c("ARO:3003294","ARO:3004335","ARO:3003926","ARO:3004334","ARO:3003287","ARO:3003285","ARO:3004059","ARO:3004060")

## Load CARD blast results
oral <- read.delim("data/DC_AMR_oral_bacteria_CARD_hits_200309.out", sep = "\t", header = FALSE)
oral.aro <- amr_2_ARO_f(oral, " ")
summary(oral.aro$evalue)
"     Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
0.000e+00 0.000e+00 0.000e+00 3.067e-10 2.000e-14 4.100e-08"

oral.best <- bestHit_f(oral.aro)
# 7213 x 19
table(oral.best$keep)
"     FALSE   multiple multiple_F multiple_T       TRUE 
       1662          4         48         32       5467"

## Manual curation of best hit

nrow(oral.best[which(oral.best$keep %in% c("multiple_T","TRUE") & oral.best$AMR.Gene.Family == "NA"), c("query.name","ARO.ID","DNA.acc")])
# 1623 hits where AMR gene family is unknown
temp.aro.nas.ids <- unique(as.character(oral.best[which(oral.best$keep %in% c("multiple_T","TRUE") & oral.best$AMR.Gene.Family == "NA"), "ARO.ID"]))
# 17 - manually add
aro.index[which(aro.index$ARO.Accession %in% temp.aro.nas.ids),]
temp.aro.nas.dna <- as.character(aro.index[which(aro.index$ARO.Accession %in% temp.aro.nas.ids),"DNA.Accession"])
temp.aro.nas.prot <- as.character(aro.index[which(aro.index$ARO.Accession %in% temp.aro.nas.ids),"Protein.Accession"])

write.csv(index[which(index$DNA.Accession %in% temp.aro.nas.dna),],"CARD_ARO_index/Rout_200309_aro_index_orphans.csv", row.names = FALSE, quote = TRUE)
index[which(index$Protein.Accession %in% temp.aro.nas.prot),]
# none in protein list
temp.aro.nas.info <- read.delim("CARD_ARO_index/Rout_200309_aro_index_orphan_manual_info.txt", header = TRUE)
temp.aro.nas.info$AMR.Gene.Family <- as.character(temp.aro.nas.info$AMR.Gene.Family)
temp.aro.nas.info$Drug.Class <- as.character(temp.aro.nas.info$Drug.Class)
temp.aro.nas.info$Resistance.Mechanism <- as.character(temp.aro.nas.info$Resistance.Mechanism)

oral.best[which(oral.best$ARO.ID %in% temp.aro.nas.ids), "AMR.Gene.Family"] <- sapply(oral.best[which(oral.best$ARO.ID %in% temp.aro.nas.ids),"ARO.ID"], function(i){
  temp.aro.nas.info[which(temp.aro.nas.info$ARO.Accession == as.character(i)), "AMR.Gene.Family"]
})
oral.best[which(oral.best$ARO.ID %in% temp.aro.nas.ids), "Drug.Class"] <- sapply(oral.best[which(oral.best$ARO.ID %in% temp.aro.nas.ids),"ARO.ID"], function(i){
  temp.aro.nas.info[which(temp.aro.nas.info$ARO.Accession == as.character(i)), "Drug.Class"]
})
oral.best[which(oral.best$ARO.ID %in% temp.aro.nas.ids), "Resistance.Mechanism"] <- sapply(oral.best[which(oral.best$ARO.ID %in% temp.aro.nas.ids),"ARO.ID"], function(i){
  temp.aro.nas.info[which(temp.aro.nas.info$ARO.Accession == as.character(i)), "Resistance.Mechanism"]
})

oral.best[which(oral.best$keep %in% c("multiple")), c("query.name","AMR.Gene.Family","Drug.Class","Resistance.Mechanism")]
# both hits to RND (abx efflux) - drug classes slightly different
oral.best[1712:1713,"keep"] <- c("multiple_T","multiple_F")
oral.best[1712,"Drug.Class"]  <- paste(unique(oral.best[1712:1713,"Drug.Class"]), sep = "", collapse = ";")
oral.best[4069:4070,"keep"] <- c("multiple_T","multiple_F")
oral.best[4069,"Drug.Class"]  <- paste(unique(oral.best[4197:4198,"Drug.Class"]), sep = "", collapse = ";")

table(oral.best$keep)
"     FALSE multiple_F multiple_T       TRUE 
       1662         50         34       5467"

grep(";",unique(oral.best[which(oral.best$keep %in% c("multiple_T","TRUE")),"AMR.Gene.Family"]), value = TRUE)
# "aminocoumarin resistant parY;aminocoumarin self resistant parY"
# "major facilitator superfamily (MFS) antibiotic efflux pump;resistance-nodulation-cell division (RND) antibiotic efflux pump"
# "vanS;glycopeptide resistance gene cluster"
# what AMR classes are actually present?
sort(unique(oral.best$AMR.Gene.Family))

unique(oral.best[which(oral.best$AMR.Gene.Family == "aminocoumarin resistant parY;aminocoumarin self resistant parY"),"ARO.ID"])
# same ARO.ID = ARO:3003318 - CARD online says aminocoumarin resistant parY
oral.best[which(oral.best$AMR.Gene.Family == "aminocoumarin resistant parY;aminocoumarin self resistant parY"),"AMR.Gene.Family"] <- "aminocoumarin resistant parY"

unique(oral.best[which(oral.best$AMR.Gene.Family == "vanS;glycopeptide resistance gene cluster"),"ARO.ID"])
# ARO:3002941 = both listed by CARD online too, and no other van/glycopeptide reistance clusters listed, leave as is.

unique(oral.best[which(oral.best$AMR.Gene.Family == "major facilitator superfamily (MFS) antibiotic efflux pump;resistance-nodulation-cell division (RND) antibiotic efflux pump"),"ARO.ID"])
# ARO:3000817 = both listed by CARD online too, so will keep as is.

oral.best$sample_species <- gsub("\\_.*","",oral.best$sample_species)

# get into amr.family format
oral.keep <- oral.best[which(oral.best$keep %in% c("TRUE","multiple_T")), c("query.name","sample_species","DNA.acc","ARO.ID","AMR.Gene.Family","Drug.Class","Resistance.Mechanism")]
# any blanks/swabs in the results?
table(oral.keep$sample_species)
# no

oral.amrfam.prop <- amr_bysample_norm_f(oral.keep, ua.amr.meta, "AMR.Gene.Family", "Reads.n.oral.bacteria", prop = TRUE)
oral.amrfam.cts <- amr_bysample_norm_f(oral.keep, ua.amr.meta, "AMR.Gene.Family", "Reads.n.oral.bacteria", prop = FALSE)

ggplot(melt(oral.amrfam.prop[,c(1:23)], by = c("SampleID")), aes(SampleID, value, fill=variable))+
  geom_col(position="stack")+
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 90))

## AMR - total and family proportions
oral.amrfam.cts$Total.AMR.load <- rowSums(oral.amrfam.cts[,c(2:23)])
oral.amrfam.cts$ARG.num <- specnumber(oral.amrfam.cts[,c(2:23)])
oral.amrfam.cts$ARG.shan <- diversity(oral.amrfam.cts[,c(2:23)], index = "shannon")
oral.amrfam.cts[is.na(oral.amrfam.cts$Total.AMR.load),"Total.AMR.load"] <- 0
oral.amrfam.cts[is.na(oral.amrfam.cts$ARG.num),"ARG.num"] <- 0
oral.amrfam.cts[is.na(oral.amrfam.cts$ARG.shan),"ARG.shan"] <- 0
oral.amrfam.cts$Total.AMR.prop <- ifelse(oral.amrfam.cts$Reads.n.oral.bacteria == 0, 0, oral.amrfam.cts$Total.AMR.load / oral.amrfam.cts$Reads.n.oral.bacteria)

# more useful names
names(oral.amrfam.prop)[2:23] <- c("ARG.AAC6p","ARG.ABCF","ARG.parY","ARG.ANT3pp",
                                   "ARG.ileS","ARG.APH3pp","ARG.APH3p","ARG.APH6","ARG.ABC",
                                   "ARG.erm","ARG.fos","ARG.porin","ARG.lipA","ARG.MFS",
                                   "ARG.MFS.RND","ARG.MATE","ARG.OXA","ARG.RND",
                                   "ARG.rpoB","ARG.SAT","ARG.sul","ARG.vanS")
names(oral.amrfam.cts)[2:23] <- c("ARG.AAC6p","ARG.ABCF","ARG.parY","ARG.ANT3pp",
                                  "ARG.ileS","ARG.APH3pp","ARG.APH3p","ARG.APH6","ARG.ABC",
                                  "ARG.erm","ARG.fos","ARG.porin","ARG.lipA","ARG.MFS",
                                  "ARG.MFS.RND","ARG.MATE","ARG.OXA","ARG.RND",
                                  "ARG.rpoB","ARG.SAT","ARG.sul","ARG.vanS")

row.names(oral.amrfam.prop) <- oral.amrfam.prop$SampleID
row.names(oral.amrfam.cts) <- oral.amrfam.cts$SampleID


##### Total AMR analysis #####
summary(oral.amrfam.cts$Total.AMR.load)
"   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
   0.00   16.00   48.00   96.51  101.00  497.00 "

## Total AMR load increases then decreases over time
tiff("figures/Figure_1b_TotalAMRload.tiff", units = "in", width = 5.5, height = 4, res = 600, compression = "zip")
ggplot(oral.amrfam.cts, aes(Spec.Period, Total.AMR.prop))+
  geom_boxplot(outlier.colour = NA)+
  geom_point(aes(fill = Spec.Period, shape = Spec.Period), size = 3, position = position_jitterdodge())+
  scale_shape_manual(values = c(21:25))+
  scale_fill_manual(values = spec.period.cols)+
  geom_vline(xintercept = 1.5, color = "grey50", linetype = 2)+
  geom_vline(xintercept = 2.5, color = "grey50", linetype = 2)+
  geom_vline(xintercept = 3.5, color = "grey50", linetype = 2)+
  geom_vline(xintercept = 4.5, color = "grey50", linetype = 2)+
  scale_x_discrete(drop = TRUE, labels = c("pre 1951","1951-1970","1971-1985","1986-2000","post 2000"))+
  labs(x = "Time period", y = "Total AMR load")+
  theme_classic() + theme(axis.text = element_text(size = 12, colour = "black"), axis.title = element_text(size = 12),
                          legend.position = "none")
dev.off()

# TABLE S1
summary(glm(Total.AMR.prop ~ Spec.Period, data = oral.amrfam.cts, family = quasibinomial(link = "logit"), 
            weights = oral.amrfam.cts$Reads.n.oral.bacteria))
"               Estimate Std. Error  t value Pr(>|t|)    
(Intercept)    -9.39090    0.06215 -151.102  < 2e-16 ***
Spec.Period.L  0.44724    0.14296    3.128  0.00288 ** 
Spec.Period.Q -0.31484    0.14514   -2.169  0.03466 *  
Spec.Period.C -0.12491    0.12856   -0.972  0.33574    
Spec.Period^4  0.29954    0.13864    2.161  0.03535 *  "

summary(glm(Total.AMR.prop ~ Spec.Period + DNA.fragment.length.median.oral, 
            data = oral.amrfam.cts, family = quasibinomial(link = "logit"), 
            weights = oral.amrfam.cts$Reads.n.oral.bacteria))
"                                  Estimate Std. Error t value Pr(>|t|)    
(Intercept)                     -11.307121   0.453542 -24.931  < 2e-16 ***
Spec.Period.L                     0.285282   0.129431   2.204  0.03205 *  
Spec.Period.Q                    -0.348063   0.124342  -2.799  0.00721 ** 
Spec.Period.C                    -0.127048   0.110147  -1.153  0.25411    
Spec.Period^4                     0.263258   0.119101   2.210  0.03159 *  
DNA.fragment.length.median.oral   0.024953   0.005806   4.298 7.77e-05 ***"

# DNA fragment length is important but period still explains a lot of the model

# other non-significant confounding factors:
batch.vars <- c("Sampling.batch","Ext.batch","LP.batch",
                "IN.batch","IN.protocol","Seq.copies.in.pool","Reads.n.raw",
                "Reads.n.nonhost","Reads.prop.host","Reads.prop.human",
                "Reads.n.oral.bacteria","DNA.fragment.length.median.oral",
                "st.sum.oral","st.sum.nonoral")

# SUPP TABLE S2
uni_var_stats_f(oral.amrfam.cts, "Total.AMR.prop", batch.vars)
# significant: DNA.fragment.length.median.oral

summary(glm(Total.AMR.prop ~ DNA.fragment.length.median.oral, data = oral.amrfam.cts, family = quasibinomial(link = "logit"), 
            weights = oral.amrfam.cts$Reads.n.oral.bacteria))
"                                 Estimate Std. Error t value Pr(>|t|)    
(Intercept)                     -11.66530    0.50082 -23.292  < 2e-16 ***
DNA.fragment.length.median.oral   0.02957    0.00633   4.671 1.97e-05 ***"

ggplot(oral.amrfam.cts, aes(DNA.fragment.length.median.oral, Total.AMR.prop))+
  geom_point(aes(fill = Spec.Period, shape = Spec.Period), size = 3)+
  annotate("text", x = 104, y = 0.00015, label = "rho: 0.457\n p: <0.001")+
  coord_cartesian(xlim = c(55, 95), ylim = c(0,0.00018), clip = "off")+
  scale_shape_manual(values = c(21:25), labels = c("pre 1951","1951-1970","1971-1985","1986-2000","post 2000"))+
  scale_fill_manual(values = spec.period.cols, labels = c("pre 1951","1951-1970","1971-1985","1986-2000","post 2000"))+
  labs(x = "Sample median DNA fragment length", y = "Total AMR load", shape = "Time period", fill = "Time period")+
  theme_classic() + theme(axis.text = element_text(size = 12, colour = "black"), axis.title = element_text(size = 13),
                          legend.position = "right")


##### Total AMR and human variables #####
# Plot AMR by geography 
ggplot(oral.amrfam.cts, aes(Spec.latitude, Spec.longitude, fill = Total.AMR.prop, shape = Spec.known.locality))+
  geom_point(size = 3)+
  scale_shape_manual(values = c(21,23))+
  scale_fill_distiller(palette = "Spectral")+
  labs(y = "Longitude", x = "Latitude", fill = "Total AMR load", shape = "Locality known")+
  theme_bw()+ 
  theme(legend.position = "right", legend.text = element_text(size = 12), legend.title = element_text(size = 13),
        axis.text = element_text(colour = "black", size = 12), axis.title = element_text(size = 13))+
  facet_wrap(~ Spec.Period, nrow = 1, labeller=labeller(
    Spec.Period = c(pre1950="pre 1951",`1950.70`="1951-1970",`1970.85`="1971-1985",`1985.00`="1986-2000",`2000.20`="post 2000"))
  )+
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "bold"))

# SUPP TABLE S1
summary(glm(Total.AMR.prop ~ Spec.latitude : Spec.longitude + Spec.Period + DNA.fragment.length.median.oral, 
            data = oral.amrfam.cts[which(oral.amrfam.cts$Spec.known.locality),],
            family = quasibinomial(link = "logit"), 
            weights = oral.amrfam.cts[which(oral.amrfam.cts$Spec.known.locality),"Reads.n.oral.bacteria"]))
"                                  Estimate Std. Error t value Pr(>|t|)    
(Intercept)                     -1.128e+01  4.969e-01 -22.707  < 2e-16 ***
Spec.Period.L                    3.367e-01  1.861e-01   1.809  0.07746 .  
Spec.Period.Q                   -4.138e-01  1.689e-01  -2.450  0.01843 *  
Spec.Period.C                   -1.035e-01  1.309e-01  -0.791  0.43347    
Spec.Period^4                    2.329e-01  1.335e-01   1.744  0.08826 .  
DNA.fragment.length.median.oral  2.223e-02  6.635e-03   3.350  0.00169 ** 
Spec.latitude:Spec.longitude     1.647e-04  3.097e-04   0.532  0.59762"

# SUPP TABLE S2
bear.vars1.stats <- uni_var_stats_f(oral.amrfam.cts, "Total.AMR.prop", c("Spec.Period","Spec.coll.year","Spec.county"))
bear.vars2.stats <- uni_var_stats_f(oral.amrfam.cts[which(oral.amrfam.cts$Spec.known.locality),], "Total.AMR.prop", 
                                    c("Spec.municipality","Spec.longitude","Spec.latitude","Spec.hist.pop.per.sqkm","Spec.HumanFootprint.mean"))
# all ns. Spec.Period p = 0.085, Spec.county = 0.051

# SUPP TABLE S1
summary(glm(Total.AMR.prop ~ Spec.hist.pop.per.sqkm, 
            data = oral.amrfam.cts[which(oral.amrfam.cts$Spec.known.locality),],
            family = quasibinomial(link = "logit"), 
            weights = oral.amrfam.cts[which(oral.amrfam.cts$Spec.known.locality),"Reads.n.oral.bacteria"]))
"                        Estimate Std. Error  t value Pr(>|t|)    
(Intercept)            -9.391385   0.070000 -134.163   <2e-16 ***
Spec.hist.pop.per.sqkm  0.005714   0.002878    1.986   0.0528 .  "

summary(glm(Total.AMR.prop ~ Spec.HumanFootprint.mean, 
            data = oral.amrfam.cts[which(oral.amrfam.cts$Spec.known.locality),],
            family = quasibinomial(link = "logit"), 
            weights = oral.amrfam.cts[which(oral.amrfam.cts$Spec.known.locality),"Reads.n.oral.bacteria"]))
"                          Estimate Std. Error t value Pr(>|t|)    
(Intercept)              -9.329108   0.106716 -87.420   <2e-16 ***
Spec.HumanFootprint.mean -0.004815   0.019757  -0.244    0.809 "

summary(glm(Total.AMR.prop ~ Spec.hist.pop.per.sqkm + Spec.Period, 
            data = oral.amrfam.cts[which(oral.amrfam.cts$Spec.known.locality),],
            family = quasibinomial(link = "logit"), 
            weights = oral.amrfam.cts[which(oral.amrfam.cts$Spec.known.locality),"Reads.n.oral.bacteria"]))
"                        Estimate Std. Error  t value Pr(>|t|)    
(Intercept)            -9.438857   0.085229 -110.746   <2e-16 ***
Spec.hist.pop.per.sqkm  0.002873   0.003372    0.852   0.3989    
Spec.Period.L           0.509902   0.203314    2.508   0.0159 *  
Spec.Period.Q          -0.309794   0.204835   -1.512   0.1376    
Spec.Period.C          -0.093822   0.149940   -0.626   0.5347    
Spec.Period^4           0.215863   0.171565    1.258   0.2150"

summary(glm(Total.AMR.prop ~ Spec.HumanFootprint.mean + Spec.Period, 
            data = oral.amrfam.cts[which(oral.amrfam.cts$Spec.known.locality),],
            family = quasibinomial(link = "logit"), 
            weights = oral.amrfam.cts[which(oral.amrfam.cts$Spec.known.locality),"Reads.n.oral.bacteria"]))
"                          Estimate Std. Error t value Pr(>|t|)    
(Intercept)              -9.425904   0.108429 -86.931   <2e-16 ***
Spec.HumanFootprint.mean  0.003567   0.018001   0.198   0.8438    
Spec.Period.L             0.511227   0.204387   2.501   0.0162 *  
Spec.Period.Q            -0.372151   0.192580  -1.932   0.0598 .  
Spec.Period.C            -0.098176   0.153061  -0.641   0.5246    
Spec.Period^4             0.287969   0.146217   1.969   0.0552 ."

## distribution of human impact scores across samples vs all Sweden
hf.dist <- read.delim("data/HumanFootprint_Sweden_raster.xyz", header = F, sep = " ")
hf.dist <- hf.dist[which(hf.dist$V3 != 100),] #100 == NA (i.e. areas outside of Sweden)
summary(hf.dist$V3)
"   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  0.000   1.000   4.000   6.164   9.000  50.000"

fig.human.a <- ggplot(oral.amrfam.cts, aes(Spec.hist.pop.per.sqkm, Total.AMR.prop))+
  geom_vline(xintercept = 21.6, linetype = 2, size = 1, color = "#006aa8")+
  annotate("text", x = 80, y = 0.00017, label = "Swedish average", color = "#006aa8")+
  geom_point(size = 3)+
  scale_x_continuous(limits = c(0,300))+
  labs(x = "Historical human population density", y = "Total AMR load")+
  theme_classic() + theme(axis.text = element_text(size = 12, colour = "black"), axis.title = element_text(size = 13),
                          legend.position = "none")
# pop density for all Sweden, year 2000
# https://www.statistikdatabasen.scb.se/pxweb/en/ssd/START__BE__BE0101__BE0101C/BefArealTathetKon/
# accessed 2020-12-14 SCB

fig.human.b <- ggplot(oral.amrfam.cts, aes(Spec.HumanFootprint.mean, Total.AMR.prop))+
  geom_vline(xintercept = 6.164, linetype = 2, size = 1, color = "#006aa8")+
  annotate("text", x = 16, y = 0.00017, label = "Swedish average", color = "#006aa8")+
  geom_point(size = 3)+
  scale_x_continuous(limits = c(0,50))+
  labs(x = "Human Footprint Index", y = "Total AMR load")+
  theme_classic() + theme(axis.text = element_text(size = 12, colour = "black"), axis.title = element_text(size = 13),
                          legend.position = "none")

tiff("figures/Figure_3_HumanImpact_vs_AMRload.tiff", units = "in", width = 5, height = 6, res = 600, compression = "zip")
plot_grid(fig.human.a, fig.human.b, nrow = 2, axis = "lr", labels = c("a","b"))
dev.off()

##### Human antibiotics use (based on publicly available data) ##### 
abx.usage <- read.csv("data/Sweden_antibiotics_usage_active_sub_per_kg.csv")
abx.usage$Period <- cut(abx.usage$Year, breaks = c(1800,1950,1970,1985,2000,2020),
                        labels = c("pre1950","1950.70","1970.85","1985.00","2000.20"),
                        ordered_result = T)

tiff("figures/Figure_1a_Huamn_antibiotics_use.tiff", units = "in", width = 5.5, height = 3, res = 600, compression = "zip")
ggplot(melt(abx.usage, id.vars = c("Year","Period"), measure.vars = c("agriculture.total","human.total")), 
       aes(Year, value, group = variable))+
  geom_vline(xintercept = 1950, color = "grey50", linetype = 2)+
  geom_vline(xintercept = 1970, color = "grey50", linetype = 2)+
  geom_vline(xintercept = 1985, color = "grey50", linetype = 2)+
  geom_vline(xintercept = 2000, color = "grey50", linetype = 2)+
  geom_line()+
  geom_point(aes(fill = Period, shape = variable), size = 3, color = "black")+
  scale_shape_manual(values = c(21:22), labels = c("Animals","Humans"))+
  scale_fill_manual(values = spec.period.cols, labels = c("1971-1985","1986-2000","post 2000"))+
  scale_x_continuous(limits = c(1940,2020), n.breaks = 9)+
  labs(x = "Year", y = "Antibiotics use \n(kg active substance)", shape = "", fill = "Time period")+
  theme_classic() + theme(axis.text = element_text(size = 12, colour = "black"), axis.title = element_text(size = 12),
                          legend.position = "none")
dev.off()

##### ARG diversity #####
# diversity by time period - heatmap
oral.div.period <- as.data.frame(t(sapply(levels(oral.amrfam.cts$Spec.Period), function(x){
  colSums(oral.amrfam.cts[which(oral.amrfam.cts$Spec.Period == x),c(2:23)], na.rm = T)
})))

temp.arg.n <- ncol(oral.div.period)
oral.div.period$Spec.Period <- factor(row.names(oral.div.period), levels = levels(oral.amrfam.cts$Spec.Period), ordered = TRUE)
oral.div.period$Total.AMR.load <- rowSums(oral.div.period[,1:temp.arg.n])
oral.div.period$Sample.ct <- sapply(as.character(oral.div.period$Spec.Period), function(x){
  nrow(oral.amrfam.cts[which(oral.amrfam.cts$Spec.Period == x),])
})
oral.div.period$ARG.num <- specnumber(oral.div.period[,1:temp.arg.n])
oral.div.period$ARG.shan=diversity(oral.div.period[,1:temp.arg.n], index = "shannon")
oral.div.period$Reads.n.oral.summed <- sapply(as.character(oral.div.period$Spec.Period), function(x){
  sum(oral.amrfam.cts[which(oral.amrfam.cts$Spec.Period == x),"Reads.n.oral.bacteria"])
})

oral.div.period.hm <- as.matrix(oral.div.period[,1:temp.arg.n] / oral.div.period$Reads.n.oral.summed)
colnames(oral.div.period.hm) <- gsub("ARG.","",colnames(oral.div.period.hm))
oral.div.period.hm.log <- oral.div.period.hm
oral.div.period.hm.log[oral.div.period.hm.log == 0] <- NA
oral.div.period.hm.log <- log10(oral.div.period.hm.log)
oral.div.period.hm.P <- pheatmap(oral.div.period.hm, cluster_rows = F, cluster_cols = T, scale = "row")
oral.div.period.hm.P.clust <- oral.div.period.hm.P$tree_col
arg.hm.ord.names <- colnames(oral.div.period.hm.log[,oral.div.period.hm.P.clust$order])
# manually tweak order
arg.hm.ord.names <- arg.hm.ord.names[c(6:9,5,10:14,3:4,2,15:16,1,17:22)]

# add annotation to ARGs
arg.anno <- data.frame(
  Mechanism=sapply(arg.hm.ord.names, function(x){
    if(x %in% c("rpoB","parY","lipA","ileS","erm","vanS")) { "Target_alteration" }
    else if(x %in% c("RND","MATE","MFS","MFS.RND","ABC")) { "Efflux_pump" }
    else if(x %in% c("ABCF")) { "Target_protection" }
    else if(x %in% c("fos","OXA","AAC6p","APH6","ANT3pp","APH3pp","APH3p","SAT")) { "Antibiotic_inactivation" }
    else if(x %in% c("sul")) { "Target_replacement" }
    else if(x %in% c("porin")) { "Reduced_permeability" }
    else { "NA" }
  }),
  Antibiotic_class=sapply(arg.hm.ord.names, function(x){
    if(x %in% c("RND","MATE","MFS","MFS.RND","ABC","porin","ABCF")) { "Multi-drug" }
    else if(x %in% c("AAC6p","APH6","ANT3pp","APH3pp","APH3p")) { "Aminoglycosides" }
    else if(x %in% c("ileS","erm","SAT")) { "MLS" }
    else if(x %in% c("rpoB")) { "Rifampicin" }
    else if(x %in% c("parY")) { "Aminocoumarin" }
    else if(x %in% c("lipA")) { "Polymyxins" }
    else if(x %in% c("fos")) { "Fosfomycin" }
    else if(x %in% c("OXA")) { "Beta-lactams" }
    else if(x %in% c("vanS")) { "Vancomycin" }
    else if(x %in% c("sul")) { "Sulfonamides" }
    else { "NA" }
  })
)
row.names(arg.anno) <- arg.hm.ord.names

arg.hm.cols1 <- brewer.pal(length(unique(arg.anno$Mechanism)),"Set2")
names(arg.hm.cols1) <- unique(arg.anno$Mechanism)
arg.hm.cols2 <- brewer.pal(length(unique(arg.anno$Antibiotic_class)),"Set3")
names(arg.hm.cols2) <- unique(arg.anno$Antibiotic_class)
arg.hm.cols <- list(Mechanism=arg.hm.cols1, Antibiotic_class=arg.hm.cols2)

tiff("figures/Figure_1d_ARGdiversity_heatmap.tiff", units = "in", width = 10.4, height = 3.9, res = 600, compression = "zip")
pheatmap(oral.div.period.hm.log[,arg.hm.ord.names], cluster_rows = F, cluster_cols = F,
         labels_row = c("pre 1951","1951-1970","1971-1985","1986-2000","post 2000"),
         labels_col = c("AAC(6')","OXA","erm","fos","sul","APH(3'')","vanS","SAT",
                        "ANT(3'')","APH(6')","GBP","MFS/RND","APH(3')","MFS","MATE",
                        "ABC","ileS","pgpB","ABC-F","parY","RND","rpoB"),
         annotation_col = arg.anno, annotation_colors = arg.hm.cols)
dev.off()

## ARG diversity subsampling to control for sequencing effort
# note that "Ua012" "Ua074" "Ua088" are zero so won't be in the subsampling
arg.readbased <- subset(oral.keep, sample_species %in% as.character(ua.amr.meta$SampleID))
arg.readbased$Spec.Period <- sapply(as.character(arg.readbased$sample_species), function(x){
  ua.amr.meta[which(ua.amr.meta$SampleID == x), "Spec.Period"]
})

table(arg.readbased$Spec.Period)
"pre1950 1950.70 1970.85 1985.00 2000.20 
    553     952     795    1558    1643 "

arg.num.subsample.READS <- arg_num_subsample_f(arg.readbased, "Spec.Period", "AMR.Gene.Family", 500, 1000, F)
arg.num.subsample.READS$Spec.Period <- factor(arg.num.subsample.READS$Variable, levels = levels(oral.amrfam.cts$Spec.Period), ordered = T)

tiff("figures/Supp_Figure_S2_ARGdiversity_subsample_reads.tiff", units = "in", width = 5.5, height = 5.5, res = 600, compression = "zip")
ggplot(arg.num.subsample.READS, aes(Spec.Period, Test))+
  geom_boxplot(outlier.colour = NA)+
  geom_point(aes(color = Spec.Period), size = 2, position = position_jitterdodge(jitter.width = 1))+
  scale_color_manual(values = spec.period.cols)+
  scale_y_continuous(limits = c(0,16))+
  scale_x_discrete(labels = c("pre 1951","1951-1970","1971-1985","1986-2000","post 2000"))+
  labs(x = "Time period", y = "ARG family number", title = "Subsampling by sequencing depth")+
  theme_classic() + theme(axis.text = element_text(size = 12, colour = "black"), axis.title = element_text(size = 13),
                          legend.position = "none", axis.text.x = element_text(angle = 0))
dev.off()

## ARG diversity subsampling to control for sample size
arg.num.subsample.SAMPLES <- arg_num_subsample_SAMPLES_f(oral.amrfam.cts, "Spec.Period", c(2:23), 8, 1000, F)
arg.num.subsample.SAMPLES$Spec.Period <- factor(arg.num.subsample.SAMPLES$Variable, levels = levels(oral.amrfam.cts$Spec.Period), ordered = T)

tiff("figures/Figure_1c_ARGdiversity_subsample_samples.tiff", units = "in", width = 3, height = 3, res = 600, compression = "zip")
ggplot(arg.num.subsample.SAMPLES, aes(Spec.Period, Test))+
  geom_boxplot(outlier.colour = NA)+
  geom_point(aes(color = Spec.Period), size = 2, position = position_jitterdodge(jitter.width = 1))+
  scale_color_manual(values = spec.period.cols)+
  geom_vline(xintercept = 1.5, color = "grey50", linetype = 2)+
  geom_vline(xintercept = 2.5, color = "grey50", linetype = 2)+
  geom_vline(xintercept = 3.5, color = "grey50", linetype = 2)+
  geom_vline(xintercept = 4.5, color = "grey50", linetype = 2)+
  scale_y_continuous(limits = c(0,16))+
  scale_x_discrete(labels = c("pre 1951","1951-1970","1971-1985","1986-2000","post 2000"))+
  labs(x = "Time period", y = "ARG family number")+
  theme_classic() + theme(axis.text = element_text(size = 10, colour = "black"), axis.title = element_text(size = 11),
                          legend.position = "none", axis.text.x = element_text(size = 9))
dev.off()

# mean
for (x in levels(arg.num.subsample.SAMPLES$Spec.Period)){
  print(mean(arg.num.subsample.SAMPLES[which(arg.num.subsample.SAMPLES$Spec.Period == x),"Test"]))
}

##### ARG family changes over time #####
# Supp Figure 3 composite - ARG family over time (excluding RND changes)

Sfig.arg.abcf <- ggplot(oral.amrfam.prop, aes(Spec.Period, ARG.ABCF))+
  geom_boxplot(outlier.colour = NA)+
  geom_point(aes(fill = Spec.Period, shape = Spec.Period), size = 2, position = position_jitterdodge())+
  scale_shape_manual(values = c(21:25))+
  scale_fill_manual(values = spec.period.cols)+
  scale_x_discrete(labels = c("pre \n1951","1951-\n1970","1971-\n1985","1986-\n2000","post\n2000"))+
  labs(x = "", y = "ARG family load", title = "ABC-F ribosomal protection protein")+
  theme_classic() + theme(axis.text = element_text(colour = "black", size = 10), plot.title = element_text(size = 11, face = "bold"), legend.position = "none")

Sfig.arg.rpoB <- ggplot(oral.amrfam.prop, aes(Spec.Period, ARG.rpoB))+
  geom_boxplot(outlier.colour = NA)+
  geom_point(aes(fill = Spec.Period, shape = Spec.Period), size = 2, position = position_jitterdodge())+
  scale_shape_manual(values = c(21:25))+
  scale_fill_manual(values = spec.period.cols)+
  scale_x_discrete(labels = c("pre \n1951","1951-\n1970","1971-\n1985","1986-\n2000","post\n2000"))+
  labs(x = "", y = "ARG family load", title = "rifamycin-resistant rpoB")+
  theme_classic() + theme(axis.text = element_text(colour = "black", size = 10), plot.title = element_text(size = 11, face = "bold"), legend.position = "none")

Sfig.arg.ileS <- ggplot(oral.amrfam.prop, aes(Spec.Period, ARG.ileS))+
  geom_boxplot(outlier.colour = NA)+
  geom_point(aes(fill = Spec.Period, shape = Spec.Period), size = 2, position = position_jitterdodge())+
  scale_shape_manual(values = c(21:25))+
  scale_fill_manual(values = spec.period.cols)+
  scale_x_discrete(labels = c("pre \n1951","1951-\n1970","1971-\n1985","1986-\n2000","post\n2000"))+
  labs(x = "", y = "ARG family load", title = "resistant ileS")+
  theme_classic() + theme(axis.text = element_text(colour = "black", size = 10), plot.title = element_text(size = 11, face = "bold"), legend.position = "none")

Sfig.arg.ABC <- ggplot(oral.amrfam.prop, aes(Spec.Period, ARG.ABC))+
  geom_boxplot(outlier.colour = NA)+
  geom_point(aes(fill = Spec.Period, shape = Spec.Period), size = 2, position = position_jitterdodge())+
  scale_shape_manual(values = c(21:25))+
  scale_fill_manual(values = spec.period.cols)+
  scale_x_discrete(labels = c("pre \n1951","1951-\n1970","1971-\n1985","1986-\n2000","post\n2000"))+
  labs(x = "", y = "ARG family load", title = "ABC efflux pump")+
  theme_classic() + theme(axis.text = element_text(colour = "black", size = 10), plot.title = element_text(size = 11, face = "bold"), legend.position = "none")

Sfig.arg.lipA <- ggplot(oral.amrfam.prop, aes(Spec.Period, ARG.lipA))+
  geom_boxplot(outlier.colour = NA)+
  geom_point(aes(fill = Spec.Period, shape = Spec.Period), size = 2, position = position_jitterdodge())+
  scale_shape_manual(values = c(21:25))+
  scale_fill_manual(values = spec.period.cols)+
  scale_x_discrete(labels = c("pre \n1951","1951-\n1970","1971-\n1985","1986-\n2000","post\n2000"))+
  labs(x = "", y = "ARG family load", title = "lipid A phosphotase (pgpB)")+
  theme_classic() + theme(axis.text = element_text(colour = "black", size = 10), plot.title = element_text(size = 11, face = "bold"), legend.position = "none")

Sfig.arg.parY <- ggplot(oral.amrfam.prop, aes(Spec.Period, ARG.parY))+
  geom_boxplot(outlier.colour = NA)+
  geom_point(aes(fill = Spec.Period, shape = Spec.Period), size = 2, position = position_jitterdodge())+
  scale_shape_manual(values = c(21:25))+
  scale_fill_manual(values = spec.period.cols)+
  scale_x_discrete(labels = c("pre \n1951","1951-\n1970","1971-\n1985","1986-\n2000","post\n2000"))+
  labs(x = "", y = "ARG family load", title = "aminocoumarin resistant parY")+
  theme_classic() + theme(axis.text = element_text(colour = "black", size = 10), plot.title = element_text(size = 11, face = "bold"), legend.position = "none")


# SUPP TABLE S1
summary(glm(ARG.ABCF ~ Spec.Period + DNA.fragment.length.median.oral, 
            data = oral.amrfam.prop, family = quasibinomial(link = "logit"),
            weights = oral.amrfam.prop$Reads.n.oral.bacteria))
"                                 Estimate Std. Error t value Pr(>|t|)    
(Intercept)                     -14.95550    1.04482 -14.314  < 2e-16 ***
Spec.Period.L                    -0.86978    0.32793  -2.652  0.01063 *  
Spec.Period.Q                    -0.88825    0.32500  -2.733  0.00860 ** 
Spec.Period.C                    -0.24320    0.22208  -1.095  0.27863    
Spec.Period^4                    -0.65217    0.29116  -2.240  0.02948 *  
DNA.fragment.length.median.oral   0.03626    0.01324   2.738  0.00849 ** "

summary(glm(ARG.parY ~ Spec.Period + DNA.fragment.length.median.oral, data = oral.amrfam.prop, family = quasibinomial(link = "logit"),
            weights = oral.amrfam.prop$Reads.n.oral.bacteria))
"                                 Estimate Std. Error t value Pr(>|t|)    
(Intercept)                     -14.43248    1.29357 -11.157 2.71e-15 ***
Spec.Period.L                    -0.78011    0.40095  -1.946   0.0572 .  
Spec.Period.Q                    -0.67079    0.40180  -1.669   0.1012    
Spec.Period.C                    -0.63435    0.29379  -2.159   0.0356 *  
Spec.Period^4                    -0.43463    0.37515  -1.159   0.2520    
DNA.fragment.length.median.oral   0.03440    0.01642   2.095   0.0411 *  "

summary(glm(ARG.ileS ~ Spec.Period + DNA.fragment.length.median.oral, data = oral.amrfam.prop, family = quasibinomial(link = "logit"),
            weights = oral.amrfam.prop$Reads.n.oral.bacteria))
"                                 Estimate Std. Error t value Pr(>|t|)    
(Intercept)                     -11.31960    1.07561 -10.524  2.2e-14 ***
Spec.Period.L                     0.25926    0.27559   0.941 0.351269    
Spec.Period.Q                     0.50284    0.34627   1.452 0.152587    
Spec.Period.C                    -1.06377    0.30045  -3.541 0.000863 ***
Spec.Period^4                    -0.13027    0.40934  -0.318 0.751598    
DNA.fragment.length.median.oral  -0.01639    0.01419  -1.154 0.253679 "

summary(glm(ARG.ABC ~ Spec.Period + DNA.fragment.length.median.oral, data = oral.amrfam.prop, family = quasibinomial(link = "logit"),
            weights = oral.amrfam.prop$Reads.n.oral.bacteria))
"                                 Estimate Std. Error t value Pr(>|t|)    
(Intercept)                     -11.65804    1.68617  -6.914 7.44e-09 ***
Spec.Period.L                    -0.51576    0.33050  -1.561   0.1248    
Spec.Period.Q                     0.99492    0.39578   2.514   0.0151 *  
Spec.Period.C                    -0.09993    0.42468  -0.235   0.8149    
Spec.Period^4                     0.70229    0.48756   1.440   0.1559    
DNA.fragment.length.median.oral  -0.03512    0.02266  -1.550   0.1273    "

summary(glm(ARG.lipA ~ Spec.Period + DNA.fragment.length.median.oral, data = oral.amrfam.prop, family = quasibinomial(link = "logit"),
            weights = oral.amrfam.prop$Reads.n.oral.bacteria))
"                                Estimate Std. Error t value Pr(>|t|)    
(Intercept)                     -9.40877    2.09072  -4.500 3.96e-05 ***
Spec.Period.L                    0.05656    0.46034   0.123    0.903    
Spec.Period.Q                    0.81644    0.50988   1.601    0.116    
Spec.Period.C                   -0.03406    0.59551  -0.057    0.955    
Spec.Period^4                    0.96650    0.61926   1.561    0.125    
DNA.fragment.length.median.oral -0.03730    0.02811  -1.327    0.190    "

summary(glm(ARG.RND ~ Spec.Period + DNA.fragment.length.median.oral, data = oral.amrfam.prop, family = quasibinomial(link = "logit"),
            weights = oral.amrfam.prop$Reads.n.oral.bacteria))
"                                 Estimate Std. Error t value Pr(>|t|)    
(Intercept)                     -13.36969    0.71720 -18.641  < 2e-16 ***
Spec.Period.L                     1.14118    0.26263   4.345 6.64e-05 ***
Spec.Period.Q                    -0.29429    0.23676  -1.243 0.219572    
Spec.Period.C                    -0.26773    0.23879  -1.121 0.267458    
Spec.Period^4                     0.77807    0.21443   3.629 0.000659 ***
DNA.fragment.length.median.oral   0.03769    0.00905   4.164 0.000120 ***"

summary(glm(ARG.rpoB ~ Spec.Period + DNA.fragment.length.median.oral, data = oral.amrfam.prop, family = quasibinomial(link = "logit"),
            weights = oral.amrfam.prop$Reads.n.oral.bacteria))
"                                  Estimate Std. Error t value Pr(>|t|)    
(Intercept)                     -12.116805   0.760018 -15.943  < 2e-16 ***
Spec.Period.L                    -0.135746   0.228205  -0.595   0.5546    
Spec.Period.Q                    -0.912240   0.211772  -4.308 7.52e-05 ***
Spec.Period.C                     0.086018   0.169641   0.507   0.6143    
Spec.Period^4                     0.070857   0.177328   0.400   0.6911    
DNA.fragment.length.median.oral   0.018613   0.009783   1.903   0.0627 .  "


##### oral bacteria ordination in bears #####
k.oral.bears <- k.oral[sample.id.byyear,]
k.oral.bears <- k.oral.bears[ ,colSums(k.oral.bears) > 0]

# Euclidean
k.oral.bears.clr <- apply(k.oral.bears, MARGIN = c(1,2), FUN = function(x) x+1)
k.oral.bears.clr <- as.data.frame(t(apply(k.oral.bears.clr, 1, function(x){log(x) - mean(log(x))})))
k.oral.bears.euc <- vegdist(k.oral.bears.clr, method = "euclidean")
k.oral.bears.euc.nmds <- metaMDS(k.oral.bears.euc, k = 2, try = 20, trymax = 150, maxit = 2000)
# solution reached run 20
k.oral.bears.euc.nmds$stress
# 0.1881088
k.oral.bears.euc.nmds.scores <- as.data.frame(scores(k.oral.bears.euc.nmds))
k.oral.bears.euc.nmds.scores$SampleID <- row.names(k.oral.bears.euc.nmds.scores)
k.oral.bears.euc.nmds.scores <- merge(k.oral.bears.euc.nmds.scores, oral.amrfam.cts, by = "SampleID")


Sfig.oral.ord.a <- ggplot(k.oral.bears.euc.nmds.scores,aes(NMDS1, NMDS2,fill=Spec.Period, shape=Spec.Period))+
  geom_vline(xintercept = 0, linetype = 3, size = 1, color = "grey42") + geom_hline(yintercept = 0, linetype = 3, size = 1, color = "grey42")+
  geom_point(size = 4)+
  scale_shape_manual(values = c(21:25), labels = c("pre 1951","1951-1970","1971-1985","1986-2000","post 2000"))+
  scale_fill_manual(values = spec.period.cols, labels = c("pre 1951","1951-1970","1971-1985","1986-2000","post 2000"))+
  scale_x_continuous(limits = c(-55,55)) + scale_y_continuous(limits = c(-55,55))+
  
  coord_fixed(ratio = 1)+
  labs(shape = "Time period", fill = "Time period")+
  theme_bw() + theme(panel.grid = element_blank(), legend.position = "right",
                     axis.text = element_text(size = 12, colour = "black"), axis.title = element_text(size = 13), panel.border = element_rect(size = 1.2))

Sfig.oral.ord.c <- ggplot(k.oral.bears.euc.nmds.scores,aes(NMDS1, NMDS2, color=Total.AMR.prop))+
  geom_vline(xintercept = 0, linetype = 3, size = 1, color = "grey42") + geom_hline(yintercept = 0, linetype = 3, size = 1, color = "grey42")+
  geom_point(size = 4)+
  scale_color_distiller(palette = "Spectral")+
  scale_x_continuous(limits = c(-55,55)) + scale_y_continuous(limits = c(-55,55))+
  coord_fixed(ratio = 1)+
  labs(color = "Total AMR load")+
  theme_bw() + theme(panel.grid = element_blank(), legend.position = "right",
                     axis.text = element_text(size = 12, colour = "black"), axis.title = element_text(size = 13), panel.border = element_rect(size = 1.2))

Sfig.oral.ord.l1 <- get_legend(
  Sfig.oral.ord.a + theme(legend.text = element_text(size = 12),
                          legend.title = element_text(size = 14, face = "bold"),
                          legend.box.margin = margin(l=5,r=5))
)
Sfig.oral.ord.l2 <- get_legend(
  Sfig.oral.ord.c + theme(legend.text = element_text(size = 12),
                          legend.title = element_text(size = 14,  face = "bold"),
                          legend.box.margin = margin(l=5,r=5))
)

adonis(k.oral.bears.euc ~ k.oral.bears.euc.nmds.scores$Spec.Period, permutations = 1000)
adonis(k.oral.bears.euc ~ k.oral.bears.euc.nmds.scores$Total.AMR.prop, permutations = 1000)

# Jaccard binary
k.oral.bears.jac <- vegdist(k.oral.bears, method = "jaccard", binary = TRUE)
k.oral.bears.jac.nmds <- metaMDS(k.oral.bears.jac, k = 2, try = 20, trymax = 150, maxit = 2000)
# solution reached run 20
k.oral.bears.jac.nmds$stress
# 0.1315104
k.oral.bears.jac.nmds.scores <- as.data.frame(scores(k.oral.bears.jac.nmds))
k.oral.bears.jac.nmds.scores$SampleID <- row.names(k.oral.bears.jac.nmds.scores)
k.oral.bears.jac.nmds.scores <- merge(k.oral.bears.jac.nmds.scores, oral.amrfam.cts, by = "SampleID")

Sfig.oral.ord.b <- ggplot(k.oral.bears.jac.nmds.scores,aes(NMDS1, NMDS2,fill=Spec.Period, shape=Spec.Period))+
  geom_vline(xintercept = 0, linetype = 3, size = 1, color = "grey42") + geom_hline(yintercept = 0, linetype = 3, size = 1, color = "grey42")+
  geom_point(size = 4)+
  scale_shape_manual(values = c(21:25), labels = c("pre 1951","1951-1970","1971-1985","1986-2000","post 2000"))+
  scale_fill_manual(values = spec.period.cols, labels = c("pre 1951","1951-1970","1971-1985","1986-2000","post 2000"))+
  scale_x_continuous(limits = c(-0.6,0.6)) + scale_y_continuous(limits = c(-0.6,0.6))+
  coord_fixed(ratio = 1)+
  labs(shape = "Time period", fill = "Time period")+
  theme_bw() + theme(panel.grid = element_blank(), legend.position = "none",
                     axis.text = element_text(size = 12, colour = "black"), axis.title = element_text(size = 13), panel.border = element_rect(size = 1.2))

Sfig.oral.ord.d <- ggplot(k.oral.bears.jac.nmds.scores,aes(NMDS1, NMDS2, color=Total.AMR.prop))+
  geom_vline(xintercept = 0, linetype = 3, size = 1, color = "grey42") + geom_hline(yintercept = 0, linetype = 3, size = 1, color = "grey42")+
  geom_point(size = 4)+
  scale_color_distiller(palette = "Spectral")+
  scale_x_continuous(limits = c(-0.6,0.6)) + scale_y_continuous(limits = c(-0.6,0.6))+
  coord_fixed(ratio = 1)+
  labs(color = "Total AMR load")+
  theme_bw() + theme(panel.grid = element_blank(), legend.position = "none",
                     axis.text = element_text(size = 12, colour = "black"), axis.title = element_text(size = 13), panel.border = element_rect(size = 1.2))

adonis(k.oral.bears.jac ~ k.oral.bears.jac.nmds.scores$Spec.Period, permutations = 1000)
adonis(k.oral.bears.jac ~ k.oral.bears.jac.nmds.scores$Total.AMR.prop, permutations = 1000)

##### oral bacterial genera correlations with ARGs #####
k.oral.genus.prop.bears <- k.oral.genus.prop[sample.id.byyear,]
k.oral.genus.prop.bears <- k.oral.genus.prop.bears[,colSums(k.oral.genus.prop.bears) > 0]

# correlations between the four most abundant ARGs and bacteria genus abundance
crepe.out <- ccrepe(x = oral.amrfam.prop[,c("ARG.ABCF","ARG.ileS","ARG.RND","ARG.rpoB")],
                    y = k.oral.genus.prop.bears,
                    sim.score = nc.score, iterations = 20, min.subj = 20, verbose = T)
crepe.out.m <- merge(melt(crepe.out$q.values, id.vars = row.names(crepe.out$q.values), value.name = "q.value"),
                     melt(crepe.out$sim.score, id.vars = row.names(crepe.out$sim.score), value.name = "sim.score"),
                     by = c("Var1","Var2"))
names(crepe.out.m)[1:2] <- c("ARG.class","Genus")
crepe.out.m[which(crepe.out.m$q.value < 0.05 & crepe.out.m$sim.score > 0),]

# Supp Figure 3 plots
k.genus.arg <- k.oral.genus.prop[,c("Neisseria","Parvimonas","Actinomyces","Gemella")]
k.genus.arg$SampleID <- row.names(k.genus.arg)
k.genus.arg <- merge(k.genus.arg, oral.amrfam.prop[,c("SampleID","ARG.ABCF","ARG.ileS","ARG.ABC","ARG.RND","ARG.rpoB","Spec.Period","DNA.fragment.length.median.oral","Reads.n.oral.bacteria")],
                     by = "SampleID")
k.genus.arg$Bracken.oral.ct <- sapply(k.genus.arg$SampleID, function(x){
  k.oral.genus[which(k.oral.genus$SampleID == x),"Bracken.oral.ct"]
})

Sfig.gen.actino <- ggplot(k.genus.arg, aes(Spec.Period, Actinomyces))+
  geom_boxplot(outlier.colour = NA)+
  geom_point(aes(fill = Spec.Period, shape = Spec.Period), size = 2, position = position_jitterdodge())+
  scale_shape_manual(values = c(21:25))+
  scale_fill_manual(values = spec.period.cols)+
  scale_x_discrete(labels = c("pre \n1951","1951-\n1970","1971-\n1985","1986-\n2000","post\n2000"))+
  labs(x = "", y = "Genus relative abundance", title = "Actinomyces")+
  theme_classic() + theme(axis.text = element_text(colour = "black", size = 10), plot.title = element_text(size = 11, face = "bold"), legend.position = "none")

Sfig.gen.parvi <- ggplot(k.genus.arg, aes(Spec.Period, Parvimonas))+
  geom_boxplot(outlier.colour = NA)+
  geom_point(aes(fill = Spec.Period, shape = Spec.Period), size = 2, position = position_jitterdodge())+
  scale_shape_manual(values = c(21:25))+
  scale_fill_manual(values = spec.period.cols)+
  scale_x_discrete(labels = c("pre \n1951","1951-\n1970","1971-\n1985","1986-\n2000","post\n2000"))+
  labs(x = "", y = "Genus relative abundance", title = "Parvimonas")+
  theme_classic() + theme(axis.text = element_text(colour = "black", size = 10), plot.title = element_text(size = 11, face = "bold"), legend.position = "none")

Sfig.gen.gamella <- ggplot(k.genus.arg, aes(Spec.Period, Gemella))+
  geom_boxplot(outlier.colour = NA)+
  geom_point(aes(fill = Spec.Period, shape = Spec.Period), size = 2, position = position_jitterdodge())+
  scale_shape_manual(values = c(21:25))+
  scale_fill_manual(values = spec.period.cols)+
  scale_x_discrete(labels = c("pre \n1951","1951-\n1970","1971-\n1985","1986-\n2000","post\n2000"))+
  labs(x = "", y = "Genus relative abundance", title = "Gemella")+
  theme_classic() + theme(axis.text = element_text(colour = "black", size = 10), plot.title = element_text(size = 11, face = "bold"), legend.position = "none")

lm.abcf <- lm(ARG.ABCF ~ Actinomyces, data=k.genus.arg)
Sfig.cor.abcf <- ggplot(k.genus.arg, aes(Actinomyces, ARG.ABCF))+
  geom_abline(slope = lm.abcf$coefficients[[2]], intercept = lm.abcf$coefficients[[1]], color = "grey66", linetype=1, size = 1)+
  geom_point(size = 2, color = "black")+
  scale_shape_manual(values = c(21:25))+
  scale_fill_manual(values = spec.period.cols)+
  labs(x = "Genus relative abundance", y = "ARG family load", title = "ABC-F vs Actinomyces")+
  theme_classic() + theme(axis.text = element_text(colour = "black", size = 10), plot.title = element_text(size = 11, face = "bold"), legend.position = "none")

lm.rpob <- lm(ARG.rpoB ~ Actinomyces, data=k.genus.arg)
Sfig.cor.rpoB <- ggplot(k.genus.arg, aes(Actinomyces, ARG.rpoB))+
  geom_abline(slope = lm.rpob$coefficients[[2]], intercept = lm.rpob$coefficients[[1]], color = "grey66", linetype=1, size = 1)+
  geom_point(size = 2, color = "black")+
  scale_shape_manual(values = c(21:25))+
  scale_fill_manual(values = spec.period.cols)+
  labs(x = "Genus relative abundance", y = "ARG family load", title = "rpoB vs Actinomyces")+
  theme_classic() + theme(axis.text = element_text(colour = "black", size = 10), plot.title = element_text(size = 11, face = "bold"), legend.position = "none")

lm.parvi <- lm(ARG.ileS ~ Parvimonas, data=k.genus.arg)
Sfig.cor.parvi <- ggplot(k.genus.arg, aes(Parvimonas, ARG.ileS))+
  geom_abline(slope = lm.parvi$coefficients[[2]], intercept = lm.parvi$coefficients[[1]], color = "grey66", linetype=1, size = 1)+
  geom_point(size = 2, color = "black")+
  scale_shape_manual(values = c(21:25))+
  scale_fill_manual(values = spec.period.cols)+
  labs(x = "Genus relative abundance", y = "ARG family load", title = "ileS vs Parvimonas")+
  theme_classic() + theme(axis.text = element_text(colour = "black", size = 10), plot.title = element_text(size = 11, face = "bold"), legend.position = "none")

lm.gem <- lm(ARG.ileS ~ Gemella, data=k.genus.arg)
Sfig.cor.gem <- ggplot(k.genus.arg, aes(Gemella, ARG.ileS))+
  geom_abline(slope = lm.gem$coefficients[[2]], intercept = lm.gem$coefficients[[1]], color = "grey66", linetype=1, size = 1)+
  geom_point(size = 2, color = "black")+
  scale_shape_manual(values = c(21:25))+
  scale_fill_manual(values = spec.period.cols)+
  labs(x = "Genus relative abundance", y = "ARG family load", title = "ileS vs Gemella")+
  theme_classic() + theme(axis.text = element_text(colour = "black", size = 10), plot.title = element_text(size = 11, face = "bold"), legend.position = "none")

# Figure 2
fig.arg.RND <- ggplot(oral.amrfam.prop, aes(Spec.Period, ARG.RND))+
  geom_boxplot(outlier.colour = NA)+
  geom_point(aes(fill = Spec.Period, shape = Spec.Period), size = 2, position = position_jitterdodge())+
  scale_shape_manual(values = c(21:25))+
  scale_fill_manual(values = spec.period.cols)+
  geom_vline(xintercept = 1.5, color = "grey50", linetype = 2)+
  geom_vline(xintercept = 2.5, color = "grey50", linetype = 2)+
  geom_vline(xintercept = 3.5, color = "grey50", linetype = 2)+
  geom_vline(xintercept = 4.5, color = "grey50", linetype = 2)+
  scale_x_discrete(labels = c("pre 1951","1951-1970","1971-1985","1986-2000","post 2000"))+
  labs(x = "Time period", y = "RND family load")+
  theme_classic() + theme(axis.text = element_text(size = 10, colour = "black"), axis.title = element_text(size = 11),
                          legend.position = "none", axis.text.x = element_text(size = 9))

fig.gen.neiss <- ggplot(k.genus.arg, aes(Spec.Period, Neisseria))+
  geom_boxplot(outlier.colour = NA)+
  geom_point(aes(fill = Spec.Period, shape = Spec.Period), size = 2, position = position_jitterdodge())+
  scale_shape_manual(values = c(21:25))+
  scale_fill_manual(values = spec.period.cols)+
  geom_vline(xintercept = 1.5, color = "grey50", linetype = 2)+
  geom_vline(xintercept = 2.5, color = "grey50", linetype = 2)+
  geom_vline(xintercept = 3.5, color = "grey50", linetype = 2)+
  geom_vline(xintercept = 4.5, color = "grey50", linetype = 2)+
  scale_x_discrete(labels = c("pre 1951","1951-1970","1971-1985","1986-2000","post 2000"))+
  labs(x = "Time period", y = "Neisseria relative abundance")+
  theme_classic() + theme(axis.text = element_text(size = 10, colour = "black"), axis.title = element_text(size = 11),
                          legend.position = "none", axis.text.x = element_text(size = 10))

lm.rnd <- lm(ARG.RND ~ Neisseria, data=k.genus.arg)
fig.cor.rnd <- ggplot(k.genus.arg, aes(Neisseria, ARG.RND))+
  geom_abline(slope = lm.rnd$coefficients[[2]], intercept = lm.rnd$coefficients[[1]], color = "grey66", linetype=1, size = 1)+
  geom_point(size = 2, color = "black")+
  scale_shape_manual(values = c(21:25))+
  scale_fill_manual(values = spec.period.cols)+
  labs(x = "Neisseria relative abundance", y = "RND family load")+
  theme_classic() + theme(axis.text = element_text(size = 10, colour = "black"), axis.title = element_text(size = 11),
                          legend.position = "none", axis.text.x = element_text(size = 10))

fig3.p1 <- plot_grid(fig.arg.RND, fig.gen.neiss, fig.cor.rnd,
                     nrow = 1, align = "hv", axis = "tblr",
                     rel_widths = c(1,1,0.75))

fig3.leg <- plot_grid(Sfig.oral.ord.l1, 
                      Sfig.oral.ord.l2, 
                      nrow = 2, align = "hv", axis = "tblr")

tiff("figures/Figure_2_ARG_family_and_oral_bacteria.tiff", units = "in", width = 12, height = 8.5, res = 600, compression = "zip")
plot_grid(fig.arg.RND, fig.gen.neiss, fig.cor.rnd,
          Sfig.oral.ord.a + theme(legend.position = "none", axis.text = element_text(size = 10, colour = "black"), axis.title = element_text(size = 11), panel.border = element_rect()),
          Sfig.oral.ord.c + theme(legend.position = "none", axis.text = element_text(size = 10, colour = "black"), axis.title = element_text(size = 11), panel.border = element_rect()),
          fig3.leg,
          nrow = 2, rel_widths = c(1,1,0.75,1,1,0.75),
          align = "hv", axis = "tblr")
dev.off()

# Supp Figure 3
Sfig.arg.all <- plot_grid(Sfig.arg.ABC, Sfig.arg.lipA, Sfig.arg.parY,
                          Sfig.arg.abcf, Sfig.arg.rpoB, Sfig.arg.ileS,
                          Sfig.gen.actino, Sfig.gen.parvi, Sfig.gen.gamella,
                          labels = c("A","B","C","D","E","F","G","H","I"),
                          ncol = 3, align = "hv", axis = "tblr")

Sfig.cor.all <- plot_grid(Sfig.cor.abcf, Sfig.cor.rpoB, Sfig.cor.parvi, Sfig.cor.gem,
                          labels = c("J","K","L","M"), nrow = 1, align = "hv", axis = "tblr")

Sfig.oral.ord.presence <- plot_grid((Sfig.oral.ord.b + theme(legend.position = "none")), Sfig.oral.ord.l1, NULL,
                                    (Sfig.oral.ord.d + theme(legend.position = "none")), Sfig.oral.ord.l2, NULL,
                                    labels = c("O","","","P","",""),
                                    nrow = 1, rel_widths = c(1,0.5,0.5,1,0.5,0.5))

tiff("figures/Supp_Figure_S3_ARG_family_and_oral_genera.tiff", units = "in", width = 12, height = 12, res = 600, compression = "zip")
plot_grid(Sfig.arg.all, Sfig.cor.all, Sfig.oral.ord.presence,
          nrow = 3, align = "hv", axis = "tblr", rel_heights = c(3,1,1.5))
dev.off()

##### END #####
