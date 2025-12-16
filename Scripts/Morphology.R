setwd("/Users/ianbrennan/Documents/GitHub/Asterophryinae")

adata <- read.csv("/Users/ianbrennan/Documents/GitHub/Asterophryinae/Data/Asterophryinae_RAWdata.csv"); nrow(adata)
# filter out non-adult individuals
adata <- adata[adata$Sex%in%c("F","M","","?"),]; nrow(adata)
# select traits of interest
adata <- dplyr::select(adata, Genus.species, SVL,HLL,HW,EY,IN,EN,THIRD); nrow(adata)
adata <- adata[complete.cases(adata),]; nrow(adata); paste(length(unique(adata$Genus.species)),": species")
adata$Genus <- strex::str_before_first(adata$Genus.species," ")

# how many samples do we have per genus
table(adata$Genus)

# how many species do we have per genus?
table(strex::str_before_first(unique(adata$Genus.species)," ")) 

############################################################################

# Inspect the raw data without log transformation or size correction

pca <- prcomp(adata[,2:8]); summary(pca); biplot(pca,cex=0)
pca <- data.frame(pca$x)
pca$Genus <- adata$Genus

ggplot() +
  geom_point(data=pca,aes(x=PC1,y=PC2,fill=Genus),shape=21,size=5) +
  theme_bw()

ggplot() +
  geom_point(data = transform(pca, Genus = NULL), color = "grey85", aes(x=PC1, y=PC2),size=3,alpha=0.5) +
  geom_point(data=pca,aes(x=PC1,y=PC2,fill=Genus),shape=21,size=3) +
  facet_wrap(~Genus) + theme_bw() + theme(legend.position = "none")

############################################################################

############################################################################

# Do size correction (LSR) and log transform data for ALL SAMPLES

# make sure your variables are all numeric
all.raw <- adata
all.raw$Genus <- strex::str_before_first(all.raw$Genus.species," ")

# Compute the geometric mean for obtaining size
sizes <- apply(all.raw[2:8], 1, prod)^(1/ncol(all.raw[2:8]))

# Compute the log shape ratios
sampLSR <- all.raw[2:8]/sizes
sampLSR$Size <- sizes
sampLSR <- log(sampLSR)

# Add the taxon, group, and family, names back onto the data frame
sampLSR$Genus.species <- all.raw$Genus.species;
sampLSR$Genus <- all.raw$Genus

# Save the log shape ratio file
# write.csv(allLSR, row.names=FALSE,
#           file="Data/Amphibolurinae_AllLSR.csv")

############################################################################

pca.samp <- prcomp(sampLSR[,1:8]); summary(pca.samp); biplot(pca.samp,cex=0)
pcaLSR.samp <- data.frame(pca.samp$x)
pcaLSR.samp$Genus <- sampLSR$Genus
pcaLSR.samp$Genus.species <- sampLSR$Genus.species

genera <- unique(sampLSR$Genus)
shapes <- NULL

des.shape <- function(x){
  if(x %in% c(genera[c(1,5,9,13,17)])){return(21)}
  if(x %in% c(genera[c(2,6,10,14)])){return(22)}
  if(x %in% c(genera[c(3,7,11,15)])){return(23)}
  if(x %in% c(genera[c(4,8,12,16)])){return(24)}
}
pcaLSR.samp$shape <- sapply(pcaLSR.samp$Genus, des.shape)

ggplot() +
  geom_point(data=pcaLSR.samp,aes(x=PC1,y=PC2,fill=Genus),shape=pcaLSR.samp$shape,size=3) +
  theme_bw()

ggplot() +
  geom_point(data=pcaLSR.samp,aes(x=PC3,y=PC4,fill=Genus),shape=21,size=3) +
  theme_bw()

pca.wrap <- ggplot() +
  geom_point(data = transform(pcaLSR.samp, Genus = NULL), color = "grey85", aes(x=PC1, y=PC2),size=3,alpha=0.5) +
  geom_point(data=pcaLSR.samp,aes(x=PC1,y=PC2,fill=Genus),shape=21,size=3) +
  facet_wrap(~Genus) + theme_bw() + theme(legend.position="none")

pc3.pc4 <- ggplot() +
  geom_point(data = transform(pcaLSR.samp, Genus = NULL), color = "grey85", aes(x=PC3, y=PC4),size=3,alpha=0.5) +
  geom_point(data=pcaLSR.samp,aes(x=PC3,y=PC4,fill=Genus),shape=21,size=3) +
  facet_wrap(~Genus) + theme_bw() + theme(legend.position="none")

library(patchwork)

pc1.pc2 + pc3.pc4

############################################################################
############################################################################
############################################################################

# Translate data to species means and do size correction and log transformation

# Create a tibble to get the species means  
sp.means <- adata %>%
  group_by(Genus.species) %>%
  summarise_at(vars(SVL,HLL,HW,EY,IN,EN,THIRD), mean)
sp.means$Genus <- strex::str_before_first(sp.means$Genus.species," ")

# Save the species mean file
write.csv(sp.means, row.names=FALSE,
          file="Data/Amphibolurinae_Morphology_spMEANS.csv")

############################################################################

# make sure your variables are all numeric
all.raw.traits <- sp.means

# Compute the geometric mean for obtaining size
allsize <- apply(sp.means[2:8], 1, prod)^(1/ncol(sp.means[2:8]))

# Compute the log shape ratios
allLSR <- sp.means[2:8]/allsize
allLSR$Size <- allsize
allLSR <- log(allLSR)

# Add the taxon, group, and family, names back onto the data frame
rownames(allLSR) <- sp.means$Genus.species;
allLSR$Genus <- sp.means$Genus

# Save the log shape ratio file
write.csv(allLSR, row.names=FALSE,
          file="Data/Amphibolurinae_AllLSR.csv")

############################################################################

pca.res <- prcomp(allLSR[,1:8]); summary(pca.res); biplot(pca.res,cex=0)
pcaLSR <- data.frame(pca.res$x)
pcaLSR$Genus <- strex::str_before_first(rownames(pcaLSR)," ")

ggplot() +
  geom_point(data=pcaLSR,aes(x=PC1,y=PC2,fill=Genus),shape=21,size=5) +
  theme_bw()

ggplot() +
  geom_point(data = transform(pcaLSR, Genus = NULL), color = "grey85", aes(x=PC1, y=PC2),size=3,alpha=0.5) +
  geom_point(data=pcaLSR,aes(x=PC1,y=PC2,fill=Genus),shape=21,size=3) +
  facet_wrap(~Genus) + theme_bw()

############################################################################

############################################################################

############################################################################

# if you want to use the uncorrected measurements
rf.data <- log(adata[,2:8]); 
rf.data$Genus <- as.factor(adata$Genus)


############################################################################

library(randomForest)

# set a seed to make this reproducible
set.seed(2025)

# run the RandomForest and visualize it
forest <- randomForest(Genus ~ ., data = rf.data, ntree = 10000, localImp = TRUE)
plot(forest)

# extract the categorizations for each taxon, map to regime designation
rf.votes <- data.frame(forest$votes)
rf.votes$Genus <- rf.data$Genus
rf.votes$Genus.species <- adata$Genus.species

# round the values
rf.votes[,1:17] <- sapply(rf.votes[,1:17],function(x) round(x,2))

# the result is categorizations without correcting for size

write.csv(forest$confusion, file="Data/RandomForests/RF_Confusion_RAW.csv")
write.csv(rf.votes, file="Data/RandomForests/RF_Results_RAW.csv")

############################################################################

# set a seed to make this reproducible
set.seed(2025)

rf.data2 <- sampLSR[,c(1:8,10)]
rf.data2$Genus <- as.factor(rf.data2$Genus)

# run the RandomForest and visualize it
forest2 <- randomForest(Genus ~ ., data = rf.data2, ntree = 10000, localImp = TRUE, mtry = 4)
plot(forest2)

# extract the categorizations for each taxon, map to regime designation
rf.votes2 <- data.frame(forest2$votes)
rf.votes2$Genus <- rf.data2$Genus
rf.votes2$Genus.species <- sampLSR$Genus.species

# round the values
rf.votes2[,1:17] <- sapply(rf.votes2[,1:17],function(x) round(x,2))

# the result is categorizations after correcting for size

write.csv(forest2$confusion, file="Data/RandomForests/RF_Confusion_LSR.csv")
write.csv(rf.votes2, file="Data/RandomForests/RF_Results_LSR.csv")



############################################################################

# 

# LDA instructions: https://eol.pages.cms.hu-berlin.de/gcg_quantitative-methods/Lab11_LDA_Model-assessment.html#:~:text=assumptions%20are%20violated.-,Classification%20error,as%20opposed%20to%20continuous%20variables.

options(repos = c(
  fawda123 = 'https://fawda123.r-universe.dev',
  CRAN = 'https://cloud.r-project.org'))
install.packages("ggord")

# TESTING OUT LDA (VERY DIFFERENT RESULTS!)
lda.test <- MASS::lda(Genus~., sampLSR[,c(1:3,10)])
ggord::ggord(lda.test, lda.test$Genus)
plot.rotations(pca.obj = lda.test, which.pca = "lda")

lda.sp <- MASS::lda(Genus_species~.,sampleLSR[,c(1:20)])
ggord::ggord(lda.sp, lda.sp$Genus)
plot.rotations(pca.obj = lda.sp, which.pca = "lda")
lda.sp$scores <- data.frame(predict(lda.sp)$x)
lda.sp$scores$Genus <- sampleLSR$Genus
lda.sp$scores$Genus_species <- sampleLSR$Genus_species

lda.gen <- MASS::lda(Genus~., sampLSR[,c(1:6,8,10)])
ggord::ggord(lda.gen, lda.gen$Genus)
#plot.rotations(pca.obj = lda.sp, which.pca = "lda")
lda.gen$scores <- data.frame(predict(lda.gen)$x)
lda.gen$scores$Genus <- sampLSR$Genus
lda.gen$scores$Genus_species <- sampLSR$Genus_species

lda.pred <- predict(lda.gen)
conf <- table(list(predicted=lda.pred$class, observed=sampLSR$Genus))

# Precision (Positive predicted value)
diag(conf) / rowSums(conf)

# Sensitivity
diag(conf) / colSums(conf)

lda.conf <- caret::confusionMatrix(conf)
write.csv(lda.conf$table, file="Data/Classification_Analyses/LDA_Confusion_LSR.csv",quote=F)


ggplot() +
  geom_point(data = transform(lda.gen$scores, Genus = NULL), 
             colour = "grey85", aes(x=LD1, y=LD2), size=2) +
  #geom_polygon(data=x.hull, aes(x=PC1, y=PC2, fill=genus), alpha=0.5) +
  geom_point(data=lda.gen$scores, aes(x=LD1, y=LD2, fill=Genus), size=3, pch=21) +
  theme_classic() + theme(legend.position = "none") + facet_wrap(~Genus)

lda.tpd <- TPD::TPDs(species = lda.gen$scores$Genus, lda.gen$scores[,1:2], alpha=0.999)
TPD::plotTPD(TPD = lda.tpd)
fun.lda <- funspace::funspace(x=lda.tpd, PCs=c(1,2), n_divisions=100)


lda.res <- data.frame(lda.pred$x); lda.res$Genus <- sampLSR$Genus
lda.wrap <- ggplot() +
  geom_point(data = transform(lda.res, Genus = NULL), color = "grey85", aes(x=LD1, y=LD2),size=3,alpha=0.5) +
  geom_point(data=lda.res,aes(x=LD1,y=LD2,fill=Genus),shape=21,size=3) +
  facet_wrap(~Genus) + theme_bw() + theme(legend.position="none")

########################


install.packages("mda")
library(mda)

fda.gen <- mda::fda(Genus~., sampLSR[,c(1:6,8,10)])
write.csv(fda.gen$confusion, file="Data/Classification_Analyses/MDA_Confusion_LSR.csv",quote=F)

fda.pred <- data.frame(plot.fda2(fda.gen))
colnames(fda.pred) <- c("FD1","FD2","FD3","FD4","FD5","FD6","FD7","FD8","FD9")
fda.pred$Genus <- sampLSR$Genus

fda.wrap <- ggplot() +
  geom_point(data = transform(fda.pred, Genus = NULL), color = "grey85", aes(x=FD1, y=FD2),size=3,alpha=0.5) +
  geom_point(data=fda.pred,aes(x=FD1,y=FD2,fill=Genus),shape=21,size=3) +
  facet_wrap(~Genus) + theme_bw() + theme(legend.position="none")


############################################################################

library(patchwork)

pca.wrap /
  lda.wrap /
  fda.wrap

############################################################################

#######################

class.error <- read.csv("Data/Classification_Analyses/Confusion_All.csv")
ggplot() +
  geom_jitter(data=class.error, aes(x=genus,y=class.error,fill=method),shape=21,size=3,width=0.5) +
  scale_x_discrete(limits=rev) + coord_flip() + theme_classic() + theme(legend.position="bottom")


###

plot.fda2 <- function (x, data, coords = c(1, 2), group = c("true","predicted"), colors, pch, mcolors = colors, mpch, pcex = 0.5, mcex = 2.5, ...){
    object <- x                         # generic/method
    group <- match.arg(group)
    if (missing(data)) {
      vars <- predict(object, type = "var")
      g <- predict(object)
      group <- "predict"
    }
    else {
      if(group=="predicted"){
        vars <- predict(object, data, type = "var")
        g <- predict(object, data)
      }
      else{
        ff <- terms(object)
        attr(ff, "intercept") <- 0
        m <- model.frame(ff, data)
        x <- model.matrix(ff, m)
        vars <- predict(object, x, type = "var")
        g <- model.extract(m, "response")
      }
    }
    means <- object$means
    if(ncol(means)==1)stop("Only one canonical variate; plot requires at least two")
    g <- as.factor(g)
    cc <- as.numeric(g)
    np <- seq(levels(g))
    mit.colors=c(Orange = "#FF9233", Cyan = "#29D0D0", Lt.Green = "#81C57A",
                 Dk.Gray = "#575757", Red = "#AD2323", Blue = "#2A4BD7", Green = "#1D6914",
                 Brown = "#814A19", Purple = "#8126C0", Lt.Gray = "#A0A0A0", Yellow = "#FFEE33",
                 Pink = "#FFCDF3")
    if (missing(colors)) colors <- mit.colors
    colors <- rep(colors, length = length(np))
    if (missing(pch))
      pch <- paste(np)
    else pch <- rep(paste(pch), length = length(np))
    mcolors <- rep(mcolors, length = length(np))
    if (missing(mpch))
      mpch <- pch
    else mpch <- rep(paste(mpch), length = length(np))
    assign <- object$assign
    if (is.null(assign))
      assign <- split(seq(np), seq(np))
    if (!is.matrix(coords)) {
      coords <- matrix(coords, length(coords), length(coords))
      tt <- lower.tri(coords)
      coords <- cbind(t(coords)[tt], coords[tt])
    }
    for (ii in seq(nrow(coords))) {
      coord.pair <- coords[ii, ]
      plot(rbind(vars[, coord.pair], means[, coord.pair]),
           ..., type = "n", xlab = paste("Canonical Var",
                                         coord.pair[1]), ylab = paste("Canonical Var",
                                                                      coord.pair[2]), main = paste("Discriminant Plot for",
                                                                                                   group, "classes"))
      for (i in np) {
        which <- cc == i
        if (any(which))
          points(vars[which, coord.pair, drop = FALSE], col = colors[i],
                 pch = pch[i], cex = pcex)
        points(means[assign[[i]], coord.pair, drop = FALSE],
               col = mcolors[i], pch = 1, cex = mcex)
        points(means[assign[[i]], coord.pair, drop = FALSE],
               col = mcolors[i], pch = mpch[i], cex = mcex/2)
      }
    }
    invisible()
    return(vars)
  }
testo <- plot.fda2(fda.gen)
