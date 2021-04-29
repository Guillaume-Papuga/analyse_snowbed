#######################################################
# Project : Analysis of snowbed vegetation dynamics
# Script : 10.community_analysis
# Analysis of community structure based solely on plant composition
# Authors : Guillaume Papuga & Thomas Masclaux
# Last update : 28 april 2021
#######################################################

#3) Structure spatiale-----



##PCOA avec BRAY distance
df_pcoa = df ### df pour l'analyse PCOA

df_pcoa$col1<- paste(df_pcoa$site, 
                     df_pcoa$placette,
                     df_pcoa$year, sep="_") # classe site*year*placette, cette colonne servira ? ?tiqueter, les lignees dans la construction du graph

rownames(df_pcoa) = make.names(paste(df_pcoa$site,df_pcoa$placette,df_pcoa$year, sep="_"), unique=TRUE) # on ?tiquette les lignes de la table

####cr?ation du mod?le d'analyse
dist= distance(df_pcoa[,4:99], "bray-curtis") # cr?ation de la matrice
PCoA=pcoa(dist)
PCoA$values/sum(PCoA$values)

barplot(PCoA$values$Relative_eig[1:10]) ### on visualise la part explicative des composantes, on choisit les deux premiers

head(PCoA)

##choix des axes et extraction des valeurs
PC1 <- PCoA$vectors[,1] #extraction des coordonn?es de points axe 1
PC2 <- PCoA$vectors[,2] #extraction des coordonn?es de points axe 2

PCs <- data.frame(cbind(PC1,PC2)) # coordonn?es des points en 1 tableau

PCs$placette <- sapply((as.character(df_pcoa$placette)), "[[", 1 ) ## ?tiquette site
PCs$site <- sapply((as.character(df_pcoa$site)), "[[", 1 ) ## ?tiquette placette
PCs$df_pcoa <- sapply(as.character(df_pcoa$col1), "[[", 1 ) ### ?tiquette points

PCs$cols <- paste(df_pcoa$site,df_pcoa$placette, sep="_") ## ?tiquette pour cr?er les ?llipses

rownames(PCs)= rownames(df_pcoa)

centroids <- aggregate(cbind(PC1,PC2)~cols,PCs,mean)

fac <-str_split_fixed(centroids$cols, "_", 2)

centroid <- cbind.data.frame(fac,centroids)

colnames(centroid$`1`) = centroid$site

## Mise en place de l'esth?tique en amont
theme<-theme(panel.background = element_blank(),
             panel.border=element_rect(fill=NA),
             panel.grid.major = element_blank(),
             panel.grid.minor = element_blank(),
             strip.background=element_blank(),
             axis.text.x = element_blank(),
             axis.text.y=element_blank(),
             axis.ticks=element_blank(),
             axis.title.y = element_text(size = 16),
             axis.title.x = element_text(size = 16),
             plot.margin=unit(c(1,1,1,1),"line"),
             legend.title=element_text(size=20) , 
             legend.text=element_text(size=20)) 

pourcentage = round((PCoA$values$Relative_eig)*100, 1)
pourcentage <- paste(colnames(PCs), "(", paste( as.character(pourcentage), "%", ")", sep=""))

gg = geom_point(data =centroid, aes(colour =centroid$`1`), alpha = 1) 

gg3=geom_text(data=centroid,aes(label=`2`), size = 5)

##construction du plot  
PCOA_PLOT <- ggplot(PCs,aes(x=PC1,y=PC2), colour = site) + geom_hline(yintercept = 0, linetype = "dotted")+
  geom_vline(xintercept=0, linetype="dotted") +
  geom_point(aes(x=PC1,y=PC2, colour = site), size = 1, alpha= 0.5) +
  theme +
  stat_ellipse(aes(x=PC1, y=PC2, colour = site, group = cols),
               geom = "path", level=0.9, size = 0.75) + 
  gg + gg3 + 
  xlab(pourcentage[1]) + ylab(pourcentage[2])

PCOA_PLOT #####affichage du graphique

adonis(formula = dist~df$site*df$placette) ## analyse de variance appliqu?e ? la matrice


####HEATMAP

rownom<- paste(df_pcoa$site,df_pcoa$placette,df_pcoa$year, sep="_")
rownames(df_pcoa) = make.names(rownom, unique=TRUE)

fac = df_pcoa[,1:3]
rownames(fac) = make.names(rownom, unique=TRUE)
fac=as.data.frame(fac)
str(fac)


df_pcoa = t(df_pcoa[,4:100])
colnames(df_pcoa)=make.names(rownom, unique=TRUE)


facDS = select(fac,site,placette)
str(facDS)
rownames(facDS) <- colnames(df_pcoa)
colnames(df_pcoa)<- rownames(facDS)

# color per level
facDS$site <- factor(facDS$site, levels = c("ARB", "CASE", "CAT", "CRE", "PLA", "POR", "RAT", "ULL"))
sitecol <- c("#FF3333", "#FFCC00", "#669900", "#00CC66", "#009999", "#0066CC", "#CC66FF", "white" )
names(sitecol) <- levels(facDS$site)

facDS$placette <- factor(facDS$placette, levels = c("1", "2", "3"))
placol <- c("#CCCCCC", "#666666", "#000000")
names(placol) = levels(facDS$placette)

# Add to a list, where names match those in factors dataframe
Colour <- list(
  site = sitecol, placette = placol)

df_pcoa = df 

dat.mva = mvabund(df_pcoa[,4:99])
plot(dat.mva)

dat.nb <- manyglm(dat.mva ~ site*placette, data = fac)
dat.aov <- anova.manyglm(dat.nb, p.uni = "unadjusted", nBoot = 500) ### record de 8 minutes de chargement
dat.aov$uni.p
summary(dat.aov)
SpeciesDiffs <- which(dat.aov$uni.p["site",] < 0.05 & dat.aov$uni.p["site:placette",] > 0.05)
xx = df_pcoa[,SpeciesDiffs ]
xx = decostand(xx[,4:73], method="hellinger")
tt = xx
rownames(tt) = rownom
tt = t(tt)


# site distance matrix (Bray-Curtis)
siteDist <- vegdist(t(tt), "bray", diag = T)
# taxa distance matrix (Bray-Curtis)
taxaDist <- vegdist(tt, "bray", diag = T)

pp = pheatmap(tt, annotation_col = facDS, 
              annotation_colors = Colour,
              clustering_distance_rows = taxaDist,
              clustering_distance_cols = siteDist,
              clustering_method = 'ward.D',
              fontsize = 8, cutree_cols = 6, cutree_rows = 4, 
              show_colnames = F)









