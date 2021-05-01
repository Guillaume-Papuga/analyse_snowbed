#######################################################
# Project : Analysis of snowbed vegetation dynamics
# Script : 10.community_analysis
# Analysis of community structure based solely on plant composition
# Authors : Guillaume Papuga & Thomas Masclaux
# Last update : 28 april 2021
#######################################################

### Load data 
df = read.csv(here::here ("data", "processed", "sb_data_cast.csv"), 
              head = T, sep = ",", dec = ".")


######################### A. Global analysis of the complete dataset ######################################
### create the table that summarize the raw analysis
dist = ecodist::distance(df[,-(1:3)], "bray-curtis") # create a distance matrix
coord.an = pcoa(dist)
round(coord.an$values, 3)

barplot(coord.an$values$Relative_eig[1:10]) ### on visualise la part explicative des composantes, on choisit les deux premiers

### create a table to plot results (per plot / per year)
df_pcoa = df [,1:3] %>%
  mutate (pc1 = coord.an$vectors[,1],
          pc2 = coord.an$vectors[,2])

### create the dataset for centroids
centroids_pcoa = df_pcoa %>%
  group_by(combe, placette) %>%
  dplyr::summarize(centro_x = mean(pc1), 
            centro_y = mean(pc2))

### create the graph
# create the aesthetic
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

# compute the inertia of each axes
percentage_1 = paste("Axis 1 (", round((coord.an$values$Relative_eig)*100, 1)[1], "%)")
percentage_2 = paste("Axis 2 (", round((coord.an$values$Relative_eig)*100, 1)[2], "%)")

# create the plot 
pcoa_plot <- ggplot(df_pcoa, aes(x = pc1, y = pc2, colour = combe)) +  # classic plot configuration
  geom_point(aes(colour = combe), size = 1, alpha= 0.5) + # add points
  geom_hline(yintercept = 0, linetype = "dotted") + # add the y axis
  geom_vline(xintercept=0, linetype="dotted") + # add the x axis
  stat_ellipse(aes(colour = combe, group = interaction(combe,placette,drop=F,sep='-')), size = 0.5) + # add the ellipses
  xlab(percentage_1) + ylab(percentage_2) + # add the inertia
  geom_text(data=centroids_pcoa, aes(x = centro_x, y = centro_y, label= placette), size = 4, colour = "black") +
  theme

pcoa_plot # show the plot

pdf(here::here("outputs", "figures", "appendix_total.community.analysis.pdf"))
print(pcoa_plot)
dev.off()

######################### B. Exclude CRE1 and run the analysis again ######################################
### create the table that summarize the raw analysis
# first we delete the placette n°1 from cre
df = df %>%
  mutate(name.pla = interaction(combe, placette)) %>%
  filter (name.pla != "cre.1") %>%
  select (-name.pla)

dist = ecodist::distance(df[,-(1:3)], "bray-curtis") # create a distance matrix
coord.an = pcoa(dist)
round(coord.an$values, 3)

barplot(coord.an$values$Relative_eig[1:10]) ### on visualise la part explicative des composantes, on choisit les deux premiers

### create a table to plot results (per plot / per year)
df_pcoa = df [,1:3] %>%
  mutate (pc1 = coord.an$vectors[,1],
          pc2 = coord.an$vectors[,2])

### create the dataset for centroids
centroids_pcoa = df_pcoa %>%
  group_by(combe, placette) %>%
  dplyr::summarize(centro_x = mean(pc1), 
                   centro_y = mean(pc2))

### create the graph
# create the aesthetic
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

# compute the inertia of each axes
percentage_1 = paste("Axis 1 (", round((coord.an$values$Relative_eig)*100, 1)[1], "%)")
percentage_2 = paste("Axis 2 (", round((coord.an$values$Relative_eig)*100, 1)[2], "%)")

# create the plot 
pcoa_plot <- ggplot(df_pcoa, aes(x = pc1, y = pc2, colour = combe)) +  # classic plot configuration
  geom_point(aes(colour = combe), size = 1, alpha= 0.5) + # add points
  geom_hline(yintercept = 0, linetype = "dotted") + # add the y axis
  geom_vline(xintercept=0, linetype="dotted") + # add the x axis
  stat_ellipse(aes(colour = combe, group = interaction(combe,placette,drop=F,sep='-')), size = 0.5) + # add the ellipses
  xlab(percentage_1) + ylab(percentage_2) + # add the inertia
  geom_text(data=centroids_pcoa, aes(x = centro_x, y = centro_y, label= placette), size = 4, colour = "black") +
  theme

pcoa_plot # show the plot


######################### C. Statistical analysis without CRE1 ############################################

### Statistical analysis
adonis(formula = dist~df_pcoa$combe + df_pcoa$placette %in% df_pcoa$combe) ## multivariate analysis of variance (nested factors)


######################### D. Heatmap graphical representation #############################################

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









