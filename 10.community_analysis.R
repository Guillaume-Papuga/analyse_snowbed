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
### Define the data
data.hm = df
rownom = paste(data.hm$combe,data.hm$placette,data.hm$year, sep="_") # unique identifier of each line
rownames(data.hm) = make.names(rownom, unique=TRUE)



### Create the factor to categorize data 
fac = data.hm[,1:2]
fac$placette <- factor(fac$placette, levels = c("1", "2", "3")) # transform placette into factor

### Create the two color vectors
# For sites (combe)
combecol <- c("#FF3333", "#FFCC00", "#669900", "#00CC66", "#009999", "#0066CC", "#CC66FF", "white" )
names(combecol) <- levels(fac$combe)

# For placettes
placol <- c("#CCCCCC", "#666666", "#000000")
names(placol) = levels(fac$placette)

# Create a list
colour.hm <- list(
  combe = combecol, placette = placol)

### Select the species that contribute to the clustering pattern (procedure based on mvabund method)
# create a mvabund object
dat.mva = mvabund(data.hm[,4:ncol(data.hm)]) 
plot(dat.mva) # plot it

# Multivariate model
dat.nb = manyglm(dat.mva ~ combe + placette%in%combe, data = data.hm) # define the model (nested variables)
dat.aov = anova.manyglm(dat.nb, p.uni = "unadjusted", nBoot = 500) # compute a significance test of variables
species.dif = which(dat.aov$uni.p["combe",] < 0.05 & # select species that contribute to inter-site difference...
                    dat.aov$uni.p["combe:placette",] > 0.05) # but show no intra-site pattern (controled for)

### Extract a subset of the dataset that matches the selection
data.hm.sign.sp = data.hm[,species.dif] # select the corresponding species in the original dataset
data.hm.sign.sp.h = decostand(data.hm.sign.sp[,-c(1,2)], method="hellinger")  # transform the dataset with Hellinger distance

### Compute distance matrices (Bray-Curtis) [for side dendrograms]
combe.dist = vegdist(data.hm.sign.sp.h, "bray", diag = T) # between-sites (combe) distance
taxa.dist = vegdist(t(data.hm.sign.sp.h), "bray", diag = T) # between-variables (species) distance

### Create the plot
hm.plot = pheatmap(t(data.hm.sign.sp.h), # dataset
              annotation_col = fac, 
              annotation_colors = colour.hm, # list of colours
              clustering_distance_rows = taxa.dist, # distance matrix btw species
              clustering_distance_cols = combe.dist, # distance matrix btw sites
              clustering_method = 'ward.D', 
              fontsize = 8, cutree_cols = 6, cutree_rows = 4, 
              show_colnames = F)


