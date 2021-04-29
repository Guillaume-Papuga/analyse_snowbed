#######################################################
# Project : Analysis of snowbed vegetation dynamics
# Script : 20.community_ecology
# Analysis of environmental drivers of plant community structure
# Authors : Guillaume Papuga & Thomas Masclaux
# Last update : 28 april 2021
#######################################################

#4) Community ecology----

species = df %>% filter(grepl('2019', year)) # on s?lectionne la flore 2019
col1 = rownames(species)
species = species[, colSums(species != 0) > 0]
species = species[ rowSums(species != 0) > 0,]
fac = select_if(species, is.factor)
species = species[,4:92]
species = decostand(species,method = "hellinger")     


clim = read.table("C:/Users/FNAC/Desktop/code/Clim_data.csv", # ouverture du JDD climat
                  sep=";",stringsAsFactors=F, header = T)

clim = clim %>% select(c("site","min_summer", "max_spring", "pet", "GDD", "ppt", "RAD_summer", "max_summer", "ppt_autumn", "w_autumn")) #s?lection des variables p?ralablement observ?es sur ACP
clim$site = as.factor(clim$site)

env = left_join(fac,clim, by = "site")

#mod?le de la RDA
rdaclim = rda(species~ min_summer +max_summer+  pet + ppt+RAD_summer, env, dist="bray")
vif.cca(rdaclim) #VIF des variables explicatives, <10 de pr?f?rence

#test stat de l'analyse 
anova.cca(rdaclim, step=5000, by="axis", permutations = 500)

#test stat pour chaque variable environnemental
anova(rdaclim, by="terms", permutations = 500)

#extract of data to ggplot figure
smry = summary(rdaclim)

eig_rdacomp = round(smry$cont$importance[,1:3]*100,2)
pourcentage2 <- paste(colnames(eig_rdacomp), "(", paste( as.character(eig_rdacomp[2,]), "%", ")", sep=""))


coef(rdaclim) # coefficient de r?ajustement

rdaclim$colsum

scrs = scores(rdaclim) # score des sites, esp?ces ...
df1 = data.frame(smry$sites[,1:2]) # s?lection des valeurs par placette
df1$site = col1
rownames(df1) = make.names(paste(fac$site,fac$placette,fac$year,sep="_"))

df1$site = fac$site
df1$placette = fac$placette

df2 = data.frame(smry$biplot[,1:2])
df3 = data.frame(smry$species[,1:2])


####class species by categories and distance contribution from 
df3$arrow_lenght = sqrt((df3$RDA1)^2+(df3$RDA2^2))
df3$species = rownames(df3)
df3 = df3 %>% arrange(desc(arrow_lenght))

df3arrow = filter(df3, arrow_lenght > 0.15) # s?lection des esp?ces repr?sentatives, distance >0.15
df3arrow = df3arrow[-3,]
rownames(df3arrow) = df3arrow$species

df33 = df3arrow # DF projections des esp?ces

#mise en place de l'esth?tique
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

gg = geom_segment(data=df2, aes(x=0, xend = RDA1, y=0, yend=RDA2), 
                  color ="grey50", arrow=arrow(length = unit(0.01, "npc")))

text = geom_text(data=df2, aes(x=RDA1, y=RDA2, label = rownames(df2),
                               hjust=0.5*(1-sign(RDA1)),
                               vjust=0.5*(1-sign(RDA2))),
                 color="black", size = 5)

sp = geom_segment(data = df33,aes(x = 0, xend = RDA1, y = 0, yend = RDA2), 
                  arrow = arrow(length = unit(0.01,"npc")),linetype=1, 
                  size=0.1,colour = "red") 

text_sp =  geom_text(data = df33,aes(x=RDA1,y=RDA2,label=row.names(df33), 
                                     hjust=0.5*(1-sign(RDA1)),
                                     vjust=0.5*(1-sign(RDA2))), 
                     colour = "grey")

dbplot = ggplot(df1,aes(x=RDA1, y=RDA2)) + gg  + sp + text_sp + 
  geom_point(aes(colour = site, shape = placette), size=3) +  
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_vline(xintercept=0, linetype="dotted") + 
  theme+ text +  xlab(pourcentage2[1]) + ylab(pourcentage2[2])

dbplot




