#######################################################
# Project : Analysis of snowbed vegetation dynamics
# Script : 20.community_ecology
# Analysis of environmental drivers of plant community structure
# Authors : Guillaume Papuga & Thomas Masclaux
# Last update : 28 april 2021
#######################################################

### Load data 
df.flo = read.csv(here::here ("data", "processed", "sb_data_cast.csv"), 
                  head = T, sep = ",", dec = ".")
df.flo = df.flo %>%
  filter(year == 2019) %>% # keep one year 
  arrange (combe) %>% # sort the table by combe name
  dplyr::filter (!(combe == "cre" & placette == 1))

df.clim = read.csv(here::here ("data", "processed", "clim.data.csv"),  # load the file
                   head = T, sep = ",", dec = ".") %>%
  mutate(combe = tolower (combe)) %>% # change the name of the combe to lower case...
  mutate(combe = replace(combe, combe == "cas", "pdlc")) %>% # and correct 'Pas de la Case'
  filter(combe != "por") # delete Portillon from the climatic dataset (no data)
  
df.microclim = read.csv(here::here ("data", "processed", "synth.microclim.rda.csv"),  # load the file
                        head = T, sep = " ", dec = ".") %>%
  dplyr::select(site, plot, n_vege, n_f.day, ht.95)

######################### A. Select ecological variables ######################################
### 1. Multivariate analysis 
res.pca = dudi.pca(df.clim %>% dplyr::select (-combe), scannf = FALSE, nf = 2) 

# Compute weight of axes and species on axes
axes.wt = round (res.pca$eig / sum(res.pca$eig)*100, 1)
var.wt = round(res.pca$co, 1)

# Draw a plot
s.corcircle(res.pca$co)

### 2. Pairwise correlation
cor.data.clim = cor (df.clim %>% dplyr::select (-combe))
c.plot = corrplot(cor.data.clim, method = "circle", 
                  tl.col="black", tl.srt=45, tl.cex = 0.4)
# Save the plot
jpeg (here::here ("outputs", "figures", "appendix_cor.clim.jpg"), 
      quality = 100, width = 480, height = 480) # Open jpeg file
corrplot(cor.data.clim, method = "circle", 
         tl.col="black", tl.srt=45, tl.cex = 0.4)
dev.off() # 3. Close the file

pdf (here::here ("outputs", "figures", "appendix_cor.clim.pdf")) # Open pdf file
corrplot(cor.data.clim, method = "circle", 
         tl.col="black", tl.srt=45, tl.cex = 0.4)
dev.off() # 3. Close the file

### 3. First step selection
cor.data.clim.2 = cor (df.clim %>% dplyr::select (GDD_0_Pyrenees, 
                                                  PET_Annual_Pyrenees, PET_Summer_Pyrenees, 
                                                  Pot_solar_rad_Summer_Pyrenees, 
                                                  Precipitation_Annual_Pyrenees,
                                                  Temperature_max_Annual_Pyrenees, Temperature_max_Summer_Pyrenees,
                                                  Temperature_mean_Annual_Pyrenees,
                                                  Temperature_min_Annual_Pyrenees, Temperature_min_Summer_Pyrenees,
                                                  Water_availability_Annual_Pyrenees))
c.plot = corrplot(cor.data.clim.2, method = "circle", 
                  tl.col="black", tl.srt=45, tl.cex = 0.6)

### 4. Final selection
clim.select = df.clim %>% dplyr::select (combe, 
                                         PET_Annual_Pyrenees, 
                                         Pot_solar_rad_Summer_Pyrenees, 
                                         Precipitation_Annual_Pyrenees, 
                                         Temperature_max_Summer_Pyrenees,
                                         Temperature_min_Summer_Pyrenees)

c.plot = corrplot(cor(clim.select %>% dplyr::select (-combe)), 
                  method = "circle", 
                  tl.col="black", tl.srt=45, tl.cex = 0.6)

# PCA
res.pca = dudi.pca(clim.select%>% dplyr::select (-combe), 
                   scannf = FALSE, nf = 2) 

jpeg (here::here ("outputs", "figures", "appendix_corcircle.clim.jpg")) # Open jpeg file
s.corcircle(res.pca$co, grid = F) # Draw a plot
dev.off() # 3. Close the file

pdf(here::here("outputs", "figures", "appendix_corcircle.clim.pdf"))
s.corcircle(res.pca$co, grid = F) # Draw a plot
dev.off()

# Compute weight of axes and species on axes
axes.wt = round (res.pca$eig / sum(res.pca$eig)*100, 1)
var.wt = round(res.pca$co, 1)

# Extend to the correct number of row (one data per quadrat)
clim.full = clim.select %>%
  full_join(df.flo %>% dplyr::select (combe, placette), by = "combe") %>%
  arrange (combe) %>% # sort data by combe name
  dplyr:: rename (pet = PET_Annual_Pyrenees,
                  rad_sum = Pot_solar_rad_Summer_Pyrenees, 
                  prec = Precipitation_Annual_Pyrenees, 
                  t_max = Temperature_max_Summer_Pyrenees, 
                  t_min = Temperature_min_Summer_Pyrenees)

######### REVIEW UPDATE ###########
######### Microclimate ############

clim.full = clim.full %>%
  inner_join(df.microclim, 
            by = c("combe" = "site", 
                   "placette" = "plot")) %>%
  dplyr:: rename (veg_per = n_vege, # Rename the variables to ease the reading
                  f_day = n_f.day, 
                  loc.t_max = ht.95)

######################### B. Run redundancy analysis ######################################
### 1. Multivariate analysis 
# Fit a first model
rda.clim = rda(df.flo %>% dplyr::select (-combe, -placette, -year) ~ 
                 pet + rad_sum + prec + t_max + t_min + veg_per + f_day + loc.t_max, 
               clim.full, dist="bray")
vif.cca(rda.clim) # Variance inflation factor

# Exclude pet and fit a second model
rda.clim = rda(df.flo %>% dplyr::select (-combe, -placette, -year) ~ 
                 rad_sum + prec + t_max + t_min + veg_per + f_day + loc.t_max, 
               clim.full, dist="bray")
vif.cca(rda.clim) # Variance inflation factor

# Test the significance of each axis
anova.cca(rda.clim, step=5000, by="axis", permutations = 500)

# Test the significance of each variable
test.rda = anova.cca(rda.clim, by="terms", permutations = 500)
test.rda
write.table(as.data.frame(test.rda), here::here ("outputs", "figures", "table_test.rda.csv"))

# Summary of the analysis
smry = summary(rda.clim)

eig_rdacomp = round(smry$cont$importance[,1:3]*100,2)
pourcentage2 <- paste(colnames(eig_rdacomp), "(", paste( as.character(eig_rdacomp[2,]), "%", ")", sep=""))

coef(rda.clim) # coefficient de reajustement

rda.clim$colsum

scrs = scores(rda.clim) # site & species scores

### 2. Plot the result of the RDS
# a. Matrix of plots' coordinates
df1 = data.frame(smry$sites[,1:2]) %>% #
  mutate (site = df.flo$combe, 
          plot = as.factor(df.flo$placette))

# b. Matrix of climatic variables' coordinate
df2 = data.frame(smry$biplot[,1:2])

# c. Matrix of species coordinates
df3 = data.frame(smry$species[,1:2]) %>%
  rownames_to_column('species') %>%
  mutate (arrow_length = sqrt((RDA1)^2+(RDA2^2))) %>%
  arrange(desc(arrow_length)) %>%
  filter(arrow_length > 0.3) %>%
  column_to_rownames('species')

# d. Create de plot
# set up the theme
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

# Draw arrows for variables
gg = geom_segment(data = df2, aes(x=0, xend = RDA1, y = 0, yend = RDA2), 
                  color ="grey50", arrow=arrow(length = unit(0.01, "npc")))

text = geom_text(data=df2, aes(x=RDA1, y=RDA2, label = rownames(df2),
                               hjust=0.5*(1-sign(RDA1)),
                               vjust=0.5*(1-sign(RDA2))),
                 color="black", size = 5)

# Draw arrows for species
sp = geom_segment(data = df3, aes(x = 0, xend = RDA1, y = 0, yend = RDA2), 
                  arrow = arrow(length = unit(0.01,"npc")),linetype=1, 
                  size=0.1,colour = "red") 

text_sp =  geom_text(data = df3, aes(x=RDA1,y=RDA2,label=row.names(df3), 
                                     hjust=0.5*(1-sign(RDA1)),
                                     vjust=0.5*(1-sign(RDA2))), 
                     colour = "grey")

# Draw the main plot (based on quadrats) and add species and variables
rda.plot = ggplot(df1, aes(x=RDA1, y=RDA2)) +
  geom_point(aes(colour = site, shape = plot), size=3) +  
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_vline(xintercept=0, linetype="dotted") + 
  theme+ text +  xlab(pourcentage2[1]) + ylab(pourcentage2[2]) + 
  gg  + 
  sp + 
  text_sp 

rda.plot

# e. Save the plot
jpeg (here::here ("outputs", "figures", "figure_rda.plot_update.jpg")) # Open jpeg file
print(rda.plot)
dev.off() # 3. Close the file

pdf(here::here("outputs", "figures", "figure_rda.plot_update.pdf"))
print(rda.plot)
dev.off()


