#######################################################
# Project : Analysis of snowbed vegetation dynamics
# Script : 32.community_temporal_dynamic_temporal.changes
# Analysis of changes through time of plant communities, change of TBI through time + decomposition
# Authors : Guillaume Papuga & Thomas Masclaux
# Last update : 28 april 2021
#######################################################

### Load data 
df = read.csv(here::here ("data", "processed", "sb_data_cast.csv"), 
              head = T, sep = ",", dec = ".")

df = df %>% # first we delete the placette n°1 from cre
  mutate(name.pla = interaction(combe, placette)) %>%
  filter (name.pla != "cre.1") %>%
  select (-name.pla)

########################################################

########################################################
##################I. Loss and gain #####################
########################################################
### Create the BCD table
# Chose the different sites
nb_obs = as.data.frame.matrix(table (df[,c("combe", "year")])) # create the table
nb_obs
site = row.names (nb_obs) [-c(1, 6)]  # store the different sites

# Create a matrix to store results
bcd = as.data.frame(matrix(nrow = 0, ncol = 5, data = NA)) %>%
  setNames (c("combe", "period", "b", "c", "d"))

# Run a loop
for (i in site) {
  # create the base vegetation releve
  start_yr = df %>%
    filter (combe == i, year == 2013) %>% # select the correct portion of the dataset
    mutate (name = interaction(combe, placette)) # %>% # create the name of each quad

  # create a vector for years with observations
  survey_yr = unique (df %>% 
                        filter (combe == i) %>% 
                        select (year) %>% 
                        filter (year != 2013))$year
  
  # Run a second loop for each pair of survey
  for (j in survey_yr) {
    # create the test vegetation releve (comparison of each given year)
    survey = df %>%
      filter (combe == i, year == j) %>% # select the correct portion of the dataset
      mutate (name = interaction(combe, placette)) #%>% # create the name of each quad

    # fit the two table
    common_quad = intersect(start_yr$name, survey$name) # list of common quadrat
    m_start = subset(start_yr, name %in% common_quad) %>%
      select (-combe, -placette, -year, -name) # delete useless columns
    
    m_year = subset(survey, name %in% common_quad) %>%
      select (-combe, -placette, -year, -name) # delete useless columns

    # compute the tbi index
    tbindex = TBI(m_start, m_year, 
                  method = "%difference", pa.tr = F, nperm = 999,
                  BCD = TRUE, save.BC = T,test.BC = TRUE, test.t.perm = T, seed. =   )

    # store
    indices = as.data.frame(matrix(nrow = 1, ncol = 5, 
                                   data = c(i, # combe name
                                            paste ("13-", substr(j, 3, 4), sep = ""),  # transition
                                            as.numeric(tbindex$BCD.summary$`mean(B/den)`), # B componant 
                                            tbindex$BCD.summary$`mean(C/den)`, # C componant
                                            tbindex$BCD.summary$`mean(D)`))) %>% # D componant
      setNames (c("combe", "period", "b", "c", "d"))
    bcd = rbind (bcd, indices)
      }
}

# Convert the dataframe into numeric values
bcd$period = as.character(bcd$period)
bcd$b = as.numeric(as.vector(bcd$b))
bcd$c = as.numeric(as.vector(bcd$c))
bcd$d = as.numeric(as.vector(bcd$d))

# Create year that haven't been surveyed for the plot
acast (bcd, period  ~ combe) # show the missing values
tot_surv = expand.grid(unique (bcd$combe), unique (bcd$period))  # dataframe of all combinations

bcd = bcd %>%
  right_join(tot_surv, by = c("combe" = "Var1", 
                              "period" = "Var2")) %>% # by adding the full df it creates NAs when no sample was done
  arrange (desc(combe)) # sort the df


########################################################
################# II. Intraplot dissimilarity ##########
########################################################
### Create the BCD table
# Chose the different sites
nb_obs = as.data.frame.matrix(table (df[,c("combe", "year")])) # create the table
nb_obs
site = row.names (nb_obs) [-c(1, 6)]  # store the different sites

# Create a matrix to store results
inter.quad.dist = as.data.frame(matrix(nrow = 0, ncol = 4, data = NA)) %>%
  setNames (c("combe", "year", "comp", "distance"))

# Run a loop
for (i in site) {
  # create the base vegetation releve
  df.comp.site = df %>%
    filter (combe == i)  # select the correct portion of the dataset

  # create a vector for years with observations
  survey_yr = unique (df.comp.site$year)
  
  # second loop to compute each index per year
  for (j in survey_yr) {
    # Create the matrix
    mat = df.comp.site %>%
      filter (year == j)
    
    # Comparison 1 - 2
    if (nrow (mat %>% filter(placette == 1)) !=0 && # filter to avoid computing missing values
        nrow (mat %>% filter(placette == 2)) !=0) {
      tab = mat %>%
        filter (placette == 1  | placette == 2) %>% # select quadrat 1 and 2
        select (-combe, -placette, -year)
      dist.quad = ecodist::distance(tab, "bray-curtis") # create a distance matrix
      
      # store the value
      index = as.data.frame(matrix(nrow = 1, ncol = 4, # create a dataframe
                                   data = c(i, # combe name
                                            j,  # year
                                            "1-2", # comparison
                                            as.numeric(dist.quad)))) %>% # D componant
        setNames (c("combe", "year", "comp", "distance"))
      
      inter.quad.dist = rbind (inter.quad.dist, index) # and bind it to the main storage DF
    }
    
    # Comparison 2 - 3
    if (nrow (mat %>% filter(placette == 2)) !=0 && # filter to avoid computing missing values
        nrow (mat %>% filter(placette == 3)) !=0) {
      tab = mat %>%
        filter (placette == 2  | placette == 3) %>% # select quadrat 1 and 2
        select (-combe, -placette, -year)
      dist.quad = ecodist::distance(tab, "bray-curtis") # create a distance matrix
      
      # store the value
      index = as.data.frame(matrix(nrow = 1, ncol = 4, # create a dataframe
                                   data = c(i, # combe name
                                            j,  # year
                                            "2-3", # comparison
                                            as.numeric(dist.quad)))) %>% # D componant
        setNames (c("combe", "year", "comp", "distance"))
      
      inter.quad.dist = rbind (inter.quad.dist, index) # and bind it to the main storage DF
    }
    
    # Comparison 1 - 3
    if (nrow (mat %>% filter(placette == 1)) !=0 && # filter to avoid computing missing values
        nrow (mat %>% filter(placette == 3)) !=0) {
      tab = mat %>%
        filter (placette == 1  | placette == 3) %>% # select quadrat 1 and 2
        select (-combe, -placette, -year)
      dist.quad = ecodist::distance(tab, "bray-curtis") # create a distance matrix
      
      # store the value
      index = as.data.frame(matrix(nrow = 1, ncol = 4, # create a dataframe
                                   data = c(i, # combe name
                                            j,  # year
                                            "1-3", # comparison
                                            as.numeric(dist.quad)))) %>% # D componant
        setNames (c("combe", "year", "comp", "distance"))
      
      inter.quad.dist = rbind (inter.quad.dist, index) # and bind it to the main storage DF
    }
    
  }
}

# Convert the dataframe into numeric values
inter.quad.dist$combe = as.character(inter.quad.dist$combe)
inter.quad.dist$year = as.numeric(as.vector(inter.quad.dist$year))
inter.quad.dist$comp = as.character(as.vector(inter.quad.dist$comp))
inter.quad.dist$distance = as.numeric(as.vector(inter.quad.dist$distance))


# Organize the main table
# Create year that haven't been surveyed for the plot
acast (inter.quad.dist, combe  ~ year) # show the missing values
tot_surv = expand.grid(unique (inter.quad.dist$combe), 
                       unique (inter.quad.dist$year))  # dataframe of all combinations




########################################################
##################III. Plot ############################
########################################################

#######################################################
### PLOT 1
### Define the theme of plots
# theme<-theme(panel.background = element_blank(),
#              panel.border=element_rect(fill=NA),
#              panel.grid.major = element_blank(),
#              panel.grid.minor = element_blank(),
#              strip.background=element_blank(),
#              axis.text.x=element_text(colour="black"),
#              axis.text.y=element_text(colour="black"),
#              axis.ticks=element_line(colour="black"),
#              plot.margin=unit(c(1,1,1,1),"line"))


theme<-theme(panel.background = element_blank(),
             panel.border=element_rect(fill=NA),
             panel.grid.major = element_blank(),
             panel.grid.minor = element_blank(),
             strip.background=element_blank(),
             axis.text.x=element_text(colour="black"),
             axis.text.y=element_text(colour="black"),
             axis.ticks=element_line(colour="black"),
             axis.title.x = element_blank(),
             axis.title.y = element_blank(),
             plot.margin=unit(c(1,1,1,1),"line"))


# axis.title.x=element_blank(),
# axis.text.x=element_blank(),
# axis.ticks.x=element_blank()

### Create the first serie of plots : BCD 
site = as.vector(unique (bcd$combe))

for (i in site) {
  q1 = ggplot() + #premiere couche
    geom_point(data = bcd[which(bcd$combe == i),], aes(x = period, y = b), colour = "red", size=3, shape = 16) + 
    scale_y_continuous(limits= c(0,.6), breaks = seq(0,0.6,0.2)) +
    geom_line(data = bcd[which(bcd$combe == i),], aes(x = period, y = b, group=1), colour = "red")+
    theme 
  
  q2 = ggplot() + # seconde couche
    geom_point(data = bcd[which(bcd$combe == i),], aes(x = period, y = c), colour = "skyblue", size=3, shape=17) + 
    scale_y_continuous(limits= c(0,.6), breaks = seq(0,0.6,0.2)) +
    geom_line(data = bcd[which(bcd$combe == i),], aes(x = period, y = c, group=2), colour = "skyblue") 
  
  q3 = ggplot() + #troisieme couche
    geom_point(data = bcd[which(bcd$combe == i),], aes(x = period, y = d), colour = "black", size=3, shape=15) + 
    scale_y_continuous(limits= c(0,.6), breaks = seq(0,0.6,0.2)) +
    geom_line(data = bcd[which(bcd$combe == i),], aes(x = period, y = d, group=3), colour = "black") 
  
  bcd.plot = q1 + q2$layers[[1]]+ q2$layers[[2]] + q3$layers[[1]] + q3$layers[[2]] # combine layers
  assign( paste("bcd.plot", i, sep = "_"), # rename the plot
          bcd.plot)
  
}

#######################################################
### PLOT 2
### Create the second serie of plots : interplot dissimilarities 
site = as.vector(unique (inter.quad.dist$combe))

for (i in site) {
  mat.plot = inter.quad.dist %>%
    filter (combe == i)
  
  # Create each plot independantly
  # Plot 1-2
  if (nrow (mat.plot %>% filter(comp == "1-2")) !=0) {
    inter.plot_12 = ggplot() + 
      geom_point(data = mat.plot[which (mat.plot$comp == "1-2"),], aes(x = year, y = distance), colour = "red", size=3, shape = 16) +  
      scale_y_continuous(limits= c(0,1), breaks = seq(0,1,.5)) +
      geom_line(data = mat.plot[which (mat.plot$comp == "1-2"),], aes(x = year, y = distance, group=1), colour = "red")+
      theme(axis.title.x = element_blank())+
      theme(axis.title.y = element_blank(), axis.text.y=element_text(size =16), axis.text.x=element_text(size =16)) +
      theme 
  } 
  
  # Plot 2-3
  if (nrow (mat.plot %>% filter(comp == "2-3")) !=0) {
    inter.plot_23 = ggplot() + 
      geom_point(data = mat.plot[which (mat.plot$comp == "2-3"),], aes(x = year, y = distance), colour = "blue", size=3, shape = 16) +  
      scale_y_continuous(limits= c(0,1), breaks = seq(0,1,.5)) +
      geom_line(data = mat.plot[which (mat.plot$comp == "2-3"),], aes(x = year, y = distance, group=1), colour = "blue")+
      theme(axis.title.x = element_blank())+
      theme(axis.title.y = element_blank(), axis.text.y=element_text(size =16), axis.text.x=element_text(size =16)) +
      theme 
  } 
  
  # Plot 1-3
  if (nrow (mat.plot %>% filter(comp == "1-3")) !=0) {
    inter.plot_13 = ggplot() + 
      geom_point(data = mat.plot[which (mat.plot$comp == "1-3"),], aes(x = year, y = distance), colour = "green", size=3, shape = 16) +  
      scale_y_continuous(limits= c(0,1), breaks = seq(0,1,.5)) +
      geom_line(data = mat.plot[which (mat.plot$comp == "1-3"),], aes(x = year, y = distance, group=1), colour = "green")+
      theme(axis.title.x = element_blank())+
      theme(axis.title.y = element_blank(), axis.text.y=element_text(size =16), axis.text.x=element_text(size =16)) +
      theme 
  } 
  
  # Draw one plot per combe
  site.plot = ggplot() +
    scale_y_continuous(limits= c(0,1), breaks = seq(0,1,.5)) +
    theme
  
  # as some plot do not exist (quadrat missing) only add layer IF the layer exists... and build the plot sequentially
  if (exists("inter.plot_12") == TRUE) {
    site.plot = site.plot +
      inter.plot_12$layers[[1]] + inter.plot_12$layers[[2]]}
  
  if (exists("inter.plot_23") == TRUE) {
    site.plot = site.plot + 
      inter.plot_23$layers[[1]] + inter.plot_23$layers[[2]]}
  
  if (exists("inter.plot_13") == TRUE) {
    site.plot = site.plot + 
      inter.plot_13$layers[[1]] + inter.plot_13$layers[[2]] }
  
  # Change the name of the plot
  assign(paste ("site.plot", i, sep = "_"), 
         site.plot)
  
  # Erase all plots to avoid re-using curves
  if (exists("inter.plot_12") == TRUE) {rm(inter.plot_12)}
  if (exists("inter.plot_23") == TRUE) {rm(inter.plot_23)}
  if (exists("inter.plot_13") == TRUE) {rm(inter.plot_13)}
}


### Assemble the global plot
temp.dissi = plot_grid(site.plot_cat, bcd.plot_cat,
                       site.plot_cre, bcd.plot_cre, 
                       site.plot_pdlc, bcd.plot_pdlc, 
                       site.plot_pla, bcd.plot_pla,
                       site.plot_rat, bcd.plot_rat, 
                       site.plot_ull, bcd.plot_ull, 
                       nrow = 6, ncol = 2) + 
  theme(plot.margin = unit(c(1,0.5,0.5,1), "cm")) + # add some extra space to write everything
  draw_label("Intraplot dissimilarity", fontface = 'bold', x = 0.13, y = 1, hjust = 0) + 
  draw_label("Loss, gain, dissimilarity", fontface = 'bold', x = 0.63, y = 1, hjust = 0) + 
  draw_label("CAT", fontface = 'bold', x = 0, y = 0.925, angle = 90, vjust = 0, size = 12) +
  draw_label("CRE", fontface = 'bold', x = 0, y = 0.75, angle = 90, vjust = 0, size = 12) +
  draw_label("PDLC", fontface = 'bold', x = 0, y = 0.59, angle = 90, vjust = 0, size = 12) +
  draw_label("PLA", fontface = 'bold', x = 0, y = 0.43, angle = 90, vjust = 0, size = 12) +
  draw_label("RAT", fontface = 'bold', x = 0, y = 0.26, angle = 90, vjust = 0, size = 12) +
  draw_label("ULL", fontface = 'bold', x = 0, y = 0.1, angle = 90, vjust = 0, size = 12)


# xlab("year of surveys") + ylab("Ruzicka dissimilarity") +  # first serie of plots