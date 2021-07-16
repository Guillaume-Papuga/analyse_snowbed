#######################################################
# Project : Analysis of snowbed vegetation dynamics
# Script : 31.community_temporal_dynamic_wl
# Comparison of species gains and losses over different period of time (via bootstrap)
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

#######################################################
### Create the matrix
# transform the original dataset
varnames = colnames(df)

start_year = 2013
end_year = 2019

df_melt = melt (df, id=varnames[1:3]) %>% # first change to long format
  filter (year == start_year | year == end_year) %>% # and keep observation based on 2 dates, 2013 and 2019
  filter (!(year == start_year & value == 0)) %>% # delete the species not observed in 2013
  dcast (... ~ year) # create columns based the two selected years

colnames (df_melt) [4:5] = c("start", "end")

df_melt = df_melt %>%
  mutate (detect.na = is.na (start)) %>% # detect species with NA in the first year (absent)
  filter (detect.na != TRUE) %>% # delete those ligne
  replace_na(list(end = 0)) %>% # remplace missing species by 0 in the end
  select (-detect.na) %>% # delete the column
  mutate (freq_change = end - start) %>% # compute the change in frequency
  dplyr::rename (species = variable)

# synthesis table
df_synth = df_melt %>%
  group_by(species) %>%
  dplyr::summarize (avg.freq.chg = mean(freq_change)) %>% # compute the mean change in frequency
  add_column (boot_low_bnd = NA,
              boot_upp_bnd = NA, 
              boot_ave = NA) # add three empty columns for lower and higher boundaries + bootstrap average
df_synth = as.data.frame(df_synth)


# compute the number of observation per species
nb_obs = as.data.frame(table(df_melt[, c("species")])) %>% #compute the number of observed changes
  dplyr::rename(species = Var1,  # rename variables
                nb_obs = Freq) %>%
  arrange (desc(nb_obs)) %>% # sort variables
  filter(nb_obs > 5) # fix the threshold to keep or delete species (how many obs?)

### Run the bootstrap procedure
# Create all the functions
# create a function basic "average"
moyenne <- function(d,w) {
  n <- length(d)
  return(sum(d[w[1:n]])/n)}

# create a function to compute the average based on bootstrap
fun.boot.mean <- function(x) {
  res = boot(x, moyenne, R=10000)
  #return(res$t0)
  return(as.vector(res$t0))}

# compute the lower boundary of the confidence interval
fun.boot.lower <- function(x) {
  res = boot(x, moyenne, R=10000)
  #return(res$t0)
  return(as.vector(quantile(res$t,probs = 0.025)))}

# compute the higher boundary of the confidence interval
fun.boot.upper <- function(x) {
  res = boot(x, moyenne, R=10000)
  #return(res$t0)
  return(as.vector(quantile(res$t,probs = 0.975)))}

# select all the species for which you want to compute the bootstrap procedure
species = as.vector(nb_obs$species) 

for (i in species){
  # select the data
  data.boot = df_melt %>%
    filter (species == i)
  vec = data.boot$freq_change
  
  # compute the indices
  bt.ave = fun.boot.mean(vec)
  bt.low = fun.boot.lower(vec)
  bt.upp = fun.boot.upper(vec)
  
  # paste in the synth dataframe
  df_synth[which(df_synth$species == i), "boot_ave"] = bt.ave
  df_synth[which(df_synth$species == i), "boot_low_bnd"] = bt.low
  df_synth[which(df_synth$species == i), "boot_upp_bnd"] = bt.upp
}

df_synth = df_synth %>%
  drop_na() %>% #remove lines without data
  arrange(desc(boot_ave)) # %>% #sort in descendind order

# Create the column "winner stable loser - WSL"¨
wsl = function (x, y){
if(x > 0){
  print("W")
} else {
  if (y < 0){
    print ("L")
  } else {
    print ("S")
  }
}}

df_synth = df_synth %>%
  mutate (wsl = NA) # create the column WSL
for (i in 1:nrow(df_synth)){ # use a loop to fill the column with the correct status (winner, stable, loser)
  x = df_synth[i,"boot_low_bnd"]
  y = df_synth[i,"boot_upp_bnd"]
  df_synth[i, "wsl"] = wsl(x, y)
  }
  
# Plot
theme = theme(panel.background = element_blank(),
             panel.border=element_rect(fill=NA),
             panel.grid.major = element_blank(),
             panel.grid.minor = element_blank(),
             strip.background=element_blank(),
             axis.title.y = element_text(size = 16),
             axis.title.x = element_text(size = 16),
             axis.text.x=element_text(colour="black"),
             axis.text.y=element_text(size =10),
             axis.ticks=element_line(colour="black"),
             plot.margin=unit(c(1,1,1,1),"line"),
             legend.position = "none")

wl = ggplot()+ 
  theme + 
  geom_point(data=df_synth, aes(x = boot_ave, y = reorder (species, boot_ave), colour = factor(wsl), size = 3, stroke = 1))  +
  geom_errorbar(data=df_synth, 
              aes( y=reorder(species, boot_ave), 
                   xmin= boot_low_bnd, 
                   xmax = boot_upp_bnd, 
                   colour = factor(wsl)), 
              size = 1, width=0.5) +
  scale_color_manual (values = c("firebrick3", "grey", "chartreuse4")) +
  geom_vline(xintercept = 0, linetype = "dotted", col = "red") + 
  xlim(-1,1) + 
  ylab("Species") + xlab("Variation") 
  
wl




