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
df_melt = melt (df, id=varnames[1:3]) %>% # first change to long format
  filter (year == 2013 | year == 2019) %>% # and keep observation based on 2 dates, 2013 and 2019
  filter (!(year == 2013 & value == 0)) %>% # delete the species not observed in 2013
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

  
  
  
### Plot
theme<-theme(panel.background = element_blank(),
             panel.border=element_rect(fill=NA),
             panel.grid.major = element_blank(),
             panel.grid.minor = element_blank(),
             strip.background=element_blank(),
             axis.title.y = element_text(size = 16),
             axis.title.x = element_text(size = 16),
             axis.text.x=element_text(colour="black"),
             axis.text.y=element_text(size =16),
             axis.ticks=element_line(colour="black"),
             plot.margin=unit(c(1,1,1,1),"line"))

  
##############################""""
selec_esp = as.data.frame.matrix(table(df_melt$variable, df_melt$combe)) %>%
  mutate (tot_plac = rowSums(...))


##W-L species variation 2013-2019

theme<-theme(panel.background = element_blank(),
             panel.border=element_rect(fill=NA),
             panel.grid.major = element_blank(),
             panel.grid.minor = element_blank(),
             strip.background=element_blank(),
             axis.title.y = element_text(size = 16),
             axis.title.x = element_text(size = 16),
             axis.text.x=element_text(colour="black"),
             axis.text.y=element_text(size =16),
             axis.ticks=element_line(colour="black"),
             plot.margin=unit(c(1,1,1,1),"line"))

wl = ggplot()+ 
  geom_point(data=df_synth, aes(x = boot_ave, y = reorder (species, boot_ave), size = 4))  +
  geom_errorbar(data=df_synth, 
              aes( y=reorder(species, boot_ave), 
                   xmin= boot_low_bnd, 
                   xmax = boot_upp_bnd), 
              size = 0.5, width=0.5) # +
  
scale_shape_manual(values = c(19, 1))+
  
  scale_color_manual(name="legends", 
                     labels = c("deprived", "winner"), 
                     values = c("deprived"="firebrick3", "winner"="chartreuse4")) +
  
  scale_fill_manual(values="white")+
  
  theme + 
  ylim(-1,1) + 
  coord_flip()+
  xlab("Species") + ylab("variation") +
  scale_x_discrete(position = "top") +
  geom_hline(yintercept = 0, linetype = "dotted")

WL_PLOT #plot abundnaces dissimilarities   





#####################################################################""


# melt 2013 and 2019
m.m1 = melt(m2013)
colnames (m.m1) [3] = "v1"
m.m2 = melt(m2019)
colnames (m.m2) [3] = "v2"

# join
j.m = full_join (m.m1, m.m2) %>%
  mutate(verif = v1 == v2) %>% #on v?rifie si la variation est r?elle sur une placette
  mutate(val = v2 - v1)

t = j.m %>% filter(verif != TRUE | v1 !=0)
t

t_selct = t %>% group_by(Var2) %>% filter(n() >= 4)
t_selct = as.data.frame(t_selct)

# test de normalit?
test = do.call("rbind", with(t_selct, 
                             tapply(val, Var2,
                                    function(x) unlist(shapiro.test(x)[c("statistic", "p.value")])))) #test de normalit?
test = as.data.frame(test) # conversion du tableau avant assignation
test$significance<- ifelse(test$p.value >0.05, "*","no_sign") # on assigne les distributions significatives
test$species = row.names(test) # on renomme les lignes


x = t %>%  dplyr::group_by(Var2) %>% dplyr::summarize(m.test = mean(val)) #moyenne des valeurs par esp?ces

ic = t %>% dplyr::group_by(Var2) %>% dplyr::summarize(m.test = 1.96*std.error(val)) # erreur standard par esp?ces

colnames(ic)[1] = "species" # renommer la colonne
colnames(ic)[2]= "sd"       # renommer la colonne
ic$species = as.character(ic$species)


x$WL<- ifelse(x$m.test >0, "winner","deprived") # assignation des esp?ces perdantes et gagnantes

x = x %>% arrange(desc(m.test)) #on range par odre d?croissant pour pr?parer le graphique

##seuil de repr?sentativit? des esp?ces sur les placettes
barplot(sort(colSums(table(t[,1:2])))) 

tt =  data.frame(sort(colSums(table(t[,1:2])))) #compter N placettes o? une esp?ce est vue
tt$species = row.names(tt) #nom esp?ce en t?te de ligne

species = which(tt$sort.colSums.table.t...1.2....>7) #seuil ? 7 ?chantillons
tt = tt[species, ]
tt$species = as.character(tt$species)
colnames(x)[1] <- "species"
x$species = as.character(x$species)


xt = semi_join(x,tt, by = "species") #fusion en base de construction du plot

xt = left_join (xt, test, "species")
xt=left_join(xt, ic, "species")


##plot W-L

WL_PLOT = ggplot()+ 
  
  geom_point(data=xt, aes( x=reorder(species, m.test),
                           y = m.test, 
                           shape = significance, 
                           color = WL), size = 4) +
  scale_shape_manual(values = c(19, 1))+
  
  scale_color_manual(name="legends", 
                     labels = c("deprived", "winner"), 
                     values = c("deprived"="firebrick3", "winner"="chartreuse4")) +
  
  scale_fill_manual(values="white")+
  geom_errorbar(data=xt, 
                aes( x=reorder(species, m.test), 
                     ymin= m.test-sd, 
                     ymax = m.test+sd, color = WL), 
                size = 0.5, width=0.5) +
  
  theme + 
  ylim(-1,1) + 
  coord_flip()+
  xlab("Species") + ylab("variation") +
  scale_x_discrete(position = "top") +
  geom_hline(yintercept = 0, linetype = "dotted")

WL_PLOT #plot abundnaces dissimilarities   

