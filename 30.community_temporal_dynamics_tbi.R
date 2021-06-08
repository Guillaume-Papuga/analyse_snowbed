#######################################################
# Project : Analysis of snowbed vegetation dynamics
# Script : 30.community_temporal_dynamic_tbi
# Analysis of changes through time of plant communities based on Legendre TBI
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


######################### Temporal Beta Diversity Index #############################################
### 


theme<-theme(panel.background = element_blank(),
             panel.border=element_rect(fill=NA),
             panel.grid.major = element_blank(),
             panel.grid.minor = element_blank(),
             strip.background=element_blank(),
             axis.text.x=element_text(colour="black"),
             axis.text.y=element_text(colour="black"),
             axis.ticks=element_line(colour="black"),
             plot.margin=unit(c(1,1,1,1),"line"))


## splitting df into two moments of survey table
bef2013<-df %>% filter(grepl('2013', year)) #2013 survey
bef2013$year <- NULL
fac <- select_if(bef2013,is.factor)
name<- paste(fac$combe,fac$placette, sep="_")
bef2013 <- cbind.data.frame(name, bef2013)
rownames(bef2013) = make.names(name, unique=TRUE)

bef2019<-df %>% filter(grepl('2019', year)) #2019 survey
bef2019$year <- NULL
fac <- select_if(bef2019,is.factor)
name<- paste(fac$combe,fac$placette, sep="_")
bef2019 <- cbind.data.frame(name, bef2019)
rownames(bef2019) = make.names(name, unique=TRUE)

#TBI 2013 to 2019
P12013<-bef2013 
P12019<-bef2019 
P12019 <- subset(P12019, combe %in% intersect(P12013$combe, P12019$combe))
P12013 <- subset(P12013, combe %in% intersect(P12019$combe, P12013$combe))
P12019 = P12019[-2,]
m2019 <- as.matrix(P12019[,4:99])
m2013 <- as.matrix(P12013[,4:99])
tbi_201319P1 <- TBI(m2013,m2019, method = "%difference",pa.tr = F, nperm = 999,BCD = TRUE, save.BC = T,test.BC = TRUE, test.t.perm = T, seed. =   )
tbi_201319P1
tbi_201319P1$p.TBI

#tbi signifance table
test_sign = rbind.data.frame(tbi_201319P1$TBI, tbi_201319P1$p.TBI, tbi_201319P1$p.adj)
colnames(test_sign) = rownames(P12019)
test_sign = t(test_sign)

#construction du B-C Plot et ?tiquettes
row.names(tbi_201319P1$BCD.mat) = row.names(P12019) #label plot on B-C matrix
col1= row.names(P12019)
col1 <- gsub("_","-",rownames(P12019))
col1 = as.data.frame(col1)
col1 <-separate(col1,col=1,into = c("combe","placette"), sep = "-") #combe and plot group labels

tbi_201319P1$BCD.mat = cbind.data.frame(tbi_201319P1$BCD.mat, col1)
tbi_201319P1$BCD.mat$combe = as.factor(tbi_201319P1$BCD.mat$combe)
tbi_201319P1$BCD.mat

plot(tbi_201319P1, s.names=rownames(tbi_201319P1$BCD.mat)  ,  #plot B-C
     col.rim = "coral", 
     pch.loss=19 , 
     pch.gain=15,
     main="Temporal change in beta diversity", xlim=c(0,0.5), ylim=c(0,0.5))


