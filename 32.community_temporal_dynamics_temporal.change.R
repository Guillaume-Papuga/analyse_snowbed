#######################################################
# Project : Analysis of snowbed vegetation dynamics
# Script : 32.community_temporal_dynamic_temporal.changes
# Analysis of changes through time of plant communities, change of TBI through time + decomposition
# Authors : Guillaume Papuga & Thomas Masclaux
# Last update : 28 april 2021
#######################################################

#5) Analyses temporelles----


theme<-theme(panel.background = element_blank(),
             panel.border=element_rect(fill=NA),
             panel.grid.major = element_blank(),
             panel.grid.minor = element_blank(),
             strip.background=element_blank(),
             axis.text.x=element_text(colour="black"),
             axis.text.y=element_text(colour="black"),
             axis.ticks=element_line(colour="black"),
             plot.margin=unit(c(1,1,1,1),"line"))


#############################################################################
#time dissimilarities
data = df %>% filter(grepl("ULL", site)) #write site name and run script , data for trajectories analysis


theme<-theme(panel.background = element_blank(),
             panel.border=element_rect(fill=NA),
             panel.grid.major = element_blank(),
             panel.grid.minor = element_blank(),
             strip.background=element_blank(),
             axis.text.x=element_text(colour="black"),
             axis.text.y=element_text(colour="black"),
             axis.ticks=element_line(colour="black"),
             plot.margin=unit(c(1,1,1,1),"line"))

# splitting DF into seven years of survey table####

bef2013<-as.data.frame(data %>% filter(grepl('2013', year)))
bef2013$year <- NULL
fac <- select_if(bef2013,is.factor)
name<-paste(fac$site,fac$placette, sep="_")
bef2013 <- cbind(name, bef2013)
rownames(bef2013) = make.names(name, unique=TRUE)


bef2014<-as.data.frame(data %>% filter(grepl('2014', year)))
bef2014$year <- NULL
fac <- select_if(bef2014,is.factor)
name<- paste(fac$site,fac$placette, sep="_")
bef2014 <- cbind(name, bef2014)
rownames(bef2014) = make.names(name, unique=TRUE)

bef2015<-as.data.frame(data %>% filter(grepl('2015', year)))
bef2015$year <- NULL
fac <- select_if(bef2015,is.factor)
name<- paste(fac$site,fac$placette, sep="_")
bef2015 <- cbind(name, bef2015)
rownames(bef2015) = make.names(name, unique=TRUE)

bef2016<-as.data.frame(data %>% filter(grepl('2016', year)))
bef2016$year <- NULL
fac <- select_if(bef2016,is.factor)
name<- paste(fac$site,fac$placette, sep="_")
bef2016 <- cbind(name, bef2016)
rownames(bef2016) = make.names(name, unique=TRUE)

bef2017<-as.data.frame(data %>% filter(grepl('2017', year)))
bef2017$year <- NULL
fac <- select_if(bef2017,is.factor)
name<- paste(fac$site,fac$placette, sep="_")
bef2017 <- cbind(name, bef2017)
rownames(bef2017) = make.names(name, unique=TRUE)

bef2018<-as.data.frame(data %>% filter(grepl('2018', year)))
bef2018$year <- NULL
fac <- select_if(bef2018,is.factor)
name<- paste(fac$site,fac$placette, sep="_")
bef2018 <- cbind(name, bef2018)
rownames(bef2018) = make.names(name, unique=TRUE)

bef2019<-as.data.frame(data %>% filter(grepl('2019', year)))
bef2019$year <- NULL
fac <- select_if(bef2019,is.factor)
name<- paste(fac$site,fac$placette, sep="_")
bef2019 <- cbind(name, bef2019)
rownames(bef2019) = make.names(name, unique=TRUE)


# see look at site getting 'one to seven' survey(s) because input matrix must have the same dimension####

##TBI between 2013&2014
m2014 <- subset(bef2014, name %in% intersect(bef2013$name, bef2014$name)) # DFP2.2 avec des lignes communes ? DFP1
m20131 <- subset(bef2013, name %in% intersect(m2014$name, bef2013$name)) # DFP1.2 avec des lignes communes ? DFP2.2
m2014 <- as.matrix(m2014[,4:99])
m20131 <- as.matrix(m20131[,4:99])
tbi_13_14 <- TBI(m20131,m2014, method = "ruzicka",pa.tr = F, nperm = 999,BCD = TRUE, save.BC = T,test.BC = TRUE, test.t.perm = T, seed. =   )
tbi_13_14 #loss and gain generated, tbi generated for 999 permutation of columns, ruzicka for incidence data

##TBI between 2013&2015
m2015 <- subset(bef2015, name %in% intersect(bef2013$name, bef2015$name)) # DFP2.2 avec des lignes communes ? DFP1
m20132 <- subset(bef2013, name %in% intersect(m2015$name, bef2013$name)) # DFP1.2 avec des lignes communes ? DFP2.2
m2015 <- as.matrix(m2015[,4:99])
m20132 <- as.matrix(m20132[,4:99])
tbi_13_15 <- TBI(m20132,m2015, method = "ruzicka",pa.tr = F, nperm = 999,BCD = TRUE, save.BC = T,test.BC = TRUE, test.t.perm = T, seed. =   )
tbi_13_15 #loss and gain generated, tbi generated for 999 permutation of columns, ruzicka for incidence data

##TBI between 2013&2016
m2016 <- subset(bef2016, name %in% intersect(bef2013$name, bef2016$name)) # DFP2.2 avec des lignes communes ? DFP1
m20133 <- subset(bef2013, name %in% intersect(m2016$name, bef2013$name)) # DFP1.2 avec des lignes communes ? DFP2.2
m2016 <- as.matrix(m2016[,4:99])
m20133 <- as.matrix(m20133[,4:99])
tbi_13_16 <- TBI(m20133,m2016, method = "ruzicka",pa.tr = F, nperm = 999,BCD = TRUE, save.BC = T,test.BC = TRUE, test.t.perm = T, seed. =   )
tbi_13_16 #loss and gain generated, tbi generated for 999 permutation of columns, ruzicka for incidence data

##TBI between 2013&2017
m2017 <- subset(bef2017, name %in% intersect(bef2013$name, bef2017$name)) # DFP2.2 avec des lignes communes ? DFP1
m20134 <- subset(bef2013, name %in% intersect(m2017$name, bef2013$name)) # DFP1.2 avec des lignes communes ? DFP2.2
m2017 <- as.matrix(m2017[,4:99])
m20134 <- as.matrix(m20134[,4:99])
tbi_13_17 <- TBI(m20134,m2017, method = "ruzicka",pa.tr = F, nperm = 999,BCD = TRUE, save.BC = T,test.BC = TRUE, test.t.perm = T, seed. =   )
tbi_13_17 #loss and gain generated, tbi generated for 999 permutation of columns, ruzicka for incidence data

##TBI between 2013&2018
m2018 <- subset(bef2018, name %in% intersect(bef2013$name, bef2018$name)) # DFP2.2 avec des lignes communes ? DFP1
m20135 <- subset(bef2013, name %in% intersect(m2018$name, bef2013$name)) # DFP1.2 avec des lignes communes ? DFP2.2
m2018 <- as.matrix(m2018[,4:99])
m20135 <- as.matrix(m20135[,4:99])
tbi_13_18 <- TBI(m20135,m2018, method = "ruzicka",pa.tr = F, nperm = 999,BCD = TRUE, save.BC = T,test.BC = TRUE, test.t.perm = T, seed. =   )
tbi_13_18 #loss and gain generated, tbi generated for 999 permutation of columns, ruzicka for incidence data

##TBI between 2013&2019
m2019 <- subset(bef2019, name %in% intersect(bef2013$name, bef2019$name)) # DFP2.2 avec des lignes communes ? DFP1
m20136 <- subset(bef2016, name %in% intersect(m2019$name, bef2013$name)) # DFP1.2 avec des lignes communes ? DFP2.2
m2019 <- as.matrix(m2019[,4:99])
m20136 <- as.matrix(m20136[,4:99])
tbi_13_19 <- TBI(m20136,m2019, method = "ruzicka",pa.tr = F, nperm = 999,BCD = TRUE, save.BC = T,test.BC = TRUE, test.t.perm = T, seed. =   )
tbi_13_19 #loss and gain generated, tbi generated for 999 permutation of columns, ruzicka for incidence data

# changes in B,c and D along the years, occurrence data####

C <- c(tbi_13_14$BCD.summary$`mean(C/den)`, #composantes gain des p?riodes
       tbi_13_15$BCD.summary$`mean(C/den)`,
       tbi_13_16$BCD.summary$`mean(C/den)`,
       tbi_13_17$BCD.summary$`mean(C/den)`,
       tbi_13_18$BCD.summary$`mean(C/den)`,
       tbi_13_19$BCD.summary$`mean(C/den)`)

B <- c(tbi_13_14$BCD.summary$`mean(B/den)`, #composantes pertes des p?riodes
       tbi_13_15$BCD.summary$`mean(B/den)`,
       tbi_13_16$BCD.summary$`mean(B/den)`,
       tbi_13_17$BCD.summary$`mean(B/den)`,
       tbi_13_18$BCD.summary$`mean(B/den)`,
       tbi_13_19$BCD.summary$`mean(B/den)`)

D <- c(tbi_13_14$BCD.summary$`mean(D)`, #dissimlarit?s des p?riodes
       tbi_13_15$BCD.summary$`mean(D)`,
       tbi_13_16$BCD.summary$`mean(D)`,
       tbi_13_17$BCD.summary$`mean(D)`,
       tbi_13_18$BCD.summary$`mean(D)`,
       tbi_13_19$BCD.summary$`mean(D)`)

BCD <- cbind.data.frame(B,C,D) #tableau des composantes et des dissimilarit?s

BCD$period = c("13-14","13-15", "13-16",  "13-17", "13-18", "13-19") #?tiquette des p?riodes...
rownames(BCD) = c("13-14","13-15", "13-16",  "13-17", "13-18", "13-19") #... en nom de ligne

Q1 <- ggplot() + #premi?re couche
  geom_point(data=BCD, aes(x=period, y=B), colour = "red", size=3, shape = 16) + 
  scale_y_continuous(limits= c(0,.6), breaks = seq(0,0.6,0.2)) +
  geom_line(data=BCD, aes(x=period, y=B, group=1), colour = "red")+
  xlab("year of surveys") + ylab("Ruzicka dissimilarity") + theme + theme(axis.title.x = element_blank()) + 
  theme(axis.title.y = element_blank(), axis.text.y=element_text(size =16), axis.text.x=element_text(size =14)) 

Q2 <- ggplot() + # seconde couche
  geom_point(data=BCD, aes(x=period, y=C), colour = "skyblue", size=3, shape=17) + 
  scale_y_continuous(limits= c(0,.6), breaks = seq(0,0.6,0.2)) +
  geom_line(data=BCD, aes(x=period, y=C, group=2), colour = "skyblue") 

Q3 <- ggplot() + #troisi?me couche
  geom_point(data=BCD, aes(x=period, y=D), colour = "black", size=3, shape=15) + 
  scale_y_continuous(limits= c(0,.6), breaks = seq(0,0.6,0.2)) +
  geom_line(data=BCD, aes(x=period, y=D, group=3), colour = "black") 

PLAPLOT = Q1 + Q2$layers[[1]]+ Q2$layers[[2]] + Q3$layers[[1]]+Q3$layers[[2]] #compilation des couches
PLAPLOT


# intraplot dissimilarities----
bef2013<-as.data.frame(df %>% filter(grepl('2013', year)))
bef2013$year <- NULL
fac <- select_if(bef2013,is.factor)
name<-paste(fac$site,fac$placette, sep="_")
bef2013 <- cbind(name, bef2013)
rownames(bef2013) = make.names(name, unique=TRUE)


bef2014<-as.data.frame(df %>% filter(grepl('2014', year)))
bef2014$year <- NULL
fac <- select_if(bef2014,is.factor)
name<- paste(fac$site,fac$placette, sep="_")
bef2014 <- cbind(name, bef2014)
rownames(bef2014) = make.names(name, unique=TRUE)

bef2015<-as.data.frame(df %>% filter(grepl('2015', year)))
bef2015$year <- NULL
fac <- select_if(bef2015,is.factor)
name<- paste(fac$site,fac$placette, sep="_")
bef2015 <- cbind(name, bef2015)
rownames(bef2015) = make.names(name, unique=TRUE)

bef2016<-as.data.frame(df %>% filter(grepl('2016', year)))
bef2016$year <- NULL
fac <- select_if(bef2016,is.factor)
name<- paste(fac$site,fac$placette, sep="_")
bef2016 <- cbind(name, bef2016)
rownames(bef2016) = make.names(name, unique=TRUE)

bef2017<-as.data.frame(df %>% filter(grepl('2017', year)))
bef2017$year <- NULL
fac <- select_if(bef2017,is.factor)
name<- paste(fac$site,fac$placette, sep="_")
bef2017 <- cbind(name, bef2017)
rownames(bef2017) = make.names(name, unique=TRUE)

bef2018<-as.data.frame(df %>% filter(grepl('2018', year)))
bef2018$year <- NULL
fac <- select_if(bef2018,is.factor)
name<- paste(fac$site,fac$placette, sep="_")
bef2018 <- cbind(name, bef2018)
rownames(bef2018) = make.names(name, unique=TRUE)

bef2019<-as.data.frame(df %>% filter(grepl('2019', year)))
bef2019$year <- NULL
fac <- select_if(bef2019,is.factor)
name<- paste(fac$site,fac$placette, sep="_")
bef2019 <- cbind(name, bef2019)
rownames(bef2019) = make.names(name, unique=TRUE)
###TBI P1/P2

##TBI between 2013 entre P1 et P2
P12013<-bef2013 %>% filter(grepl('1', placette))
P22013<-bef2013 %>% filter(grepl('2', placette))
P12013 <- subset(P12013, site %in% intersect(P22013$site, P12013$site))
P22013 <- subset(P22013, site %in% intersect(P12013$site, P22013$site))
m20132 <- as.matrix(P22013[,4:99])
m20131 <- as.matrix(P12013[,4:99])
tbi_2013P1P2 <- TBI(m20131,m20132, method = "ruzicka",pa.tr = F, nperm = 99,BCD = TRUE, save.BC = T,test.BC = TRUE, test.t.perm = T, seed. =   )

##TBI between 2014 entre P1 et P2
P12014<-bef2014 %>% filter(grepl('1', placette))
P22014<-bef2014 %>% filter(grepl('2', placette))
P12014 <- subset(P12014, site %in% intersect(P22014$site, P12014$site))
P22014 <- subset(P22014, site %in% intersect(P12014$site, P22014$site))
m20142 <- as.matrix(P22014[,4:99])
m20141 <- as.matrix(P12014[,4:99])
tbi_2014P1P2 <- TBI(m20141,m20142, method = "ruzicka",pa.tr = F, nperm = 99,BCD = TRUE, save.BC = T,test.BC = TRUE, test.t.perm = T, seed. =   )


##TBI between 2015 entre P1 et P2
P12015<-bef2015 %>% filter(grepl('1', placette))
P22015<-bef2015 %>% filter(grepl('2', placette))
P12015 <- subset(P12015, site %in% intersect(P22015$site, P12015$site))
P22015 <- subset(P22015, site %in% intersect(P12015$site, P22015$site))
m20152 <- as.matrix(P22015[,4:99])
m20151 <- as.matrix(P12015[,4:99])
tbi_2015P1P2 <- TBI(m20151,m20152, method = "ruzicka",pa.tr = F, nperm = 99,BCD = TRUE, save.BC = T,test.BC = TRUE, test.t.perm = T, seed. =   )


##TBI between 2016 entre P1 et P2
P12016<-bef2016 %>% filter(grepl('1', placette))
P22016<-bef2016 %>% filter(grepl('2', placette))
P12016 <- subset(P12016, site %in% intersect(P22016$site, P12016$site))
P22016 <- subset(P22016, site %in% intersect(P12016$site, P22016$site))
m20162 <- as.matrix(P22016[,4:99])
m20161 <- as.matrix(P12016[,4:99])
tbi_2016P1P2 <- TBI(m20161,m20162, method = "ruzicka",pa.tr = F, nperm = 99,BCD = TRUE, save.BC = T,test.BC = TRUE, test.t.perm = T, seed. =   )


##TBI between 2017 entre P1 et P2
P12017<-bef2017 %>% filter(grepl('1', placette))
P22017<-bef2017 %>% filter(grepl('2', placette))
P12017 <- subset(P12017, site %in% intersect(P22017$site, P12017$site))
P22017 <- subset(P22017, site %in% intersect(P12017$site, P22017$site))
m20172 <- as.matrix(P22017[,4:99])
m20171 <- as.matrix(P12017[,4:99])
tbi_2017P1P2 <- TBI(m20171,m20172, method = "ruzicka",pa.tr = F, nperm = 99,BCD = TRUE, save.BC = T,test.BC = TRUE, test.t.perm = T, seed. =   )

##TBI between 2018 entre P1 et P2
P12018<-bef2018 %>% filter(grepl('1', placette))
P22018<-bef2018 %>% filter(grepl('2', placette))
P12018 <- subset(P12018, site %in% intersect(P22018$site, P12018$site))
P22018 <- subset(P22018, site %in% intersect(P12018$site, P22018$site))
m20182 <- as.matrix(P22018[,4:99])
m20181 <- as.matrix(P12018[,4:99])
tbi_2018P1P2 <- TBI(m20181,m20182, method = "ruzicka",pa.tr = F, nperm = 99,BCD = TRUE, save.BC = T,test.BC = TRUE, test.t.perm = T, seed. =   )

##TBI between 2019 entre P1 et P2
P12019<-bef2019 %>% filter(grepl('1', placette))
P22019<-bef2019 %>% filter(grepl('2', placette))
P12019 <- subset(P12019, site %in% intersect(P22019$site, P12019$site))
P22019 <- subset(P22019, site %in% intersect(P12019$site, P22019$site))
m20192 <- as.matrix(P22019[,4:99])
m20191 <- as.matrix(P12019[,4:99])
tbi_2019P1P2 <- TBI(m20191,m20192, method = "ruzicka",pa.tr = F, nperm = 99,BCD = TRUE, save.BC = T,test.BC = TRUE, test.t.perm = T, seed. =   )


ARBP1P2 <- c(tbi_2014P1P2$BCD.mat$`D=(B+C)/(A+B+C)`[1],
             tbi_2015P1P2$BCD.mat$`D=(B+C)/(A+B+C)`[1],
             tbi_2016P1P2$BCD.mat$`D=(B+C)/(A+B+C)`[1],
             tbi_2017P1P2$BCD.mat$`D=(B+C)/(A+B+C)`[1],
             tbi_2018P1P2$BCD.mat$`D=(B+C)/(A+B+C)`[1],
             tbi_2019P1P2$BCD.mat$`D=(B+C)/(A+B+C)`[1])

CASEP1P2 <- c(tbi_2013P1P2$BCD.mat$`D=(B+C)/(A+B+C)`[1], 
              tbi_2014P1P2$BCD.mat$`D=(B+C)/(A+B+C)`[2],
              tbi_2015P1P2$BCD.mat$`D=(B+C)/(A+B+C)`[2],
              tbi_2016P1P2$BCD.mat$`D=(B+C)/(A+B+C)`[2],
              tbi_2017P1P2$BCD.mat$`D=(B+C)/(A+B+C)`[2],
              tbi_2018P1P2$BCD.mat$`D=(B+C)/(A+B+C)`[2],
              tbi_2019P1P2$BCD.mat$`D=(B+C)/(A+B+C)`[2])

CATP1P2 <- c(tbi_2013P1P2$BCD.mat$`D=(B+C)/(A+B+C)`[2], 
             tbi_2014P1P2$BCD.mat$`D=(B+C)/(A+B+C)`[3],
             tbi_2015P1P2$BCD.mat$`D=(B+C)/(A+B+C)`[3],
             tbi_2016P1P2$BCD.mat$`D=(B+C)/(A+B+C)`[3],
             tbi_2017P1P2$BCD.mat$`D=(B+C)/(A+B+C)`[3],
             tbi_2018P1P2$BCD.mat$`D=(B+C)/(A+B+C)`[3],
             tbi_2019P1P2$BCD.mat$`D=(B+C)/(A+B+C)`[3])

PLAP1P2 <- c(tbi_2013P1P2$BCD.mat$`D=(B+C)/(A+B+C)`[3], 
             tbi_2014P1P2$BCD.mat$`D=(B+C)/(A+B+C)`[4],
             tbi_2015P1P2$BCD.mat$`D=(B+C)/(A+B+C)`[4],
             tbi_2016P1P2$BCD.mat$`D=(B+C)/(A+B+C)`[4],
             tbi_2017P1P2$BCD.mat$`D=(B+C)/(A+B+C)`[4],
             tbi_2018P1P2$BCD.mat$`D=(B+C)/(A+B+C)`[4],
             tbi_2019P1P2$BCD.mat$`D=(B+C)/(A+B+C)`[4])

RATP1P2 <- c(tbi_2013P1P2$BCD.mat$`D=(B+C)/(A+B+C)`[5],
             tbi_2016P1P2$BCD.mat$`D=(B+C)/(A+B+C)`[6],
             tbi_2017P1P2$BCD.mat$`D=(B+C)/(A+B+C)`[6],
             tbi_2018P1P2$BCD.mat$`D=(B+C)/(A+B+C)`[6],
             tbi_2019P1P2$BCD.mat$`D=(B+C)/(A+B+C)`[5])

ULLP1P2 <- c(tbi_2013P1P2$BCD.mat$`D=(B+C)/(A+B+C)`[6],
             tbi_2016P1P2$BCD.mat$`D=(B+C)/(A+B+C)`[7],
             tbi_2017P1P2$BCD.mat$`D=(B+C)/(A+B+C)`[7],
             tbi_2018P1P2$BCD.mat$`D=(B+C)/(A+B+C)`[7],
             tbi_2019P1P2$BCD.mat$`D=(B+C)/(A+B+C)`[6])

###TBI P2/P3

##TBI between 2013 entre P2 et P3
P22013<-bef2013 %>% filter(grepl('2', placette))
P32013<-bef2013 %>% filter(grepl('3', placette))
P22013 <- subset(P22013, site %in% intersect(P32013$site, P22013$site))
P32013 <- subset(P32013, site %in% intersect(P22013$site, P32013$site))
m20132 <- as.matrix(P32013[,4:99])
m20131 <- as.matrix(P22013[,4:99])
tbi_2013P2P3 <- TBI(m20131,m20132, method = "ruzicka",pa.tr = F, nperm = 99,BCD = TRUE, save.BC = T,test.BC = TRUE, test.t.perm = T, seed. =   )

##TBI between 2014 entre P2 et P3
P22014<-bef2014 %>% filter(grepl('2', placette))
P32014<-bef2014 %>% filter(grepl('3', placette))
P22014 <- subset(P22014, site %in% intersect(P32014$site, P22014$site))
P32014 <- subset(P32014, site %in% intersect(P22014$site, P32014$site))
m20142 <- as.matrix(P32014[,4:99])
m20141 <- as.matrix(P22014[,4:99])
tbi_2014P2P3 <- TBI(m20141,m20142, method = "ruzicka",pa.tr = F, nperm = 99,BCD = TRUE, save.BC = T,test.BC = TRUE, test.t.perm = T, seed. =   )


##TBI between 2015 entre P2 et P3
P22015<-bef2015 %>% filter(grepl('2', placette))
P32015<-bef2015 %>% filter(grepl('3', placette))
P22015 <- subset(P22015, site %in% intersect(P32015$site, P22015$site))
P32015 <- subset(P32015, site %in% intersect(P22015$site, P32015$site))
m20152 <- as.matrix(P32015[,4:99])
m20151 <- as.matrix(P22015[,4:99])
tbi_2015P2P3 <- TBI(m20151,m20152, method = "ruzicka",pa.tr = F, nperm = 99,BCD = TRUE, save.BC = T,test.BC = TRUE, test.t.perm = T, seed. =   )


##TBI between 2016 entre P2 et P3
P22016<-bef2016 %>% filter(grepl('2', placette))
P32016<-bef2016 %>% filter(grepl('3', placette))
P22016 <- subset(P22016, site %in% intersect(P32016$site, P22016$site))
P32016 <- subset(P32016, site %in% intersect(P22016$site, P32016$site))
m20162 <- as.matrix(P32016[,4:99])
m20161 <- as.matrix(P22016[,4:99])
tbi_2016P2P3 <- TBI(m20161,m20162, method = "ruzicka",pa.tr = F, nperm = 99,BCD = TRUE, save.BC = T,test.BC = TRUE, test.t.perm = T, seed. =   )


##TBI between 2017 entre P2 et P3
P22017<-bef2017 %>% filter(grepl('2', placette))
P32017<-bef2017 %>% filter(grepl('3', placette))
P22017 <- subset(P22017, site %in% intersect(P32017$site, P22017$site))
P32017 <- subset(P32017, site %in% intersect(P22017$site, P32017$site))
m20172 <- as.matrix(P32017[,4:99])
m20171 <- as.matrix(P22017[,4:99])
tbi_2017P2P3 <- TBI(m20171,m20172, method = "ruzicka",pa.tr = F, nperm = 99,BCD = TRUE, save.BC = T,test.BC = TRUE, test.t.perm = T, seed. =   )

##TBI between 2018 entre P2 et P3
P22018<-bef2018 %>% filter(grepl('2', placette))
P32018<-bef2018 %>% filter(grepl('3', placette))
P22018 <- subset(P22018, site %in% intersect(P32018$site, P22018$site))
P32018 <- subset(P32018, site %in% intersect(P22018$site, P32018$site))
m20182 <- as.matrix(P32018[,4:99])
m20181 <- as.matrix(P22018[,4:99])
tbi_2018P2P3 <- TBI(m20181,m20182, method = "ruzicka",pa.tr = F, nperm = 99,BCD = TRUE, save.BC = T,test.BC = TRUE, test.t.perm = T, seed. =   )

##TBI between 2019 entre P2 et P3
P22019<-bef2019 %>% filter(grepl('2', placette))
P32019<-bef2019 %>% filter(grepl('3', placette))
P22019 <- subset(P22019, site %in% intersect(P32019$site, P22019$site))
P32019 <- subset(P32019, site %in% intersect(P22019$site, P32019$site))
m20192 <- as.matrix(P32019[,4:99])
m20191 <- as.matrix(P22019[,4:99])
tbi_2019P2P3 <- TBI(m20191,m20192, method = "ruzicka",pa.tr = F, nperm = 99,BCD = TRUE, save.BC = T,test.BC = TRUE, test.t.perm = T, seed. =   )



ARBP2P3 <- c(tbi_2015P2P3$BCD.mat$`D=(B+C)/(A+B+C)`[1],
             tbi_2017P2P3$BCD.mat$`D=(B+C)/(A+B+C)`[1],
             tbi_2018P2P3$BCD.mat$`D=(B+C)/(A+B+C)`[1])

CASEP2P3 <- c(tbi_2013P2P3$BCD.mat$`D=(B+C)/(A+B+C)`[1], 
              tbi_2014P2P3$BCD.mat$`D=(B+C)/(A+B+C)`[1],
              tbi_2015P2P3$BCD.mat$`D=(B+C)/(A+B+C)`[2],
              tbi_2016P2P3$BCD.mat$`D=(B+C)/(A+B+C)`[1],
              tbi_2017P2P3$BCD.mat$`D=(B+C)/(A+B+C)`[2],
              tbi_2018P2P3$BCD.mat$`D=(B+C)/(A+B+C)`[2],
              tbi_2019P2P3$BCD.mat$`D=(B+C)/(A+B+C)`[1])

CATP2P3 <- c(tbi_2013P2P3$BCD.mat$`D=(B+C)/(A+B+C)`[2], 
             tbi_2014P2P3$BCD.mat$`D=(B+C)/(A+B+C)`[2],
             tbi_2015P2P3$BCD.mat$`D=(B+C)/(A+B+C)`[3],
             tbi_2016P2P3$BCD.mat$`D=(B+C)/(A+B+C)`[2],
             tbi_2017P2P3$BCD.mat$`D=(B+C)/(A+B+C)`[3],
             tbi_2018P2P3$BCD.mat$`D=(B+C)/(A+B+C)`[3],
             tbi_2019P2P3$BCD.mat$`D=(B+C)/(A+B+C)`[2])

PLAP2P3 <- c(tbi_2013P2P3$BCD.mat$`D=(B+C)/(A+B+C)`[4], 
             tbi_2014P2P3$BCD.mat$`D=(B+C)/(A+B+C)`[3],
             tbi_2015P2P3$BCD.mat$`D=(B+C)/(A+B+C)`[4],
             tbi_2016P2P3$BCD.mat$`D=(B+C)/(A+B+C)`[4],
             tbi_2017P2P3$BCD.mat$`D=(B+C)/(A+B+C)`[5],
             tbi_2018P2P3$BCD.mat$`D=(B+C)/(A+B+C)`[5],
             tbi_2019P2P3$BCD.mat$`D=(B+C)/(A+B+C)`[4])

CREP2P3 <- c(tbi_2013P2P3$BCD.mat$`D=(B+C)/(A+B+C)`[3],
             tbi_2016P2P3$BCD.mat$`D=(B+C)/(A+B+C)`[3],
             tbi_2017P2P3$BCD.mat$`D=(B+C)/(A+B+C)`[4],
             tbi_2018P2P3$BCD.mat$`D=(B+C)/(A+B+C)`[4],
             tbi_2019P2P3$BCD.mat$`D=(B+C)/(A+B+C)`[3])

RATP2P3 <- c(tbi_2013P2P3$BCD.mat$`D=(B+C)/(A+B+C)`[6],
             tbi_2016P2P3$BCD.mat$`D=(B+C)/(A+B+C)`[6],
             tbi_2017P2P3$BCD.mat$`D=(B+C)/(A+B+C)`[7],
             tbi_2018P2P3$BCD.mat$`D=(B+C)/(A+B+C)`[7],
             tbi_2019P2P3$BCD.mat$`D=(B+C)/(A+B+C)`[5])

ULLP2P3 <- c(tbi_2013P2P3$BCD.mat$`D=(B+C)/(A+B+C)`[7],
             tbi_2016P2P3$BCD.mat$`D=(B+C)/(A+B+C)`[7],
             tbi_2017P2P3$BCD.mat$`D=(B+C)/(A+B+C)`[8],
             tbi_2018P2P3$BCD.mat$`D=(B+C)/(A+B+C)`[8],
             tbi_2019P2P3$BCD.mat$`D=(B+C)/(A+B+C)`[6])


###TBI P1/P3

P12013<-bef2013 %>% filter(grepl('1', placette))
P32013<-bef2013 %>% filter(grepl('3', placette))
P12013 <- subset(P12013, site %in% intersect(P32013$site, P12013$site))
P32013 <- subset(P32013, site %in% intersect(P12013$site, P32013$site))
m20132 <- as.matrix(P32013[,4:99])
m20131 <- as.matrix(P12013[,4:99])
tbi_2013P1P3 <- TBI(m20131,m20132, method = "ruzicka",pa.tr = F, nperm = 99,BCD = TRUE, save.BC = T,test.BC = TRUE, test.t.perm = T, seed. =   )

##TBI between 2014 entre P1 et P3
P12014<-bef2014 %>% filter(grepl('1', placette))
P32014<-bef2014 %>% filter(grepl('3', placette))
P12014 <- subset(P12014, site %in% intersect(P32014$site, P12014$site))
P32014 <- subset(P32014, site %in% intersect(P12014$site, P32014$site))
m20142 <- as.matrix(P32014[,4:99])
m20141 <- as.matrix(P12014[,4:99])
tbi_2014P1P3 <- TBI(m20141,m20142, method = "ruzicka",pa.tr = F, nperm = 99,BCD = TRUE, save.BC = T,test.BC = TRUE, test.t.perm = T, seed. =   )


##TBI between 2015 entre P1 et P3
P12015<-bef2015 %>% filter(grepl('1', placette))
P32015<-bef2015 %>% filter(grepl('3', placette))
P12015 <- subset(P12015, site %in% intersect(P32015$site, P12015$site))
P32015 <- subset(P32015, site %in% intersect(P12015$site, P32015$site))
m20152 <- as.matrix(P32015[,4:99])
m20151 <- as.matrix(P12015[,4:99])
tbi_2015P1P3 <- TBI(m20151,m20152, method = "ruzicka",pa.tr = F, nperm = 99,BCD = TRUE, save.BC = T,test.BC = TRUE, test.t.perm = T, seed. =   )


##TBI between 2016 entre P1 et P3
P12016<-bef2016 %>% filter(grepl('1', placette))
P32016<-bef2016 %>% filter(grepl('3', placette))
P12016 <- subset(P12016, site %in% intersect(P32016$site, P12016$site))
P32016 <- subset(P32016, site %in% intersect(P12016$site, P32016$site))
m20162 <- as.matrix(P32016[,4:99])
m20161 <- as.matrix(P12016[,4:99])
tbi_2016P1P3 <- TBI(m20161,m20162, method = "ruzicka",pa.tr = F, nperm = 99,BCD = TRUE, save.BC = T,test.BC = TRUE, test.t.perm = T, seed. =   )


##TBI between 2017 entre P1 et P3
P12017<-bef2017 %>% filter(grepl('1', placette))
P32017<-bef2017 %>% filter(grepl('3', placette))
P12017 <- subset(P12017, site %in% intersect(P32017$site, P12017$site))
P32017 <- subset(P32017, site %in% intersect(P12017$site, P32017$site))
m20172 <- as.matrix(P32017[,4:99])
m20171 <- as.matrix(P12017[,4:99])
tbi_2017P1P3 <- TBI(m20171,m20172, method = "ruzicka",pa.tr = F, nperm = 99,BCD = TRUE, save.BC = T,test.BC = TRUE, test.t.perm = T, seed. =   )

##TBI between 2018 entre P1 et P3
P12018<-bef2018 %>% filter(grepl('1', placette))
P32018<-bef2018 %>% filter(grepl('3', placette))
P12018 <- subset(P12018, site %in% intersect(P32018$site, P12018$site))
P32018 <- subset(P32018, site %in% intersect(P12018$site, P32018$site))
m20182 <- as.matrix(P32018[,4:99])
m20181 <- as.matrix(P12018[,4:99])
tbi_2018P1P3 <- TBI(m20181,m20182, method = "ruzicka",pa.tr = F, nperm = 99,BCD = TRUE, save.BC = T,test.BC = TRUE, test.t.perm = T, seed. =   )

##TBI between 2019 entre P1 et P3
P12019<-bef2019 %>% filter(grepl('1', placette))
P32019<-bef2019 %>% filter(grepl('3', placette))
P12019 <- subset(P12019, site %in% intersect(P32019$site, P12019$site))
P32019 <- subset(P32019, site %in% intersect(P12019$site, P32019$site))
m20192 <- as.matrix(P32019[,4:99])
m20191 <- as.matrix(P12019[,4:99])
tbi_2019P1P3 <- TBI(m20191,m20192, method = "ruzicka",pa.tr = F, nperm = 99,BCD = TRUE, save.BC = T,test.BC = TRUE, test.t.perm = T, seed. =   )


ARBP1P3 <- c(tbi_2015P1P3$BCD.mat$`D=(B+C)/(A+B+C)`[1],
             tbi_2017P1P3$BCD.mat$`D=(B+C)/(A+B+C)`[1],
             tbi_2018P1P3$BCD.mat$`D=(B+C)/(A+B+C)`[1])

CASEP1P3 <- c(tbi_2013P1P3$BCD.mat$`D=(B+C)/(A+B+C)`[1], 
              tbi_2014P1P3$BCD.mat$`D=(B+C)/(A+B+C)`[1],
              tbi_2015P1P3$BCD.mat$`D=(B+C)/(A+B+C)`[2],
              tbi_2016P1P3$BCD.mat$`D=(B+C)/(A+B+C)`[1],
              tbi_2017P1P3$BCD.mat$`D=(B+C)/(A+B+C)`[2],
              tbi_2018P1P3$BCD.mat$`D=(B+C)/(A+B+C)`[2],
              tbi_2019P1P3$BCD.mat$`D=(B+C)/(A+B+C)`[1])

CATP1P3 <- c(tbi_2013P1P3$BCD.mat$`D=(B+C)/(A+B+C)`[2], 
             tbi_2014P1P3$BCD.mat$`D=(B+C)/(A+B+C)`[2],
             tbi_2015P1P3$BCD.mat$`D=(B+C)/(A+B+C)`[3],
             tbi_2016P1P3$BCD.mat$`D=(B+C)/(A+B+C)`[2],
             tbi_2017P1P3$BCD.mat$`D=(B+C)/(A+B+C)`[3],
             tbi_2018P1P3$BCD.mat$`D=(B+C)/(A+B+C)`[3],
             tbi_2019P1P3$BCD.mat$`D=(B+C)/(A+B+C)`[2])

PLAP1P3 <- c(tbi_2013P1P3$BCD.mat$`D=(B+C)/(A+B+C)`[3], 
             tbi_2014P1P3$BCD.mat$`D=(B+C)/(A+B+C)`[3],
             tbi_2015P1P3$BCD.mat$`D=(B+C)/(A+B+C)`[4],
             tbi_2016P1P3$BCD.mat$`D=(B+C)/(A+B+C)`[3],
             tbi_2017P1P3$BCD.mat$`D=(B+C)/(A+B+C)`[4],
             tbi_2018P1P3$BCD.mat$`D=(B+C)/(A+B+C)`[4],
             tbi_2019P1P3$BCD.mat$`D=(B+C)/(A+B+C)`[3])

RATP1P3 <- c(tbi_2013P1P3$BCD.mat$`D=(B+C)/(A+B+C)`[5],
             tbi_2016P1P3$BCD.mat$`D=(B+C)/(A+B+C)`[5],
             tbi_2017P1P3$BCD.mat$`D=(B+C)/(A+B+C)`[6],
             tbi_2018P1P3$BCD.mat$`D=(B+C)/(A+B+C)`[6],
             tbi_2019P1P3$BCD.mat$`D=(B+C)/(A+B+C)`[4])

ULLP1P3 <- c(tbi_2013P1P3$BCD.mat$`D=(B+C)/(A+B+C)`[6],
             tbi_2016P1P3$BCD.mat$`D=(B+C)/(A+B+C)`[6],
             tbi_2017P1P3$BCD.mat$`D=(B+C)/(A+B+C)`[7],
             tbi_2018P1P3$BCD.mat$`D=(B+C)/(A+B+C)`[7],
             tbi_2019P1P3$BCD.mat$`D=(B+C)/(A+B+C)`[5])


theme<-theme(panel.background = element_blank(),
             panel.border=element_rect(fill=NA),
             panel.grid.major = element_blank(),
             panel.grid.minor = element_blank(),
             strip.background=element_blank(),
             axis.text.x=element_text(colour="black"),
             axis.text.y=element_text(colour="black"),
             axis.ticks=element_line(colour="black"),
             plot.margin=unit(c(1,1,1,1),"line"))


###CASE SITUATION
TBICASE <- cbind.data.frame(CASEP1P2,CASEP1P3, CASEP2P3)
TBICASE$period = c("2013","2014","2015","2016","2017","2018","2019")
rownames(TBICASE) =c("2013","2014","2015","2016","2017","2018","2019")

Q1 <- ggplot() + 
  geom_point(data=TBICASE, aes(x=period, y=CASEP1P2), colour = "#d1495b", size=3, shape = 16) +  theme +
  scale_y_continuous(limits= c(0,1), breaks = seq(0,1,.5)) +
  geom_line(data=TBICASE, aes(x=period, y=CASEP1P2, group=1), colour = "#d1495b")+
  theme(axis.title.x = element_blank())+
  theme(axis.title.y = element_blank(), axis.text.y=element_text(size =16), axis.text.x=element_text(size =16)) 

Q2 <- ggplot() +
  geom_point(data=TBICASE, aes(x=period, y=CASEP1P3), colour = "#edae49", size=3, shape=17) + 
  scale_y_continuous(limits= c(0,1),  breaks = seq(0,1,.5)) +
  geom_line(data=TBICASE, aes(x=period, y=CASEP1P3, group=2), colour = "#edae49") 

Q3 <- ggplot() +
  geom_point(data=TBICASE, aes(x=period, y=CASEP2P3), colour = "#66a182", size=3, shape=15) + 
  scale_y_continuous(limits= c(0,1),  breaks = seq(0,1,.5)) +
  geom_line(data=TBICASE, aes(x=period, y=CASEP2P3, group=3), colour = "#66a182") 

CASE_p = Q1 + Q2$layers[[1]]+ Q2$layers[[2]] + Q3$layers[[1]]+Q3$layers[[2]]

###CATAPERDIS SITUATION
TBICAT <- cbind.data.frame(CATP1P2,CATP1P3, CATP2P3)
TBICAT$period = c("2013","2014","2015","2016","2017","2018","2019")
rownames(TBICAT) =c("2013","2014","2015","2016","2017","2018","2019")

Q1 <- ggplot() + 
  geom_point(data=TBICAT, aes(x=period, y=CATP1P2), colour = "#d1495b", size=3, shape = 16) + theme +
  scale_y_continuous(limits= c(0,1), breaks = seq(0,1,.5)) +
  geom_line(data=TBICAT, aes(x=period, y=CATP1P2, group=1), colour = "#d1495b")+
  theme(axis.title.x = element_blank())+ 
  theme(axis.title.y = element_blank(), axis.text.y=element_text(size =16), axis.text.x=element_text(size =16)) 

Q2 <- ggplot() +
  geom_point(data=TBICAT, aes(x=period, y=CATP1P3), colour = "#edae49", size=3, shape=17) + 
  scale_y_continuous(limits= c(0,1),  breaks = seq(0,1,.5)) +
  geom_line(data=TBICAT, aes(x=period, y=CATP1P3, group=2), colour = "#edae49") 

Q3 <- ggplot() +
  geom_point(data=TBICAT, aes(x=period, y=CATP2P3), colour = "#66a182", size=3, shape=15) + 
  scale_y_continuous(limits= c(0,1),  breaks = seq(0,1,.5)) +
  geom_line(data=TBICAT, aes(x=period, y=CATP2P3, group=3), colour = "#66a182") 

CAT_p = Q1 + Q2$layers[[1]]+ Q2$layers[[2]] + Q3$layers[[1]]+Q3$layers[[2]]

###CREUSSANS SITUATION
TBICRE <- cbind.data.frame(CREP1P2,CREP1P3, CREP2P3)
TBICRE$period = c("2013","2016","2017","2018","2019")
rownames(TBICRE) = c("2013","2016","2017","2018","2019")

Q1 <- ggplot() + 
  geom_point(data=TBICRE, aes(x=period, y=CREP1P2), colour = "#d1495b", size=3, shape = 16) +  theme +
  scale_y_continuous(limits= c(0,1), breaks = seq(0,1,.5)) +
  geom_line(data=TBICRE, aes(x=period, y=CREP1P2, group=1), colour = "#d1495b")+
  theme(axis.title.x = element_blank())+
  theme(axis.title.y = element_blank(), axis.text.y=element_text(size =16), axis.text.x=element_text(size =16)) 

Q2 <- ggplot() +
  geom_point(data=TBICRE, aes(x=period, y=CREP1P3), colour = "#edae49", size=3, shape=17) + 
  scale_y_continuous(limits= c(0,1),  breaks = seq(0,1,.5)) +
  geom_line(data=TBICRE, aes(x=period, y=CREP1P3, group=2), colour = "#edae49") 

Q3 <- ggplot() +
  geom_point(data=TBICRE, aes(x=period, y=CREP2P3), colour = "#66a182", size=3, shape=15) + 
  scale_y_continuous(limits= c(0,1),  breaks = seq(0,1,.5)) +
  geom_line(data=TBICRE, aes(x=period, y=CREP2P3, group=3), colour = "#66a182") 

CRE_p = Q1 + Q2$layers[[1]]+ Q2$layers[[2]] + Q3$layers[[1]]+Q3$layers[[2]]

###PLANES SITUATION
TBIPLA <- cbind.data.frame(PLAP1P2,PLAP1P3, PLAP2P3)
TBIPLA$period = c("2013","2014","2015","2016","2017","2018","2019")
rownames(TBIPLA) =c("2013","2014","2015","2016","2017","2018","2019")

Q1 <- ggplot() + 
  geom_point(data=TBIPLA, aes(x=period, y=PLAP1P2), colour = "#d1495b", size=3, shape = 16) +  theme +
  scale_y_continuous(limits= c(0,1), breaks = seq(0,1,.5)) +
  geom_line(data=TBIPLA, aes(x=period, y=PLAP1P2, group=1), colour = "#d1495b")+
  theme(axis.title.x = element_blank())+
  theme(axis.title.y = element_blank(), axis.text.y=element_text(size =16), axis.text.x=element_text(size =16)) 

Q2 <- ggplot() +
  geom_point(data=TBIPLA, aes(x=period, y=PLAP1P3), colour = "#edae49", size=3, shape=17) + 
  scale_y_continuous(limits= c(0,1),  breaks = seq(0,1,.5)) +
  geom_line(data=TBIPLA, aes(x=period, y=PLAP1P3, group=2), colour = "#edae49") 

Q3 <- ggplot() +
  geom_point(data=TBIPLA, aes(x=period, y=PLAP2P3), colour = "#66a182", size=3, shape=15) + 
  scale_y_continuous(limits= c(0,1),  breaks = seq(0,1,.5)) +
  geom_line(data=TBIPLA, aes(x=period, y=PLAP2P3, group=3), colour = "#66a182") 

PLA_p = Q1 + Q2$layers[[1]]+ Q2$layers[[2]] + Q3$layers[[1]]+Q3$layers[[2]]

###RATERA SITUATION
TBIRAT <- cbind.data.frame(RATP1P2,RATP1P3, RATP2P3)
TBIRAT$period = c("2013","2016","2017","2018","2019")
rownames(TBIRAT) = c("2013","2016","2017","2018","2019")

Q1 <- ggplot() + 
  geom_point(data=TBIRAT, aes(x=period, y=RATP1P2), colour = "#d1495b", size=3, shape = 16) +  theme +
  scale_y_continuous(limits= c(0,1), breaks = seq(0,1,.5)) +
  geom_line(data=TBIRAT, aes(x=period, y=RATP1P2, group=1), colour = "#d1495b")+
  theme(axis.title.x = element_blank())+
  theme(axis.title.y = element_blank(), axis.text.y=element_text(size =16), axis.text.x=element_text(size =16)) 

Q2 <- ggplot() +
  geom_point(data=TBIRAT, aes(x=period, y=RATP1P3), colour = "#edae49", size=3, shape=17) + 
  scale_y_continuous(limits= c(0,1),  breaks = seq(0,1,.5)) +
  geom_line(data=TBIRAT, aes(x=period, y=RATP1P3, group=2), colour = "#edae49") 

Q3 <- ggplot() +
  geom_point(data=TBIRAT, aes(x=period, y=RATP2P3), colour = "#66a182", size=3, shape=15) + 
  scale_y_continuous(limits= c(0,1),  breaks = seq(0,1,.5)) +
  geom_line(data=TBIRAT, aes(x=period, y=RATP2P3, group=3), colour = "#66a182") 

RAT_p = Q1 + Q2$layers[[1]]+ Q2$layers[[2]] + Q3$layers[[1]]+Q3$layers[[2]]

###Ulldeter SITUATION
TBIULL <- cbind.data.frame(ULLP1P2,ULLP1P3, ULLP2P3)
TBIULL$period = c("2013","2016","2017","2018","2019")
rownames(TBIULL) = c("2013","2016","2017","2018","2019")

Q1 <- ggplot() + 
  geom_point(data=TBIULL, aes(x=period, y=ULLP1P2), colour = "#d1495b", size=3, shape = 16) +  theme +
  scale_y_continuous(limits= c(0,1), breaks = seq(0,1,.5)) +
  geom_line(data=TBIULL, aes(x=period, y=ULLP1P2, group=1), colour = "#d1495b")+
  theme(axis.title.x = element_blank())+
  theme(axis.title.y = element_blank(), axis.text.y=element_text(size =16), axis.text.x=element_text(size =16)) 

Q2 <- ggplot() +
  geom_point(data=TBIULL, aes(x=period, y=ULLP1P3), colour = "#edae49", size=3, shape=17) + 
  scale_y_continuous(limits= c(0,1),  breaks = seq(0,1,.5)) +
  geom_line(data=TBIULL, aes(x=period, y=ULLP1P3, group=2), colour = "#edae49") 

Q3 <- ggplot() +
  geom_point(data=TBIULL, aes(x=period, y=ULLP2P3), colour = "#66a182", size=3, shape=15) + 
  scale_y_continuous(limits= c(0,1),  breaks = seq(0,1,.5)) +
  geom_line(data=TBIULL, aes(x=period, y=ULLP2P3, group=3), colour = "#66a182") 

ULL_p = Q1 + Q2$layers[[1]]+ Q2$layers[[2]] + Q3$layers[[1]]+Q3$layers[[2]]


plot_grid(CASE_p,CASEPLOT, CAT_p, CATPLOT, CRE_p, CREPLOT,  PLA_p, PLAPLOT, RAT_p, RATPLOT,  ULL_p, ULLPLOT,  nrow = 6, ncol = 2) ## graphique g??nral des intraplots

