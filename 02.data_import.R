#######################################################
# Project : Analysis of snowbed vegetation dynamics
# Script : 02.data_import
# Load data
# Authors : Guillaume Papuga & Thomas Masclaux
# Last update : 29 april 2021
#######################################################

#####
# 0. Preliminary loading
#####
# taxonomic referential
ref.taxo = read.csv (here::here("data", "raw", "ref.taxo2.csv"), 
                     head = T, row.names = NULL, sep = ";", dec = ",")
full.sp.list = unique(ref.taxo$nom.short)

#####
# 1. CBN Med
#####

# raw floristic datasets
cbnmed = read.csv (here::here("data", "raw", "cbnmed", "planes.pdlc3.csv"), 
                   head = T, row.names = NULL, sep = ";", dec = ",", 
                   check.names = FALSE) # carefull : conserve the full species name to concatenate with the referential (column name can be invalid)

####
# We used the melt-cast process from reshape2 to transform the original dataset into its definitive form
####
# melt the dataset
varnames = colnames(cbnmed)
cbnmed.melt = melt(cbnmed, id=varnames[1:8]) # melt the dataset : the value (pheno & abundance) is in the "value" column

# change the phenology to binary data
cbnmed.melt = cbnmed.melt %>%
  mutate(value = replace(value, value != "",1)) %>% # transform the column "values" to build a 1-0 metric 
  mutate(value = as.numeric (value)) %>% # transform the column to numeric
  replace_na(list(value = 0)) # transform NA to 0

# change species name for the correct name (nom.short) in the referential
ref.taxo.cbnmed = ref.taxo[which(ref.taxo$source=="cbnmed"),] # create a sub-referential containing only the names linked to the dataset

cbnmed.melt = cbnmed.melt %>%
  left_join(ref.taxo.cbnmed [,c("nom.saisi", "nom.short")], by = c("variable" = "nom.saisi"))  %>% # left join the correct name
  mutate(variable = nom.short) 
cbnmed.melt = cbnmed.melt[,!colnames(cbnmed.melt) == "nom.short"]

# Summerize by combe * placette * quad * variable *year if a species has been recorded under different name between several quadrats in a year
cbnmed.melt.pa = cbnmed.melt %>%
  group_by(combe, placette, quad, variable, year) %>% # define the combination of factors 
  dplyr::summarize(value.pa = sum(value)) %>% # compute the sum due to several observations/year
  mutate(value.pa = replace(value.pa, value.pa >0, 1)) %>% # transform to presence absence
  mutate (placette = as.numeric(placette), 
          year = as.numeric(year)) # ensure that placette & year are numeric

#####
# 2. UB (university of Barcelona)
#####

# Contains the sites CRE, RAT, ULL

### load data
raw.ub = read.csv (here::here("data", "raw", "ub", "UB_data.csv"), 
               head = T, row.names = 1, sep = ";", dec = ",", 
               check.names = FALSE) # carefull : conserve the full species name to concatenate with the referential (column name can be invalid)

raw.ub = as.data.frame(t(raw.ub)) # change direction of the table

# Create the data on site/obs by splitting the name of the line
ub <- raw.ub %>% 
  rownames_to_column(var = "row_name") %>%  # create a variable based on row.names
  separate(row_name,sep = "_",into = c("combe","placette", "quad", "year")) %>%  # split that variable
  relocate (combe, quad, placette, year) # and reorder

### Melt the dataset and homogeneize species name
# melt the dataset
varnames = colnames(ub)
ub.melt = melt(ub, id=varnames[1:4]) # melt the dataset : the value (pheno & abundance) is in the "value" column

# change the phenology to binary data IF NEEDED
# ub.melt = ub.melt %>%
#   mutate(value = replace(value, value != "",1)) %>% # transform the column "values" to build a 1-0 metric 
#   mutate(value = as.numeric (value)) %>% # transform the column to numeric
#   replace_na(list(value = 0)) # transform NA to 0

# change species name for the correct name (nom.short) in the referential
ref.taxo.ub = ref.taxo[which(ref.taxo$source=="ub"),] # create a sub-referential containing only the names linked to the dataset

ub.melt = ub.melt %>%
  left_join(ref.taxo.ub [,c("nom.saisi", "nom.short")], by = c("variable" = "nom.saisi"))  %>% # left join the correct name
  mutate(variable = nom.short) 
ub.melt = ub.melt[,!colnames(ub.melt) == "nom.short"]

# Summerize by combe * placette * quad * variable *year if a species has been recorded under different name between several quadrats in a year
ub.melt.pa = ub.melt %>%
  group_by(combe, placette, quad, variable, year) %>% # define the combination of factors 
  dplyr::summarize(value.pa = sum(value)) %>% # compute the sum due to several observations/year
  mutate(value.pa = replace(value.pa, value.pa >0, 1)) %>% # transform to presence absence
  mutate (placette = as.numeric(placette), 
          year = as.numeric(year)) # ensure that placette & year are numeric

#####
# 3. portillon
#####

# Contains the sites por

### load data
raw.por = read.csv (here::here("data", "raw", "por", "DATAPOR.csv"), 
                   head = T, row.names = 1, sep = ";", dec = ",", 
                   check.names = FALSE) # carefull : conserve the full species name to conporenate with the referential (column name can be invalid)

raw.por = as.data.frame(t(raw.por)) # change direction of the table

# Create the data on site/obs by splitting the name of the line
por <- raw.por %>% 
  rownames_to_column(var = "row_name") %>%  # create a variable based on row.names
  separate(row_name,sep = "_",into = c("combe","placette", "quad", "year")) %>%  # split that variable
  relocate (combe, quad, placette, year) # and reorder

### Melt the dataset and homogeneize species name
# melt the dataset
varnames = colnames(por)
por.melt = melt(por, id=varnames[1:4]) # melt the dataset : the value (pheno & abundance) is in the "value" column

# change the phenology to binary data IF NEEDED
por.melt = por.melt %>%
  mutate(value = replace(value, value != "",1)) %>% # transform the column "values" to build a 1-0 metric
  mutate(value = as.numeric (value)) %>% # transform the column to numeric
  replace_na(list(value = 0)) # transform NA to 0

# change species name for the correct name (nom.short) in the referential
ref.taxo.por = ref.taxo[which(ref.taxo$source=="por"),] # create a spor-referential containing only the names linked to the dataset

por.melt = por.melt %>%
  left_join(ref.taxo.por [,c("nom.saisi", "nom.short")], by = c("variable" = "nom.saisi"))  %>% # left join the correct name
  mutate(variable = nom.short) 
por.melt = por.melt[,!colnames(por.melt) == "nom.short"]

# Summerize by combe * placette * quad * variable *year if a species has been recorded under different name between several quadrats in a year
por.melt.pa = por.melt %>%
  group_by(combe, placette, quad, variable, year) %>% # define the combination of factors 
  dplyr::summarize(value.pa = sum(value)) %>% # compute the sum due to several observations/year
  mutate(value.pa = replace(value.pa, value.pa >0, 1)) %>% # transform to presence absence
  mutate (placette = as.numeric(placette), 
          year = as.numeric(year)) # ensure that placette & year are numeric

#####
# 4. Cat (Cataperdis)
#####

# Contains the sites cat

### load data
raw.cat = read.csv (here::here("data", "raw", "cat", "CAT2013-2019.csv"), 
                    head = T, row.names = 1, sep = ";", dec = ",", 
                    check.names = FALSE) # carefull : conserve the full species name to concatenate with the referential (column name can be invalid)

raw.cat = as.data.frame(t(raw.cat)) # change direction of the table

# Create the data on site/obs by splitting the name of the line
cat <- raw.cat %>% 
  rownames_to_column(var = "row_name") %>%  # create a variable based on row.names
  separate(row_name,sep = "_",into = c("combe","placette", "quad", "year")) %>%  # split that variable
  relocate (combe, quad, placette, year) # and reorder

### Melt the dataset and homogeneize species name
# melt the dataset
varnames = colnames(cat)
cat.melt = melt(cat, id=varnames[1:4]) # melt the dataset : the value (pheno & abundance) is in the "value" column

# change the phenology to binary data IF NEEDED
cat.melt = cat.melt %>%
  mutate(value = replace(value, value != "",1)) %>% # transform the column "values" to build a 1-0 metric
  mutate(value = as.numeric (value)) %>% # transform the column to numeric
  replace_na(list(value = 0)) # transform NA to 0

# change species name for the correct name (nom.short) in the referential
ref.taxo.cat = ref.taxo[which(ref.taxo$source=="cat"),] # create a scat-referential containing only the names linked to the dataset

cat.melt = cat.melt %>%
  left_join(ref.taxo.cat [,c("nom.saisi", "nom.short")], by = c("variable" = "nom.saisi"))  %>% # left join the correct name
  mutate(variable = nom.short) 
cat.melt = cat.melt[,!colnames(cat.melt) == "nom.short"]

# Summerize by combe * placette * quad * variable *year if a species has been recorded under different name between several quadrats in a year
cat.melt.pa = cat.melt %>%
  group_by(combe, placette, quad, variable, year) %>% # define the combination of factors 
  dplyr::summarize(value.pa = sum(value)) %>% # compute the sum due to several observations/year
  mutate(value.pa = replace(value.pa, value.pa >0, 1)) %>% # transform to presence absence
  mutate (placette = as.numeric(placette), 
          year = as.numeric(year)) # ensure that placette & year are numeric

#####
# 5. Arb (Arbella)
#####

# Contains the sites arb

### load data
raw.arb = read.csv (here::here("data", "raw", "arb", "ARB2013-2018.csv"), 
                    head = T, row.names = 1, sep = ";", dec = ",", 
                    check.names = FALSE) # carefull : conserve the full species name to conarbenate with the referential (column name can be invalid)

raw.arb = as.data.frame(t(raw.arb)) # change direction of the table

# Create the data on site/obs by splitting the name of the line
arb <- raw.arb %>% 
  rownames_to_column(var = "row_name") %>%  # create a variable based on row.names
  separate(row_name,sep = "_",into = c("combe","placette", "quad", "year")) %>%  # split that variable
  relocate (combe, quad, placette, year) # and reorder

### Melt the dataset and homogeneize species name
# melt the dataset
varnames = colnames(arb)
arb.melt = melt(arb, id=varnames[1:4]) # melt the dataset : the value (pheno & abundance) is in the "value" column

# change the phenology to binary data IF NEEDED
arb.melt = arb.melt %>%
  mutate(value = replace(value, value != "",1)) %>% # transform the column "values" to build a 1-0 metric
  mutate(value = as.numeric (value)) %>% # transform the column to numeric
  replace_na(list(value = 0)) # transform NA to 0

# change species name for the correct name (nom.short) in the referential
ref.taxo.arb = ref.taxo[which(ref.taxo$source=="arb"),] # create a sarb-referential containing only the names linked to the dataset

arb.melt = arb.melt %>%
  left_join(ref.taxo.arb [,c("nom.saisi", "nom.short")], by = c("variable" = "nom.saisi"))  %>% # left join the correct name
  mutate(variable = nom.short) 
arb.melt = arb.melt[,!colnames(arb.melt) == "nom.short"]

# Summerize by combe * placette * quad * variable *year if a species has been recorded under different name between several quadrats in a year
arb.melt.pa = arb.melt %>%
  group_by(combe, placette, quad, variable, year) %>% # define the combination of factors 
  dplyr::summarize(value.pa = sum(value)) %>% # compute the sum due to several observations/year
  mutate(value.pa = replace(value.pa, value.pa >0, 1)) %>% # transform to presence absence
  mutate (placette = as.numeric(placette), 
          year = as.numeric(year)) # ensure that placette & year are numeric

#####
# 6. Assembling the dataset
#####

### Presentation
head (cbnmed.melt.pa)
head (ub.melt.pa)
head (por.melt.pa)
head (cat.melt.pa)
head (arb.melt.pa)

### Joining the tables
sb_data = cbnmed.melt.pa %>%
  bind_rows(ub.melt.pa, por.melt.pa, cat.melt.pa, arb.melt.pa)

### Cleaning the dataset
sb_data = sb_data %>%
  mutate (combe = str_to_lower(combe), 
          quad = str_to_lower(quad),
          variable = str_to_lower(variable)) # simplify the writing to lower case

sb_data = sb_data %>%
  filter(variable != "delete") # delete all the taxa we excluded from the analysis with "delete"

missing_sp = sb_data %>%
  group_by(variable) %>%
  dplyr::summarize (n_obs = sum(value.pa)) %>%
  filter(n_obs == 0) # give the list of species for which there is no observation
missing_sp = as.vector(missing_sp$variable)
`%notin%` <- Negate(`%in%`) # creating a not in operator
sb_data = sb_data %>%
  filter(variable %notin% missing_sp) # delete all the taxa we showed no data

### Compute the frequency for each species (= variable) on each site (= combe) per plot (= placette) for each year
sb_data_freq = sb_data %>%
  group_by(combe, placette, year, variable) %>%
  dplyr::summarize (n_obs = sum(value.pa)/12)

### Cast the dataset into a large format with species as variables
sb_data_cast = dcast(sb_data_freq, ... ~ variable)# %>% # structure the dataset
sb_data_cast[is.na(sb_data_cast)] = 0 # replace NA with 0

### Delete empty rows = quadrat that were not followed due to snow, etc.
sb_data_cast = sb_data_cast %>%
  mutate (obs_quad = rowSums(sb_data_cast[, -(1:4)])) %>%
  filter(obs_quad != 0) %>%
  dplyr::select (-obs_quad)

as.data.frame(table(sb_data_cast[,1:3])) %>%
  arrange(combe, placette)

### save the dataset
write_csv( sb_data_cast, 
           here::here ("data", "processed", "sb_data_cast.csv"))



## SP name revision
write_csv( full.sp.list %>% mutate(sp = str_replace(sp, "_", " ")), 
           here::here ("data", "processed", "species.list.csv"))

colnames(full.sp.list)  = "sp"
t = full.sp.list %>%
  mutate(sp = str_replace(sp, "_", " "))
TPL(t)



