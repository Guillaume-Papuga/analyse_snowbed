#######################################################
# Project : Analysis of snowbed vegetation dynamics
# Script : 04.data_import_snow
# Load data from temperature loggers n compute indexes
# Authors : Guillaume Papuga
# Last update : 1st July 2022
#######################################################

#####
# 0. Loading data
#####
# Raw data
arb = read.csv (here::here("data", "raw", "temp", "arb.csv"), 
                    head = T, row.names = NULL, sep = ";", dec = ",")
cat = read.csv (here::here("data", "raw", "temp", "cat.csv"), 
                    head = T, row.names = NULL, sep = ";", dec = ",")
cre = read.csv (here::here("data", "raw", "temp", "cre.csv"), 
                    head = T, row.names = NULL, sep = ";", dec = ",")
pdlc = read.csv (here::here("data", "raw", "temp", "pdlc.csv"), 
                     head = T, row.names = NULL, sep = ";", dec = ",")
pla = read.csv (here::here("data", "raw", "temp", "pla.csv"), 
                    head = T, row.names = NULL, sep = ";", dec = ",")
rat = read.csv (here::here("data", "raw", "temp", "rat.csv"), 
                    head = T, row.names = NULL, sep = ";", dec = ",")
ull = read.csv (here::here("data", "raw", "temp", "ull.csv"), 
                    head = T, row.names = NULL, sep = ";", dec = ",")

# Convert time and date
arb$date = as.Date(arb$date, format =  "%d/%m/%Y")
cat$date = as.Date(cat$date, format =  "%d/%m/%Y")
cre$date = as.Date(cre$date, format =  "%d/%m/%Y")
pdlc$date = as.Date(pdlc$date, format =  "%d/%m/%Y")
pla$date = as.Date(pla$date, format =  "%d/%m/%Y")
rat$date = as.Date(rat$date, format =  "%d/%m/%Y")
ull$date = as.Date(ull$date, format =  "%d/%m/%Y")

# Create a list of DF
data.lst = lst(arb, cat, cre, pdlc, pla, rat, ull)

#####
# 1. Combler les trous
#####


#####
# 2. Find the correct melting/freezing date for each logger
#####
### A. Prepare a synthesis table
synth = data.frame ()

# Loop on the site
site.name = c("arb", "cat", "cre", "pdlc", "pla", "rat", "ull")

for (p in site.name) {
  # Daily dataset
  site.data = data.lst[[p]]
  data.sum = site.data %>%
    group_by(plot, year, date) %>% # group by different factors
    dplyr::summarize(max.t = max(temperature), # compute indexes
                     min.t = min (temperature),
                     var = var(temperature)) %>%
    mutate (mean.t = (max.t + min.t)/2) # compute the mean temperature for GDD
  
  time.lag = 7 # duration (in days) of continuous snow cover
  data.sum$max.var = frollapply(data.sum$var, time.lag, max) # cumulative value of variance to avoid "rare snow event"

  # Save each dataset for GDD
  assign(paste ("daily.data_", p, sep = ""), data.sum) 
  
  # Boucle for sur une année
  for (i in 2012:2020) {
    # Define the middle of the growing season
    tab.year = data.sum %>%
      dplyr::filter(year == i) # boucle sur l'année
    
    for (j in c("Early", "Intermediate", "Late")) {
      tab = tab.year %>%
        dplyr::filter (plot == j) %>%
        dplyr::filter (var == max(var))
      
      midyear = tab[1, "date"] # approximate peak of the summer season, cut the year in two
      
      # Find the melting day
      m.day = data.sum %>%
        dplyr::filter(date >= paste( i, "-01-01", sep = "") & date <= midyear) %>%
        dplyr::filter(plot == j) %>%
        dplyr::filter (max.var < 1) %>% # delete values > 1 (melting threshold)
        dplyr::filter (date == max(date)) %>%
        dplyr::mutate (site = p, 
                       time = "m.day")
      
      # Find the freezing day
      f.day = data.sum %>%
        dplyr::filter(date >= midyear & date <= paste( i, "-12-31", sep = "")) %>%
        dplyr::filter(plot == j) %>%
        dplyr::filter (max.var < 1) %>% # delete values > 1 (melting threshold)
        dplyr::filter (date == min(date)) %>%
        dplyr::mutate (site = p, 
                       time = "f.day", 
                       date = date - time.lag) 
      
      # Assigne the correct name
      assign (paste ("m.day_", j, sep = ""), m.day) # Assign the correct name to each result
      assign (paste ("f.day_", j, sep = ""), f.day)
    }
    # Merge
    synth = bind_rows(synth, m.day_Early, f.day_Early, 
                      m.day_Intermediate, f.day_Intermediate, 
                      m.day_Late, f.day_Late)
  }
}

# Rearrange the table
synth = synth %>%
  relocate (site, time, plot, year) %>%
  arrange (site, desc(time), plot, year)
synth

#####
# 3. Calculer les indices
#####
# Table de synthèse
# A. Basic index
tab.data.synth = dcast (synth %>% dplyr::select (-max.t, -min.t, -var, -mean.t, -max.var), # delete useless variables
                        site + year + plot ~ time) %>%
  mutate (m.day = as.Date(m.day), 
          f.day = as.Date(f.day)) %>% # transform as Date
  relocate (site, year, plot, m.day, f.day) %>% # reorder the columns
  mutate (first.day = as.Date(ISOdate(year, 1, 1)), 
          n_m.day = as.numeric(difftime (m.day, first.day)), 
          n_f.day = as.numeric(difftime (f.day, first.day)), 
          n_vege = as.numeric(difftime (f.day, m.day)))

# B. Snow length during the previous winter
for (i in 1:nrow(tab.data.synth)){
  # Characteristics
  st = tab.data.synth[i, "site"]
  yr = tab.data.synth[i, "year"]
  plt = tab.data.synth[i, "plot"]
  
  # Compute data
  val = tab.data.synth %>% 
    filter (site == st, 
            year == yr - 1, 
            plot == plt) %>%
    dplyr::select (f.day)
  
  tab.data.synth[i, "prev_f.day"] = (val$f.day)
}

tab.data.synth = tab.data.synth %>%
  mutate (past.w_snow_day = as.numeric(difftime (m.day, prev_f.day)))
          
# C. Growing degree days GDD & highest mean-t (95%) & absolute highest temperature
daily.data = lst (daily.data_arb,  # create a list of dataset to import them automatically in the loop
                  daily.data_cat, 
                  daily.data_cre, 
                  daily.data_pdlc, 
                  daily.data_pla, 
                  daily.data_rat, 
                  daily.data_ull)

for (i in 1:nrow(tab.data.synth)){
  # Characteristics
  st = tab.data.synth[i, "site"]
  yr = tab.data.synth[i, "year"]
  plt = tab.data.synth[i, "plot"]
  m.day.yr = tab.data.synth[i, "m.day"]
  f.day.yr = tab.data.synth[i, "f.day"]
  
  # Select the corresponding daily data
  d.data = daily.data[[paste("daily.data_", st, sep = "")]] # based on the name of the site 'st'
  
  # Compute the GDD index
  gdd = d.data %>%
    filter (plot == plt, # filter the correct plot
            date >= m.day.yr & date <= f.day.yr) %>% # select the correct vegetative period
    dplyr::summarise (gdd = sum (mean.t))
  
  # Compute the 95% mean temperature
  ht.95 = d.data %>%
    filter (plot == plt, # filter the correct plot
            date >= m.day.yr & date <= f.day.yr) %>% # select the correct vegetative period
    dplyr::summarise (ht.95 = as.numeric(quantile (mean.t, 0.95)))
  
  # Compute the absolute max temperature
  max.t = d.data %>%
    filter (plot == plt, # filter the correct plot
            date >= m.day.yr & date <= f.day.yr) %>% # select the correct vegetative period
    dplyr::summarise (max.t = max (max.t))
  
  # Type the data into the dataset
  if (length (gdd$gdd) != 0) {tab.data.synth[i, "gdd"] = gdd$gdd}
  if (length (gdd$gdd) != 0) {tab.data.synth[i, "ht.95"] = ht.95$ht.95}
  if (length (gdd$gdd) != 0) {tab.data.synth[i, "max.t"] = max.t$max.t}
}

head(tab.data.synth)

#####
# 4. Production graph et tab
#####
### A. Synthetic table
write_csv( tab.data.synth, 
           here::here ("data", "processed", "tab.microclimate.detailed.csv"))

### B. Plot temp ~ time
site.name
for (i in site.name) {
  # Dataset
  dtst = data.frame(data.lst[[i]]) 
  
  onset = tab.data.synth %>% filter (site == i) %>% dplyr::select (plot, m.day)
  end = tab.data.synth %>% filter (site == i) %>% dplyr::select (plot, f.day)
  
  # Plot
  plt = ggplot(data = dtst, aes(x = date, y = temperature, group = plot))+
    geom_line() + 
    facet_grid(plot~.) +
    theme_bw() +
    ggtitle(paste (i)) +
    geom_vline(data = onset, aes(xintercept = onset$m.day, group = plot), color="red", linetype="dashed", size=0.5) + 
    geom_vline(data = end, aes(xintercept = end$f.day, group = plot), color="blue", linetype="dashed", size=0.5)
  
  # Save
  ggsave(paste("snow.plot_", i, ".pdf", sep = ""),
         path = here::here("figures/"),
         units = "px",
         width = 600, height = 350,
         dpi = 100,
         limitsize = TRUE)
  }

### C. Plot variable ~ year
# Table with mean values per year (after sorting)
synth.yr = tab.data.synth %>%
  dplyr::filter (!(site == "pdlc" & plot == "Early"), # exclude Early pdlc due to ibutton issue
                 !(site == "cre" & plot == "Early")) %>% # same issue (CRE1 already excluded from floristic analysis)
  dplyr::group_by(site, year) %>%
    dplyr::summarize(n_vege.yr = mean(n_vege, na.rm = TRUE), 
                     n_m.day.yr = mean(n_m.day, na.rm = TRUE), 
                     n_f.day.yr = mean(n_f.day, na.rm = TRUE), 
                     past.w_snow_day.yr = mean(past.w_snow_day, na.rm = TRUE), 
                     gdd.yr = mean(gdd, na.rm = TRUE), 
                     ht.95.yr = mean(ht.95, na.rm = TRUE))

# Plot
var = colnames (synth.yr)[-(1:2)]
for (i in var) {
  # Dataset
  dtst = synth.yr %>% dplyr::select (, site, year, i)
  colnames (dtst) [3] = "variable"

  # Plot
  plt = ggplot(data = dtst, aes(x = year, y = variable, colour = site))+
    geom_line() + 
    theme_bw() +
    ggtitle(paste (i)) 

  plt
  
  # Save
  ggsave(paste("snow.trend_", i, ".pdf", sep = ""),
         path = here::here("figures/microclim_trend/"),
         units = "px",
         width = 600, height = 350,
         dpi = 100,
         limitsize = TRUE)
  }

### D. Data for the RDA analysis
# Dataset
synth.rda = tab.data.synth %>%
  dplyr::group_by(site, plot) %>%
  dplyr::summarise_all(.funs = mean, na.rm = TRUE) %>% # Average over the study period
  dplyr::filter (!(site == "cre" & plot == "Early")) %>% # Exclude CRE 1
  dplyr::select (site, plot, n_m.day, n_f.day, n_vege, past.w_snow_day, gdd, ht.95) %>%
  dplyr::mutate (plot = as.character(plot)) %>%
  dplyr::mutate(plot = replace(plot, plot == "Early", "1"), 
                plot = replace(plot, plot == "Intermediate", "2"), 
                plot = replace(plot, plot == "Late", "3")) %>%
  dplyr::mutate(plot = as.integer(plot))
  
# Autocorrelation in the variables
corrplot(cor(synth.rda[, -(1:2)]), 
         method = 'number')

# Save the dataset
write.table(synth.rda, 
            here::here("data", "processed", "synth.microclim.rda.csv"))


########################################################################################################