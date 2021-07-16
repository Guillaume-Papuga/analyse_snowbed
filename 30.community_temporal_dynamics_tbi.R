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

### Create the 2 matrices
start_year = 2013
end_year = 2019

# For start
m_start = df %>% 
  filter (year == start_year) %>% # select the appropriate year
  unite(name_quad, c(combe, placette), sep = "_", remove = FALSE) %>% # create the rowname..
  select (-combe, -placette, -year) # delete other variables

# For end
m_end = df %>% 
  filter (year == end_year) %>% # select the appropriate year
  unite(name_quad, c(combe, placette), sep = "_", remove = FALSE) %>% # create the rowname..
  select (-combe, -placette, -year) # delete other variables

# Find the quadrat that have been sensused the two years
quad.list = as.vector (unique (inner_join(m_start, m_end, by = "name_quad")$name_quad))

# Flter the two datasets
# start
m_start = m_start %>%
  filter (name_quad %in% quad.list) %>%
  column_to_rownames(var = "name_quad")  # and add it to row

# end
m_end = m_end %>%
  filter (name_quad %in% quad.list) %>%
  column_to_rownames(var = "name_quad") # and add it to row

### TBI analysis
# compute the index
tbi_st.end <- TBI(m_start,
                    m_end, 
                    method = "%difference",
                    pa.tr = F, 
                    nperm = 999,
                    BCD = TRUE, 
                    save.BC = T,
                    test.BC = TRUE, 
                    test.t.perm = T, 
                    seed. =   )

# add columns to BCD.mat
row.names(tbi_st.end$BCD.mat) = row.names(m_start) # change the label for the B-C plot
col1 = row.names(m_start) # create 2 variables based on plot name (combe & quad)
col1 = as.data.frame(col1)
col1 = separate(col1, col=1, into = c("combe","placette"), sep = "_") # combe and plot group labels
tbi_st.end$BCD.mat = cbind.data.frame(tbi_st.end$BCD.mat, col1)
tbi_st.end$BCD.mat$combe = as.factor(tbi_st.end$BCD.mat$combe)

# build the signifance table for each site
sp_site = rbind.data.frame(tbi_st.end$TBI,
                           tbi_st.end$p.TBI, 
                           tbi_st.end$p.adj)

colnames(sp_site) = rownames(m_start)
sp_site = t(sp_site)

# significance table for global analysis
bc_sign = cbind(tbi_st.end$BCD.summary, # 
                tbi_st.end$t.test_B.C)


### Plot
# Define the panel
theme<-theme(panel.background = element_blank(),
             panel.border=element_rect(fill=NA),
             panel.grid.major = element_blank(),
             panel.grid.minor = element_blank(),
             strip.background=element_blank(),
             axis.text.x=element_text(colour="black"),
             axis.text.y=element_text(colour="black"),
             axis.ticks=element_line(colour="black"),
             plot.margin=unit(c(1,1,1,1),"line"))

# Plot
tbi.plot = plot(tbi_st.end, s.names=rownames(tbi_st.end$BCD.mat)  ,  #plot B-C
                col.rim = "coral", 
                pch.loss=19 , 
                pch.gain=15, 
                main="Temporal change in beta diversity", 
                xlim=c(0,0.4), ylim=c(0,0.4))


# Save the plot

