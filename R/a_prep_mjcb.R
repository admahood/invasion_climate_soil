## libraries -------------------------------------------------------------------
library(tidyverse)

# functions ====================================================================
pct_to_kg_ha <- function(pct, bulk_density, depth){
  #1 % = 1000 mg.kg-1 = 10 g.kg-1= 0.001 kg.kg-1
  ha= 10^8 #cm2
  soil_kg_per_ha = depth * ha * bulk_density / 1000 
  return(soil_kg_per_ha * pct/100)
}

# look up tables =========================================================
lut_variables <- c("Bromus_TN_pct" = "Bromus N (%)", "Bromus_TC_pct" = "Bromus C (%)",
                   "Bromus_CN" = "Bromus C:N", "Other_TN_pct" = "Other N (%)",
                   "Other_TC_pct" = "Other C (%)", "Other_CN" = "Other C:N", 
                   "Poa_TN_pct" = "Poa N (%)", "Poa_TC_pct" = "Poa C (%)",
                   "Poa_CN" = "Poa C:N", "Litter_TN_pct" = "Litter N (%)",
                   "Litter_TC_pct" = "Litter C (%)", "Litter_CN" = "Litter C:N",
                   "SOIL_SurSo4_kg_ha" = "Soil SO<sub>4</sub> (kg/ha)",
                   "SOIL_Ca_kg_ha" = "Soil Ca (kg/ha)", 
                   "SOIL_Mg_kg_ha" = "Soil Mg (kg/ha)",
                   "SOIL_CN" = "Soil C:N", "soil_n_kg_ha" = "Soil Total N (kg/ha)",
                   "soil_c_kg_ha" = "Soil Total C (kg/ha)",
                   "total_mineral_n" = "Soil Total Mineral N (kg/ha)",
                   "NO3_kg_ha" = "Soil Nitrate (kg/ha)",
                   "NH4_kg_ha" = "Soil Ammonium (kg/ha)",
                   "ja_ju_def" = "CWD",
                   "ja_ju_aet" = "AET",
                   "tmin" = "T_min",
                   "Annuals" = "Annuals",
                   "Perennials" = "Perennials",
                   "Forbs" = "Forbs",
                   "Grasses" = "Grasses")

lut_variables_nounit <- c("Bromus_TN_pct" = "Bromus N", "Bromus_TC_pct" = "Bromus C",
                          "Bromus_CN" = "Bromus C:N", "Other_TN_pct" = "Other N",
                          "Other_TC_pct" = "Other C", "Other_CN" = "Other C:N", 
                          "Poa_TN_pct" = "Poa N", "Poa_TC_pct" = "Poa C",
                          "Poa_CN" = "Poa C:N", "Litter_TN_pct" = "Litter N",
                          "Litter_TC_pct" = "Litter C", "Litter_CN" = "Litter C:N",
                          "SOIL_SurSo4_kg_ha" = "Soil SurSo4",
                          "SOIL_Ca_kg_ha" = "Soil Ca", 
                          "SOIL_Mg_kg_ha" = "Soil Mg",
                          "SOIL_OM_pct" = "Soil Total C",
                          "SOIL_TN_pct" = "Soil Total N",
                          "SOIL_CN" = "Soil C:N", "soil_n_kg_ha" = "Soil Total N",
                          "soil_c_kg_ha" = "Soil Total C",
                          "total_mineral_n" = "Soil Mineral N",
                          "NO3_kg_ha" = "Soil Nitrate",
                          "NH4_kg_ha" = "Soil Ammonium",
                          "ja_ju_def" = "CWD",
                          "ja_ju_aet" = "AET",
                          "tmin" = "T_min",
                          "Annuals" = "Annuals",
                          "Perennials" = "Perennials",
                          "Forbs" = "Forbs",
                          "Grasses" = "Grasses")

site_colors = c("springgreen4","deepskyblue2",
           "darkgoldenrod2","chocolate4")
## data prep -------------------------------------------------------------------
# adding precip with previously generated climate data from prism
# and year as a factor, then joining the two data frames

ppt <- read.csv("data/ppt_sums.csv") %>%
  dplyr::select(-X, 
                ppt2011_2013 = X2011.2013_ppt,
                ppt2014_2016 = X2014.2016_ppt) 

# adding more climate variables
climate_sums <- read.csv("data/climate_summaries.csv") %>%
  dplyr::select(-X)

# cleaned and matched data frames
glmm13 <- read.csv("data/Jones_2013_transect_dec17.csv") %>%
  mutate(Year = 2013) %>%
  left_join(ppt,by=c("Site_number"="plot")) %>%
  left_join(climate_sums, by=c("Site_number"="plot")) %>%
  dplyr::select(-ppt2014_2016, -ppt14_15,-ppt15_16,-tmax16,-tmin_16,-vpdmax_16,
                ppt_2yrtot = ppt2011_2013,
                ppt_2yr = ppt11_12,
                ppt_1yr = ppt12_13,
                tmax = tmax_13,
                tmin = tmin_13,
                vpdmx = vpdmax_13)

glmm16 <- read.csv("data/jones_all_nov_2017.csv") %>%
  mutate(Year = 2016) %>%
  left_join(ppt, by = c("Site_number"="plot")) %>%
  left_join(climate_sums, by=c("Site_number"="plot")) %>%
  dplyr::select(-ppt2011_2013, -ppt11_12,-ppt12_13,-tmax_13,-tmin_13,-vpdmax_13,
                # -ANF, -PNF,
                ppt_2yrtot = ppt2014_2016,
                ppt_2yr = ppt14_15,
                ppt_1yr = ppt15_16,
                tmax = tmax16,
                tmin = tmin_16,
                vpdmx = vpdmax_16)


# making Shrubs use the same data from 2013 (that was line
# intercept, I just included artr in my quadrats)

glmm16$Shrubs <- filter(glmm13,
                        Site_number != 4, 
                        Site_number != 5,
                        Site_number != 19,
                        Site_number != 20,
                        Site_number != 21)$Shrubs

# # center and scale soil SO4 by year

# relative cover of Crypto by year
# there was a discrepancy in measurement techinque (2016 was counting some
# types of moss as crypto)
glmm16$Crypto <- glmm16$Crypto / max(glmm16$Crypto)
glmm13$Crypto <- glmm13$Crypto / max(glmm13$Crypto)

# grabbing some soil stuff that's inconsistently collected but still data

soil_stuff16 <- read_csv("data/jones_all_oct5.csv") %>%
  dplyr::select(Plot_TP, NO3_ppm = "NO3.N..ppm.", 
                bulk_density = "bulkDensity.g.cm3.", NH4_ppm = "NH4.N..ppm.",
                net_mineralization = net.mineralization,
                soil_moisture = soilMoisture) %>%
  mutate(Year = "2016",
         NO3_kg_ha = NO3_ppm * bulk_density,
         NH4_kg_ha = NH4_ppm * bulk_density,
         Site_number = as.character(as.numeric(substr(Plot_TP, 2,3))),
         Transect = substr(Plot_TP,4,4))

soil_stuff13 <- read_csv("data/data_2013/Soils 10-23-13.csv") %>%
  dplyr::select(Site_number = `Site number`, Transect, NO3_kg_ha = `Surno3 (kg/ha)`, 
                NH4_kg_ha = `NH4 (kg/ha)`) %>%
  mutate(net_mineralization = NA,
         Year = "2013")

bulk_dens<- soil_stuff16%>%
  dplyr::select(Site_number, Transect, bulk_density) %>%
  group_by(Site_number) %>%
  dplyr::summarise(bulk_d = mean(bulk_density)) %>%
  ungroup()%>%
  mutate(Site_number = as.numeric(Site_number))

# couldn't find the bd measurements from '13, using '16, but 5 sites weren't 
# repeated in '16 so just using the mean across sites for those
bd_mean <- mean(bulk_dens$bulk_d)
bulk_dens[21, "Site_number"] <- 19; bulk_dens[21, "bulk_d"] <- bd_mean
bulk_dens[22, "Site_number"] <- 20; bulk_dens[22, "bulk_d"] <- bd_mean
bulk_dens[23, "Site_number"] <- 21; bulk_dens[23, "bulk_d"] <- bd_mean
bulk_dens[24, "Site_number"] <- 4; bulk_dens[24, "bulk_d"] <- bd_mean
bulk_dens[25, "Site_number"] <- 5; bulk_dens[25, "bulk_d"] <- bd_mean

dp = 10

lut_st <- c("M" = "Invaded Sagebrush", "I" = "Intact Sagebrush",
            "C" = "Cheatgrass-dominated", "D" = "Cheatgrass Dieoff")

all <- rbind(glmm13, glmm16) %>%
  left_join(bulk_dens) %>%
  mutate(Year = as.factor(Year),
         shrub_b = as.factor(ifelse(Site.type == "I" | Site.type == "M", "Shrub", "Grass")),
         Site_number = as.factor(Site_number),
         soil_n_kg_ha = pct_to_kg_ha(pct = SOIL_TN_pct, depth = dp, bulk_density = bulk_d),
         soil_c_kg_ha = pct_to_kg_ha(pct = SOIL_OM_pct/1.724, depth = dp, bulk_density = bulk_d),
         Site.type = factor(lut_st[as.character(Site.type)], 
                            levels = c("Intact Sagebrush","Invaded Sagebrush",
                                       "Cheatgrass-dominated", "Cheatgrass Dieoff")))%>%
  as_tibble()



# grabbing some individual species ---------------------------------------------

soil_stuff <- rbind(soil_stuff13, 
                    dplyr::select(soil_stuff16,
                                  Site_number,
                                  Transect,
                                  Year,
                                  NO3_kg_ha,
                                  NH4_kg_ha,
                                  net_mineralization)) %>%
  mutate(total_mineral_n = NO3_kg_ha +NH4_kg_ha)


all <- all %>%
  left_join(dplyr::select(soil_stuff, Site_number,
                          Transect,
                          Year,
                          NO3_kg_ha,
                          NH4_kg_ha,
                          total_mineral_n), by = c("Year", "Site_number", "Transect"))



all_p <- all %>%
  dplyr::select(-Site.type, -shrub_b) %>%
  group_by(Site_number,Year) %>%
  summarise_all(mean, na.rm=T) %>%
  left_join(unique(dplyr::select(all, Site_number, Site.type, shrub_b))) %>%
  mutate(Site.type = as.character(Site.type),
         shrub_b = as.character(shrub_b)) %>%
  dplyr::select(-Transect) %>%
  ungroup()%>%
  left_join(read_csv("data/climate_summaries_def_aet_tmin.csv") %>%
              mutate(Year = as.character(year), Site_number = as.character(plot)) %>%
              dplyr::select(-year,-plot)
            ) %>%
  left_join(read_csv("data/growing_degree_days.csv", col_types = "ffd")) %>%
  ungroup() %>%
  as_tibble()

all_se <- all %>%
  dplyr::select(-Site.type, -shrub_b) %>%
  group_by(Site_number,Year) %>%
  summarise_all(function(x,...) sd(x, na.rm=T)/sqrt(length(x)), na.rm=T) %>%
  left_join(unique(dplyr::select(all, Site_number, Site.type, shrub_b))) %>%
  mutate(Site.type = as.character(Site.type),
         shrub_b = as.character(shrub_b)) %>%
  dplyr::select(-Transect) %>%
  ungroup()

all_sd <- all %>%
  dplyr::select(-Site.type, -shrub_b) %>%
  group_by(Site_number,Year) %>%
  summarise_all(function(x,...) sd(x, na.rm=T)) %>%
  left_join(unique(dplyr::select(all, Site_number, Site.type, shrub_b))) %>%
  mutate(Site.type = as.character(Site.type),
         shrub_b = as.character(shrub_b)) %>%
  dplyr::select(-Transect) %>%
  ungroup()

