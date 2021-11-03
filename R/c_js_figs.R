# setup ==========
source("R/a_prep_mjcb.R")
libs <- c("viridis", "raster", "rcompanion","lsmeans","vegan", "broom", "ggtext",
          "multcomp", "emmeans", "cowplot", "scales", "ggrepel")
iini <-function(x){
  # stands for install if not installed
  if (!x %in% rownames(installed.packages())) install.packages(x)
}
lapply(libs, iini)
lapply(libs, library, character.only = TRUE, verbose = FALSE)




# Figure 1: climographs --------------------------------------------------------

## raster data 
plots <- read.csv("data/locations.csv")
coordinates(plots) <- ~Longitude+Latitude
gm_crs <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
crs(plots) <- gm_crs
pts <- SpatialPoints(coordinates(plots)[,c(1,2)])

if(!file.exists("data/def_aet_tmin.csv")){
  output_data <- data.frame(plot = NA, 
                           date = NA, 
                           variable = NA, 
                           value = NA)
  
  variables <- c("def", "aet")
  counter=1
  for(v in 1:length(variables)){
    aet_files <- list.files(paste0("/home/a/data/climate/",variables[v]), 
                            full.names = TRUE)[33:38]
    for (i in 1:length(aet_files)){
      rr <- stack(aet_files[i])
      print(paste(variables[v], substr(aet_files[i],30,33)))
      ex <- raster::extract(rr, pts) %>% as.data.frame()
      for(m in 1:12){
        for(p in 1:25){
          nn <- ex[p,m]
        
          output_data[counter,1] <- plots@data$Site.number[p]
          output_data[counter,2] <- paste(substr(aet_files[i],30,33), 
                                          if(m<9){paste0("0",as.character(m))}else{m}, 
                                          "01",
                                          sep="-")
          output_data[counter,3] <- variables[v]
          output_data[counter,4] <- nn
          counter<-counter+1
          
        }
      }
    }
  }
  
  output_data$date <- as.Date(output_data$date)
  
  mean_def <- mean(output_data[output_data$variable=="def",]$value)
  sd_def <- sd(output_data[output_data$variable=="def",]$value)
  mean_aet <- mean(output_data[output_data$variable=="aet",]$value)
  sd_aet <- sd(output_data[output_data$variable=="aet",]$value)
  
  norms1<- data.frame(variable = c("def", "aet"),
                      mean = c(mean_def, mean_aet),
                      sd = c(sd_def, sd_aet))
  output_data <- left_join(output_data,norms1)
  output_data <- mutate(output_data, z = (value-mean)/sd)
  
  write.csv(output_data, "data/def_aet.csv")
  
  
  prism<-read.csv("data/monthly_climate.csv") %>%
    select(-X) %>%
    mutate(date = paste0(date, "01"))%>%
    mutate(date = as.Date(date, "%Y%m%d")) %>%
    filter(variable == "tmin")
    
  mean_tmin <- mean(prism[prism$variable=="tmin",]$value)
  sd_tmin <- sd(prism[prism$variable=="tmin",]$value)
  mean_tmax <- mean(prism[prism$variable=="tmax",]$value)
  sd_tmax <- sd(prism[prism$variable=="tmax",]$value)
  
  norms2<- data.frame(variable = c("tmin", "tmax"),
                      mean = c(mean_tmin, mean_tmax),
                      sd = c(sd_tmin, sd_tmax))
  prism <- left_join(prism,norms2)
  prism <- mutate(prism, z = (value-mean)/sd)
  
  output_data <- rbind(output_data,prism)
  
  write.csv(output_data, "data/def_aet_tmin.csv")
}

def_aet_norms <- read_csv("data/def_aet_norms.csv") %>%
  mutate(plot = as.character(plot)) %>%
  dplyr::select(-`...1`) %>%
  dplyr::rename(norm=value)
  
output_data <- read_csv("data/def_aet_tmin.csv")%>%
  dplyr::select(-`...1`)%>%
  mutate(month = as.numeric(substr(date, 6,7)), 
         plot = as.character(plot)) %>%
  left_join(def_aet_norms, by=c("plot","variable","month")) %>%
  mutate(anomaly = value - norm)

mean_tmin <- mean(output_data[output_data$variable=="tmin",]$anomaly)
sd_tmin <- sd(output_data[output_data$variable=="tmin",]$anomaly)
mean_def <- mean(output_data[output_data$variable=="def",]$anomaly)
sd_def <- sd(output_data[output_data$variable=="def",]$anomaly)
mean_aet <- mean(output_data[output_data$variable=="aet",]$anomaly)
sd_aet <- sd(output_data[output_data$variable=="aet",]$anomaly)

norms3<- data.frame(variable = c("tmin", "def", "aet"),
                    amean = c(mean_tmin, mean_def, mean_aet),
                    asd = c(sd_tmin, sd_def, sd_aet))

output_data <- left_join(output_data, norms3, by = "variable") %>%
  mutate(z_anom = (anomaly - amean)/asd)

# sidetrack, writing out a climate summary for sems
# 2013
sample13 <- as.Date("2013-06-01")
sample16 <- as.Date("2013-06-01")

rbind(
output_data %>%
  filter(date < sample13 & date > sample13-180) %>%
  group_by(plot,variable) %>%
  summarise(mx = max(z_anom),
            mn = min(z_anom),
            mean = mean(z_anom))%>%
  ungroup() %>%
  mutate(Year=2013)
,
output_data %>%
  filter(date < sample16 & date > sample16-180) %>%
  group_by(plot,variable) %>%
  summarise(mx = max(z_anom),
            mn = min(z_anom),
            mean = mean(z_anom)) %>%
  ungroup()%>%
  mutate(Year=2016)
)%>%pivot_wider(names_from ="variable", values_from = c(mx, mn, mean))  %>%
  dplyr::rename(Site_number=plot) %>%
  write.csv("data/climate_anomalies_summary.csv")

## weather station data 
norms <- read_csv("data/winnemucca_avgs.csv")

weather <- read.csv("data/winnemucca_airport_weather_10_17.csv", stringsAsFactors = FALSE)
weather <- weather[,c(1:7,9,12,14,17:19,21,24:25)]
weather$date = as.Date(weather$DATE, "%m/%d/%y %H")
weather$year_month = substr(as.character(weather$date), 1,7)
for(i in 8:16){
  weather[,i] = as.numeric(weather[,i])
}
weather$SVP = 610.7* (10^((7.5*weather$HOURLYDRYBULBTEMPC)/(237.3 + weather$HOURLYDRYBULBTEMPC)))
weather$VPD = ((100 - weather$HOURLYRelativeHumidity)/100) * weather$SVP

library(doBy)
monthly_precip_sums = summaryBy(HOURLYPrecip ~ year_month, data = na.omit(weather), FUN = sum)
monthly_precip_sums$date = paste0(monthly_precip_sums$year_month, "-01")
monthly_precip_sums$date = as.Date(monthly_precip_sums$date)
monthly_precip_sums$month = substr(monthly_precip_sums$year_month,6,7)
monthly_precip_sums$year = substr(monthly_precip_sums$year_month,1,4)

yearly_precip_sums = summaryBy(HOURLYPrecip.sum ~ year, data = monthly_precip_sums, FUN = sum)
yearly_precip_sums$Precip_mm_year <- yearly_precip_sums$HOURLYPrecip.sum.sum *25.4
yearly_precip_sums$HOURLYPrecip.sum.sum = NULL
yearly_precip_sums <-yearly_precip_sums[c(2:7),] #gettting rid of incomplete years
yearly_precip_sums$year = as.numeric(yearly_precip_sums$year)

monthly_precip_sums$Precip_mm <- monthly_precip_sums$HOURLYPrecip.sum *25.4
monthly_precip_sums$year = as.numeric(substr(monthly_precip_sums$year_month,1,4))

monthly_precip_sums <- left_join(monthly_precip_sums, yearly_precip_sums) 
monthly_precip_sums <- monthly_precip_sums[5:76,]
monthly_precip_sums$avg_precip <- as.numeric(norms[3,2:13])*25.4 %>% rep(6)
monthly_precip_sums$sd_precip <- sd(monthly_precip_sums$Precip_mm)
monthly_precip_sums <- mutate(monthly_precip_sums, z = (Precip_mm - avg_precip)/sd_precip) 



plines <- data.frame(means = (monthly_precip_sums$Precip_mm_year-210)/
                              monthly_precip_sums$sd_precip/12,
                            years = monthly_precip_sums$year,
                            month = monthly_precip_sums$month,
                            date = monthly_precip_sums$date) %>%
  filter(month == "01" | month == "12")

starts <- plines %>% filter(month =="01") %>% dplyr::rename(x = date, y = means)
ends <- plines %>% filter(month == "12") %>% dplyr::rename(xend = date, yend=means)
segs <- left_join(starts, ends, by = "years")


## plotting climographs 
pp1 <- ggplot(data=monthly_precip_sums, aes(x=date, y=z)) +
  geom_vline(xintercept = as.numeric(as.Date("2013-04-01")), linetype=5) +
  geom_vline(xintercept = as.numeric(as.Date("2016-03-25")), linetype=5) +
  geom_vline(xintercept = as.numeric(as.Date("2013-06-01")), linetype=3) +
  geom_vline(xintercept = as.numeric(as.Date("2016-06-20")), linetype=3) +
  geom_bar(stat="identity") +
  geom_segment(data = segs, aes(x=x,y=y,yend=yend, xend=xend,group=years),
               col = "red", lwd = 1.5)+
  ylab(NULL) +
  xlab(NULL) +
  annotate("text",x= as.Date("2011-09-01"), y=1.5, 
           label = "Annual Means",
           hjust = "left") +
  geom_segment(x = as.Date("2011-06-01"),
               xend = as.Date("2011-08-01"), 
               y = 1.5, yend=1.5,
               col = "red", lwd = 1.5)+
  theme_pubr() +
  scale_x_date(date_labels = "%b %Y",date_breaks = "1 year")+
  theme(legend.position = c(0,1),
        legend.justification = c(0,1),
        legend.title=element_blank(),
        legend.background = element_rect(fill = "transparent"));pp1
  
oc1<-ggplot(data=monthly_precip_sums, aes(x=date, y=anomaly)) +
  geom_vline(xintercept = as.numeric(as.Date("2013-04-01")), linetype=5) +
  geom_vline(xintercept = as.numeric(as.Date("2016-03-25")), linetype=5) +
  geom_vline(xintercept = as.numeric(as.Date("2013-06-01")), linetype=3) +
  geom_vline(xintercept = as.numeric(as.Date("2016-06-20")), linetype=3) +
  geom_hline(yintercept = 0)+
  geom_line(data = output_data %>% filter(variable == "def"), 
            aes(y=z_anom, x=date, color = variable, group = paste(variable, plot)),
              alpha =0.2)+
  ylab(NULL) +
  xlab(NULL) +
  #annotate("text",x= 2013, y=213, label = "30 year normal = 210 mm") +
  theme_pubr() +
  # scale_y_continuous(breaks = c(-1,0,2,4))+
  scale_x_date(date_labels = "%Y") +
  scale_color_manual(values = "red", labels="Climatic Water Deficit") +
  theme(axis.text.x = element_blank()) +
  theme(legend.position = "none",
        legend.justification = c(0,1),
        legend.title=element_blank(),
        legend.background = element_rect(fill = "transparent"))+
  guides(colour = guide_legend(override.aes = list(alpha = 1)))

oc2<-ggplot(data=monthly_precip_sums, aes(x=date, y=anomaly)) +
  geom_vline(xintercept = as.numeric(as.Date("2013-04-01")), linetype=5) +
  geom_vline(xintercept = as.numeric(as.Date("2016-03-25")), linetype=5) +
  geom_vline(xintercept = as.numeric(as.Date("2013-06-01")), linetype=3) +
  geom_vline(xintercept = as.numeric(as.Date("2016-06-20")), linetype=3) +
  geom_hline(yintercept = 0)+
  geom_line(data = output_data%>% filter(variable == "aet"), 
            aes(y=z_anom,
                                    x=date, color = variable, group = paste(variable, plot)),
            alpha =0.2)+
  ylab(NULL) +
  xlab(NULL) +
  #annotate("text",x= 2013, y=213, label = "30 year normal = 210 mm") +
  theme_pubr() +
  # scale_y_continuous(breaks = c(-1,0,2,4))+
  scale_x_date(date_labels = "%Y") +
  scale_y_continuous(breaks=c(-2,0,2,4))+
  scale_color_manual(values = "blue", labels= "Actual Evapotranspiration") +
  theme(axis.text.x = element_blank()) +
  theme(legend.position = "none",
        legend.justification = c(0,1),
        legend.title=element_blank(),
        legend.background = element_rect(fill = "transparent"))+
  guides(colour = guide_legend(override.aes = list(alpha = 1)))

oc3<-ggplot(data=monthly_precip_sums, aes(x=date, y=anomaly)) +
  geom_vline(xintercept = as.numeric(as.Date("2013-04-01")), linetype=5) +
  geom_vline(xintercept = as.numeric(as.Date("2016-03-25")), linetype=5) +
  geom_vline(xintercept = as.numeric(as.Date("2013-06-01")), linetype=3) +
  geom_vline(xintercept = as.numeric(as.Date("2016-06-20")), linetype=3) +
  # geom_smooth(data = output_data, aes(y=z,
  #                                   x=date, color = variable, group = paste(variable, plot)),
  #           method = "lm", formula=y~poly(x,21),
  #           se=F,
  #           alpha = 0.95)+
  geom_hline(yintercept = 0)+
  geom_line(data = output_data%>% filter(variable == "tmin"), aes(y=z_anom,
                                    x=date, color = variable, group = paste(variable, plot)),
            alpha =0.2)+
  ylab(NULL) +
  xlab(NULL) +
  #annotate("text",x= 2013, y=213, label = "30 year normal = 210 mm") +
  theme_pubr() +
  # scale_y_continuous(breaks = c(-1,0,2,4))+
  scale_x_date(date_labels = "%Y") +
  scale_color_manual(values = "black", labels= "Minimum Temperature") +
  theme(axis.text.x = element_blank()) +
  theme(legend.position = "none",
        legend.justification = c(0,1),
        legend.title=element_blank(),
        legend.background = element_rect(fill = "transparent"))+
  guides(colour = guide_legend(override.aes = list(alpha = 1)))

panel <- ggarrange(oc1,oc2,oc3,pp1, nrow=4, 
                   labels = c("(a) Climatic Water Deficit",
                              "(b) Actual Evapotranspiration", 
                              "(c) Minimum Temperature", 
                              "(d) Precipitation"),
                   label.x = 0.05,hjust = "left") %>%
  annotate_figure(left= text_grob("Standarized Values", rot = 90))

ggsave("figures/figure_1_climate.pdf",plot = panel, limitsize = F, dpi = 300,
       width = 10, height = 7.5) 



# Figure 2: nmds - fg by plot ==================================================

fgs <- all_p %>%
  dplyr::select(AIG, AIF, PNG, PNF, ANF) %>%
  wisconsin()

nms <- metaMDS(fgs, trymax = 10000, wascores = TRUE) %>%
  metaMDS(fgs, previous.best = ., trymax = 10000, wascores = TRUE)

scores0 <- as.data.frame(scores(nms)) %>%
  mutate(site_number = all_p$Site_number,
         year = all_p$Year,
         site_type = factor(all_p$Site.type, order = c(3,4,2,1)))

fgs <- fgs %>%
  mutate(Annuals = AIF+AIG+ANF,
         Perennials = PNG+PNF,
         Forbs = AIF+ANF+PNF,
         Grasses = PNG+AIG)

ef <- envfit(nms, fgs, na.rm = T, permutations = 9999)
env <- envfit(nms,  all_p%>%
                dplyr::select(starts_with("SOIL"), 
                              starts_with("Litter_"),
                              starts_with("Other"),
                              starts_with("soil")), 
              permutations = 9999)

poa <- envfit(nms,  all_p %>% dplyr::select(starts_with("Poa_")),
              permutations = 9999, na.rm=T)
bro <- envfit(nms,  all_p %>% dplyr::select(starts_with("Bromus_")),
              permutations = 9999, na.rm=T)

sp <-as.data.frame(ef$vectors$arrows*sqrt(ef$vectors$r))
spe <-as.data.frame(env$vectors$arrows*sqrt(env$vectors$r)) %>%
  rbind(as.data.frame(poa$vectors$arrows*sqrt(poa$vectors$r))) %>%
  rbind(as.data.frame(bro$vectors$arrows*sqrt(bro$vectors$r)))


species <- as.data.frame(cbind(sp, p=ef$vectors$pvals)) %>%
  tibble::rownames_to_column("species")  %>%
  filter(p < 0.05)

climscores <- as.data.frame(cbind(as.data.frame(env$vectors$arrows*sqrt(env$vectors$r)),
                                  p=env$vectors$pvals)) %>%
  tibble::rownames_to_column("species") %>%
  slice(13:16) %>%
  mutate(species = recode_factor(species,`ppt_2yr` = 'P[ant]', 
                                `sd_def` = 'sigma[CWD]',
                                `mean_aet` = 'AET',
                                `mean_tmin` = 'T[min]'))

vars <- rbind(as.data.frame(cbind(spe, p=c(env$vectors$pvals,
                                           poa$vectors$pvals,
                                           bro$vectors$pvals),
                                  r2 =c(env$vectors$r %>% as.numeric(),
                                        poa$vectors$r %>% as.numeric(),
                                        bro$vectors$r %>% as.numeric())))) %>%
  tibble::rownames_to_column("Variables")  %>%
  filter(p < 0.05)%>%
  mutate(Variables = lut_variables_nounit[Variables])



ef_table <- data.frame(Variables = ef$vectors$arrows %>% rownames(),
                       NMDS1 = ef$vectors$arrows[,1] %>% as.numeric(),
                       NMDS2 = ef$vectors$arrows[,2] %>% as.numeric(),
                       r2 = ef$vectors$r %>% as.numeric(),
                       p = ef$vectors$pvals) %>%
  filter(p<0.05) 


write_csv(rbind(ef_table, vars), "figures/envfit.csv")


scores <- mutate(scores0, i_y = str_c(site_type, ", ", year),
                 site_type = factor(site_type, 
                                    levels = c("Intact Sagebrush",
                                               "Invaded Sagebrush",
                                               "Cheatgrass-dominated", 
                                               "Cheatgrass Dieoff"), 
                                    labels = c("I. Intact Sagebrush",
                                               "II. Invaded Sagebrush",
                                               "III. Cheatgrass-dominated", 
                                               "IV. Cheatgrass Dieoff")))

mean_change <- scores %>%
  group_by(site_type,year) %>%
  dplyr::summarise(NMDS1_mean = mean(NMDS1),
                   NMDS2_mean = mean(NMDS2)) %>%
  ungroup()

# legend only
op1<- ggplot(data=scores) +
  geom_point(aes(x=NMDS1, y=NMDS2, color = site_type, shape = site_type),
             size= 4, stroke =1.5) +
  scale_shape_manual(values = c(17,18, 16,15),name = "Invasion Stage, Year") +
  scale_color_manual(values = c("springgreen4","deepskyblue2",
                                "darkgoldenrod2","chocolate4"),
                     name = "Invasion Stage, Year") +
  theme_bw() +
  theme(legend.title = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  theme(legend.justification=c(0,0), legend.position=c(0,0),
        legend.background = element_rect(fill = 'transparent'),
        legend.text = element_text(size=20,face="plain"));op1
nmds_legend<- get_legend(op1)

# mean change
p1 <- ggplot(data=scores) +
  geom_point(aes(x=NMDS1, y=NMDS2, color = i_y, shape = i_y),
             size= 4, stroke =1.5) +
  geom_path(data = mean_change,
            aes(group=site_type, x=NMDS1_mean, y=NMDS2_mean, color = site_type),
            arrow = arrow(angle = 20, length=unit(0.2, "inches")), lwd = 2,show.legend = F)+
  scale_shape_manual(values = c(15,0,16,1,17,2,18,5), name = "Invasion Stage") +
  scale_color_manual(values = c( "chocolate4", "chocolate4","darkgoldenrod2",
                                 "darkgoldenrod2","springgreen4","deepskyblue2",
                                 "darkgoldenrod2", "springgreen4","springgreen4",
                                 "deepskyblue2", "deepskyblue2", "chocolate4"),
                     name = "Invasion Stage") +
  xlim(c(-1.05,1.05))+
  ylim(c(-0.95,1.05))+
  theme_pubr() +
  coord_fixed() +
  annotate("text", x=-1.05, y=-0.9, fontface = "bold",
           label="Closed Symbols: 2013\nOpen Symbols: 2016",
           hjust = "left", size = 5.5)+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1),
        legend.justification=c(0,1), 
        axis.title = element_text(size=17, face="bold"),
        legend.position="none",
        legend.background = element_rect(fill = 'transparent'));p1

# climate fig
p_clim <- ggplot(data=scores) +
  geom_point(aes(x=NMDS1, y=NMDS2, color = i_y, shape = i_y),
             size= 4, stroke =1.5, alpha = 0.25) +
  geom_segment(data = climscores[1,],x=0,y=0,arrow = arrow(), color = "blue",
               aes(yend = NMDS2, xend = NMDS1), lwd=1)+
  geom_segment(data = climscores[2,],x=0,y=0,arrow = arrow(), color = "black",
               aes(yend = NMDS2, xend = NMDS1), lwd=1)+
  geom_segment(data = climscores[3,],x=0,y=0,arrow = arrow(), color = "red",
               aes(yend = NMDS2, xend = NMDS1),lwd=1)+
  geom_segment(data = climscores[4,],x=0,y=0,arrow = arrow(), color = "darkgoldenrod2",
               aes(yend = NMDS2, xend = NMDS1),lwd=1)+
  geom_text_repel(data=climscores, parse = T,
            aes(x=NMDS1,y=NMDS2,label=species),# hjust = "left",
            size=5,nudge_y = c(0.05,0,0,0), nudge_x = -.1,
            fontface = "bold")+
  scale_shape_manual(values = c(15,0,16,1,17,2,18,5),name = "Invasion Stage") +
  scale_color_manual(values = c( "chocolate4", "chocolate4",
                                 "darkgoldenrod2","darkgoldenrod2",
                                 "springgreen4", "springgreen4",
                                 "deepskyblue2", "deepskyblue2"),
                     name = "Invasion Stage") +
  xlim(c(-1.05,1.05))+
  ylim(c(-0.95,1.05))+
  theme_pubr()+
  coord_fixed() +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1),
        legend.justification=c(0,1), 
        axis.title = element_text(size=17, face="bold"),
        legend.position="none",
        legend.background = element_rect(fill = 'transparent'));p_clim

p3<-ggplot(data=scores) +
  geom_point(aes(x=NMDS1, y=NMDS2, color = i_y, shape = i_y),
             size= 4, stroke =1.5, alpha = 0.25) +
  geom_segment(data = species,x=0,y=0, color = "black",arrow = arrow(),
               aes(yend = NMDS2, xend = NMDS1), lwd=1)+
  geom_text_repel(data=species %>% filter(NMDS2>0),aes(x=NMDS1,y=NMDS2,label=species),
            size=5, color = "black",nudge_y = 0.05, fontface = "bold")+
  geom_text_repel(data=species %>% filter(NMDS2<=0),aes(x=NMDS1,y=NMDS2,label=species),
            size=5, color = "black",nudge_y = -0.05, fontface = "bold")+
  scale_shape_manual(values = c(15,0,16,1,17,2,18,5),name = "Invasion Stage") +
  scale_color_manual(values = c( "chocolate4", "chocolate4",
                                 "darkgoldenrod2","darkgoldenrod2",
                                 "springgreen4", "springgreen4",
                                 "deepskyblue2", "deepskyblue2"),
                     name = "Invasion Stage") +
  theme_pubr() +
  xlim(c(-1.05,1.05))+
  ylim(c(-0.95,1.05))+
  coord_fixed()+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1),
        legend.justification=c(0,1), 
        axis.title = element_text(size=17, face="bold"),
        legend.position="none",
        legend.background = element_rect(fill = 'transparent'));p3

# soil and plant nutrients
p4<-ggplot(data=scores) +
  geom_point(aes(x=NMDS1, y=NMDS2, color = i_y, shape = i_y),
             size= 4, stroke =1.5, alpha = 0.25) +
  geom_segment(data = vars,x=0,y=0, color = "black", arrow = arrow(),lwd=1,
               aes(yend = NMDS2, xend = NMDS1))+
  geom_richtext(data=vars %>%
                  mutate(Variables = str_replace_all(Variables,c("Bromus"="*B. tectorum*",
                                                             "Poa" = "*P. secunda*"))),
                  aes(x=NMDS1,y=NMDS2,label=Variables), hjust = "left",
                  label.color = "transparent", # remove background and outline
                  label.padding = grid::unit(rep(0, 4), "pt"),
                  size=5, color = "black",
                 nudge_y = c(  0, .1,   0.05, 0,   0,         0, -.05, -.05, -.02), 
                 nudge_x = c(-.5,  0, -.4, .05, 0.05, 0.05,  .05, -.7,  0.05),
                  fontface = "bold")+
  scale_shape_manual(values = c(15,0,16,1,17,2,18,5),name = "Invasion Stage") +
  scale_color_manual(values = c( "chocolate4", "chocolate4",
                                 "darkgoldenrod2","darkgoldenrod2",
                                 "springgreen4", "springgreen4",
                                 "deepskyblue2", "deepskyblue2"),
                     name = "Invasion Stage") +
  theme_pubr() +
  xlim(c(-1.05,1.05))+
  ylim(c(-0.95,1.05))+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1),
        legend.justification=c(0,1), 
        legend.position="none",
        axis.title = element_text(size=17, face="bold"),
        legend.background = element_rect(fill = 'transparent'))+
  coord_fixed();p4

p2<-ggdraw(p1) +
  draw_plot(nmds_legend, x=0.12, y=.76, height = 0.5, width = 0.4)

big_p<-ggarrange(p2,p3, p4,font.label = list(size=30), 
          nrow = 1, ncol=3,
          labels = c("(a)", "(b)", "(c)", "(d)"),
          label.x = 0.85, label.y=0.95,
          widths = c(2.5,2.5, 2.5,2.5)) 
ggsave(big_p,filename = "figures/fig_2_nmds_3panel_final.pdf",
         width = 19, height = 6.5, bg="white", dpi=300)

# Figures 3 & 4: site means by year and site type =========================
# 

# first we loop through building the linear mixed model to look at whether 
# differences in means are significant and do a post-hoc test

resp <- c("Bromus_TN_pct",    
          "Bromus_TC_pct", "Bromus_CN", "Other_TN_pct", "Other_TC_pct",     
          "Other_CN", "Poa_TN_pct", "Poa_TC_pct", "Poa_CN",    
          "Litter_TN_pct", "Litter_TC_pct", "Litter_CN", "SOIL_SurSo4_kg_ha",
          "SOIL_Ca_kg_ha", "SOIL_Mg_kg_ha" ,
          "SOIL_CN", "soil_n_kg_ha", "soil_c_kg_ha", "total_mineral_n",
          "NO3_kg_ha","NH4_kg_ha")
ss=list()
yy=list()
sy=list()
ssyy <-list()
mods = list()
coefs <- list()
ints<-c()
aovs<- list()
for(i in 1:length(resp)){
  # getting rid of NAs in the response variable
  d <- filter(all, is.na(dplyr::select(all, resp[i]))==F )#;dim(d)
 
  # making the mixed mods
  f<-formula(paste(resp[i],"~ Year*Site.type + (1|Site_number)"))
 
  mods[[i]] <- lmer(f,
             data=d, na.action = na.omit,# contrasts = T,
             REML = TRUE)
  aovs[[i]] <- anova(mods[[i]], ddf="Kenward-Roger") %>%
      as_tibble(rownames = "variable") %>%
      mutate(response = lut_variables[resp[i]])

  mods[[i]]@call$formula <- f
  
  if(aovs[[i]][3,"Pr(>F)"] < 0.05) {
    int <- TRUE
    ints[i] <- TRUE
  }else{
    int <- FALSE
    ints[i] <- FALSE
  }

  # making an output table
  cc <- summary(mods[[i]])

  coefs[[i]] <- cc$coefficients %>%
    as.data.frame() %>%
    tibble::rownames_to_column("variable")%>%
    dplyr::rename(p = "Pr(>|t|)", t="t value", se = "Std. Error")%>%
    dplyr::mutate(pr = round(p, 3),
                  response = resp[i],
                  ps = ifelse(.$p<0.001,"***", 
                              ifelse(.$p<0.01& .$p>0.001, "**",
                                     ifelse(.$p<0.05 & .$p>0.01, "*", "")))) %>%
    mutate(val_se = paste0(round(Estimate, 2)," (",
                           round(se, 2),")",
                           ps))
  
  # performing the post-hoc test in different ways
  marginal_y <- lsmeans(mods[[i]], ~ Year) %>% cld(alpha   = 0.05,
                                                               Letters = letters,
                                                               adjust  = "tukey",
                                                   decreasing = T)
  marginal_s <- lsmeans(mods[[i]], ~ Site.type) %>% cld(alpha   = 0.05,
                                                        Letters = letters,
                                                        adjust  = "tukey",
                                                        decreasing = T)
  marginal_sy <- lsmeans(mods[[i]], ~ Site.type:Year) %>% cld(alpha   = 0.05,
                                                              Letters = letters,
                                                              adjust  = "sidak",
                                                              decreasing = T)
  marginal_z1 <- lsmeans(mods[[i]], ~ Site.type|Year) %>% cld(alpha   = 0.05,
                                                              Letters = letters,
                                                              adjust  = "sidak",
                                                              decreasing = T)
  marginal_z2 <- lsmeans(mods[[i]], ~ Year|Site.type) %>% cld(alpha   = 0.05,
                                                              Letters = letters,
                                                              adjust  = "sidak",
                                                              decreasing = T)
  yy[[i]] <- marginal_y %>% as.data.frame() %>% mutate(variable = resp[i], int =int)
  ss[[i]]<- marginal_s %>% as.data.frame()  %>% mutate(variable = resp[i], int=int)
  sy[[i]]<- marginal_sy %>% as.data.frame()  %>% mutate(variable = resp[i], int=int)
  
  # for the interactions, we're just gonna do the post-hoc test on each
  # year separately from the beginning, so that way we dont have a million 
  # group letters
  
  if(!int){
    ssyy[[i]] <- marginal_sy %>% as.data.frame()  %>% mutate(variable = resp[i], int=int)
  }else{
    fi<-formula(paste(resp[i],"~0 + Site.type + (1|Site_number)"))
    
    s13 <- lmer(fi,
               data=d %>% filter(Year == "2013"), 
               na.action = na.omit, contrasts = T,
               REML = TRUE) %>%
      lsmeans(~ Site.type) %>% 
      cld(alpha   = 0.05,
          Letters = letters,
          adjust  = "tukey")%>% 
      as.data.frame() %>% 
      mutate(variable = resp[i], int =int, Year = "2013",
             .group = str_to_upper(.group))
    s16 <- lmer(fi,
               data=d %>% filter(Year == "2016"), 
               na.action = na.omit, contrasts = T,
               REML = T) %>%
      lsmeans( ~ Site.type) %>% 
      cld(alpha   = 0.05,
          Letters = letters,
          adjust  = "tukey")%>% 
      as.data.frame() %>% 
      mutate(variable = resp[i], int =int, Year = "2016")
    
    ssyy[[i]] <- rbind(s13, s16)
    

  }
  print(i)
}

do.call("rbind",aovs[1:9]) %>%
  saveRDS("data/aov_tab_plnt.RDS")
do.call("rbind",aovs[10:21]) %>%
  saveRDS("data/aov_tab_soil.RDS")

saveRDS(mods[ints], "data/mods_i.RDS")
saveRDS(mods[!ints], "data/mods_ni.RDS")

cc <- do.call("rbind", coefs) %>%
  dplyr::select(variable, val_se, response, df) %>%
  spread(key = variable, value = val_se, fill = "") %>%
  dplyr::select(-Year2016, -Year2013) %>%
  mutate(response = lut_variables[response])

do.call("rbind", yy) %>%
  mutate(mean_se =str_c(round(lsmean,2)," (",round(SE, 2),") ", trimws(.group))) %>%
  dplyr::select(Year, mean_se, variable) %>%
  pivot_wider(id_cols = variable, names_from = Year, values_from = mean_se) %>%
  write_csv("figures/year_emmeans.csv")

site_means <- do.call("rbind", ss)

site_means %>% 
  group_by(variable) %>%
  mutate(n_groups = length(unique(.group))) %>%
  ungroup() %>%
  mutate(.group = ifelse(n_groups ==1, "", .group)) %>%
  mutate(mean_se =str_c(round(lsmean,2)," (",round(SE, 2),") ", trimws(.group))) %>%
  dplyr::select(Site.type, mean_se, variable) %>%
  pivot_wider(id_cols = variable, names_from = Site.type, values_from = mean_se) %>%
  write_csv("figures/site_emmeans.csv")

site_x_year_means <- do.call("rbind", ssyy) %>%
  group_by(variable) %>%
  mutate(group = str_to_lower(.group) %>% trimws,
         keep = ifelse(length(unique(group))==1,F,T)) %>%
  ungroup%>%
  filter(keep) %>%
  dplyr::select(-keep, -group)
  

site_x_year_means%>%
  dplyr::filter(int == TRUE) %>%
  mutate(mean_se =str_c(round(lsmean,2)," (",round(SE, 2),") ", trimws(.group))) %>%
  dplyr::select(Year, Site.type, mean_se, variable) %>%
  pivot_wider(names_from = variable, id_cols = c(Site.type,Year), values_from = mean_se) %>%
  write_csv("figures/site_x_year_emmeans_split.csv")


site_x_year_means_long <- do.call("rbind", sy) %>%
  dplyr::filter(int == TRUE) %>%
  mutate(mean_se =str_c(round(lsmean,2)," (",round(SE, 2),") ", trimws(.group))) %>%
  dplyr::select(Year, Site.type, mean_se, variable) %>%
  pivot_wider(names_from = variable, id_cols = c(Site.type,Year), values_from = mean_se) %>%
  write_csv("figures/site_x_year_emmeans.csv")


site_groups <- site_means %>%
  dplyr::select(Site.type, variable, s_group = .group, int)

sxy_groups <- site_x_year_means %>%
  dplyr::select(Site.type, Year, variable, sy_group = .group)

write_csv(cc, "/home/a/projects/Jones_Study/figures/site_means_mods.csv")

# now plotting =============================================================
lut_variables <- c("Bromus_TN_pct" = "(b) Bromus N (%)",
                   "Bromus_TC_pct" = "(a) Bromus C (%)",
                   "Bromus_CN" = "(c) Bromus C:N",
                   "Poa_TC_pct" = "(d) Poa C (%)",
                   "Poa_TN_pct" = "(e) Poa N (%)",
                   "Poa_CN" = "(f) Poa C:N",
                   "Other_TN_pct" = "(h) Other N (%)",
                   "Other_TC_pct" = "(g) Other C (%)",
                   "Other_CN" = "(i) Other C:N",
                   "Litter_TN_pct" = "(h) Litter N (%)",
                   "Litter_TC_pct" = "(g) Litter C (%)",
                   "Litter_CN" = "(i) Litter C:N",
                   "SOIL_SurSo4_kg_ha" = "Soil SurSo4 (kg/ha)",
                   "SOIL_Ca_kg_ha" = "(f) Soil Ca (kg/ha)",
                   "SOIL_Mg_kg_ha" = "(e) Soil Mg (kg/ha)",
                   "SOIL_CN" = "(c) Soil C:N",
                   "soil_n_kg_ha" = "(b) Soil Total N (kg/ha)",
                   "soil_c_kg_ha" = "(a) Soil Total C (kg/ha)",
                   "total_mineral_n" = "(d) Soil Mineral N (kg/ha)",
                   "NO3_kg_ha" = "Soil Nitrate (kg/ha)",
                   "NH4_kg_ha" = "Soil Ammonium (kg/ha)")

lut_col <- c("springgreen4", "deepskyblue2", "darkgoldenrod2", "chocolate4")

# lll <- do.call('rbind', ll)%>%
#   dplyr::select(group = .group, Site.type, Year, variable)

# code to fix the scales for the individual facets
# https://stackoverflow.com/questions/18046051/setting-individual-axis-limits-with-facet-wrap-and-scales-free-in-ggplot2#21585521


plant_groups <- all_p %>%
  dplyr::select(Site_number, Year, Site.type,
                Bromus_TN_pct, Bromus_TC_pct, Bromus_CN,
                Poa_TN_pct, Poa_TC_pct, Poa_CN,
                Other_TN_pct, Other_TC_pct, Other_CN) %>%
  gather(key =variable,value = value, -Site.type, -Year, -Site_number) %>%
  # left_join(site_groups, by = c("Site.type", "variable")) %>%
  left_join(sxy_groups, by = c("Site.type", "variable", "Year"))%>%
  # mutate(gg_group = ifelse(int == TRUE, sy_group, s_group)%>%trimws) %>%
  mutate(gg_group =  sy_group%>%trimws) %>%
  mutate(variable = factor(lut_variables[variable],
                           levels = c("(a) Bromus C (%)",
                                      "(b) Bromus N (%)",
                                      "(c) Bromus C:N",
                                      "(d) Poa C (%)",
                                      "(e) Poa N (%)",
                                      "(f) Poa C:N",
                                      "(g) Other C (%)",
                                      "(h) Other N (%)",
                                      "(i) Other C:N"),
                           labels=c("(a)~italic(B.~tectorum)~C~('%')",
                                    "(b)~italic(B.~tectorum)~N~('%')",
                                    "(c)~italic(B.~tectorum)~C:N",
                                    "(d)~italic(P.~secunda)~C~('%')",
                                    "(e)~italic(P.~secunda)~N~('%')",
                                    "(f)~italic(P.~secnuda)~C:N",
                                    "(g)~Other~C~('%')"  ,
                                    "(h)~Other~N~('%')"  ,
                                    "(i)~Other~C:N")),
         Site.type = factor(Site.type, 
                            levels = c("Intact Sagebrush",
                                       "Invaded Sagebrush",
                                       "Cheatgrass-dominated",
                                       "Cheatgrass Dieoff"),
                            labels = c("I. Intact Sagebrush", 
                                       "II. Invaded Sagebrush",
                                       "III. Cheatgrass-dominated", 
                                       "IV. Cheatgrass Dieoff")))%>%
  group_by(variable, Site.type,Year) %>%
  dplyr::summarise(vjust = boxplot.stats(na.omit(value))$stats[5],
            value = median(value, na.rm=T),
            gg_group = first(gg_group)
            ) %>%
  ungroup() %>%
  mutate(x_position = rep(c(0.63,1.63, 0.82,1.82, 1.01,2.01, 1.2,2.2), 
                          length.out = nrow(.))) %>%
  group_by(variable) %>%
  mutate(num_grps = length(unique(gg_group))) %>%
  ungroup() %>%
  group_by(variable) %>%
  mutate(gg_group = replace(gg_group, num_grps ==1, "")) %>%
  ungroup()

vjusts<- c(40, 40.8, 42, 41.7, #brom c
           43.3, 44, 43.5, 43.7,
           1.55,  0.9,  1.9,  1.4, #poa n
           0.58,  0.64,  1,  0.7,
           41.8, 43.5, 44.6, 42.2, #oc
           42.1, 43.5, 46.4, 46.3,
           1.6,  1.5,  2.7,  3.1, #on
           0.9,  0.85,  2.5,  2.2,
           37, 44, 37.5, 22, #ocn
           58, 80, 51, 37)

letter_positions <- plant_groups %>%
  dplyr::select(variable, Year, x_position, vjust,gg_group, Site.type) %>%
  filter(gg_group != "") %>%
  arrange(variable, Year) %>%
  mutate(x_position = replace(x_position, Site.type=="I. Intact Sagebrush" &
                                Year == "2013", 0.66))%>%
  mutate(x_position = replace(x_position, Site.type=="II. Invaded Sagebrush" &
                                Year == "2013", 0.86))%>%
  mutate(x_position = replace(x_position, Site.type=="III. Cheatgrass-dominated" &
                                Year == "2013", 1.035))%>%
  mutate(x_position = replace(x_position, Site.type=="IV. Cheatgrass Dieoff" &
                                Year == "2013", 1.235))%>%
  mutate(x_position = replace(x_position, Site.type=="I. Intact Sagebrush" &
                                Year == "2016", 1.67))%>%
  mutate(x_position = replace(x_position, Site.type=="II. Invaded Sagebrush" &
                                Year == "2016", 1.86))%>%
  mutate(x_position = replace(x_position, Site.type=="III. Cheatgrass-dominated" &
                                Year == "2016", 2.04))%>%
  mutate(x_position = replace(x_position, Site.type=="IV. Cheatgrass Dieoff" &
                                Year == "2016", 2.24))%>%
  mutate(vjust = vjusts)

letter_positions[27, "x_position"] <- 1.03
# letter_positions[11, "x_position"] <- 1

dummy <- letter_positions %>%
   filter(variable == "(a)~italic(B.~tectorum)~C~('%')" | 
            variable == "(e)~italic(P.~secunda)~N~('%')"|
            variable == "(g)~Other~C~('%')" | variable== "(h)~Other~N~('%')"|
            variable == "(i)~Other~C:N")
dummy[1,"vjust"] <- 46
dummy[9, "vjust"] <- 2
dummy[17, "vjust"] <- 48
dummy[25, "vjust"] <- 4
dummy[33, "vjust"] <- 90

all_p %>%
  dplyr::select(Site_number, Year, Site.type,
                Bromus_TN_pct, Bromus_TC_pct, Bromus_CN,
                Poa_TN_pct, Poa_TC_pct, Poa_CN,
                Other_TN_pct, Other_TC_pct, Other_CN) %>%
  gather(key =variable,value = value, -Site.type, -Year, -Site_number) %>%
  mutate(variable = factor(lut_variables[variable],
                           levels = c("(a) Bromus C (%)",
                                      "(b) Bromus N (%)",
                                      "(c) Bromus C:N",
                                      "(d) Poa C (%)",
                                      "(e) Poa N (%)",
                                      "(f) Poa C:N",
                                      "(g) Other C (%)",
                                      "(h) Other N (%)",
                                      "(i) Other C:N"),
                           labels=c("(a)~italic(B.~tectorum)~C~('%')",
                                    "(b)~italic(B.~tectorum)~N~('%')",
                                    "(c)~italic(B.~tectorum)~C:N",
                                    "(d)~italic(P.~secunda)~C~('%')",
                                    "(e)~italic(P.~secunda)~N~('%')",
                                    "(f)~italic(P.~secnuda)~C:N",
                                    "(g)~Other~C~('%')"  ,
                                    "(h)~Other~N~('%')"  ,
                                    "(i)~Other~C:N")),
         Site.type = factor(Site.type, 
                            levels = c("Intact Sagebrush", 
                                       "Invaded Sagebrush",
                                       "Cheatgrass-dominated", 
                                       "Cheatgrass Dieoff"),
                            labels = c("I. Intact Sagebrush", 
                                       "II. Invaded Sagebrush",
                                       "III. Cheatgrass-dominated", 
                                       "IV. Cheatgrass Dieoff")))%>%
  # mutate(variable = factor(variable,labels=c("(a)~italic(B.~tectorum)~C~('%')",
  #                                            "(b)~italic(B.~tectorum)~N~('%')",
  #                                            "(c)~italic(B.~tectorum)~C:N",
  #                                            "(d)~italic(P.~secunda)~C~('%')",
  #                                            "(e)~italic(P.~secunda)~N~('%')",
  #                                            "(f)~italic(P.~secnuda)~C:N",
  #                                            "(g)~Other~C~('%')"  ,
  #                                            "(h)~Other~N~('%')"  ,
  #                                            "(i)~Other~C:N")))%>%
  ggplot() +
    geom_boxplot(aes(x=Year, y = value, fill = Site.type),
                 outlier.shape = NA) +
    scale_fill_manual(name="Invasion Stage", values = lut_col)+
    geom_text(data = letter_positions, vjust = 0,hjust = "left",
              aes(x=x_position, y = vjust,
                  group = Site.type, label = gg_group)) +
    facet_wrap(~variable, scales = "free", labeller=label_parsed) +
    theme_pubr()+
    theme(legend.position = "bottom",
          legend.title = element_blank()) +
    geom_blank(data=dummy, aes(y=vjust))+
    xlab(NULL)+
    ylab(NULL)+
    guides(fill=guide_legend(nrow=2,byrow=TRUE)) +
    ggsave("figures/figure_4_tissue_raw.pdf", width = 8.25, height = 8.25) 

# soil data monger ==========
soil_groups <- all_p %>%
  dplyr::select(Site_number, Year, Site.type,
                soil_n_kg_ha, soil_c_kg_ha, SOIL_CN,
                total_mineral_n, SOIL_Mg_kg_ha, SOIL_Ca_kg_ha,
                Litter_TN_pct, Litter_TC_pct, Litter_CN) %>%
  gather(key =variable,value = value, -Site.type, -Year, -Site_number) %>%
  left_join(site_groups, by = c("Site.type", "variable")) %>%
  left_join(sxy_groups, by = c("Site.type", "variable", "Year"))%>%
  mutate(gg_group = ifelse(int == TRUE, sy_group, s_group)%>%trimws) %>%
  mutate(variable = factor(lut_variables[variable],
                           levels= c( "(a) Soil Total C (kg/ha)", 
                                      "(b) Soil Total N (kg/ha)",
                                      "(c) Soil C:N",
                                      "(d) Soil Mineral N (kg/ha)",
                                      "(e) Soil Mg (kg/ha)",
                                      "(f) Soil Ca (kg/ha)",
                                      "(g) Litter C (%)",  
                                      "(h) Litter N (%)",
                                      "(i) Litter C:N")),
         Site.type = factor(Site.type, 
                            levels = c("Intact Sagebrush", 
                                       "Invaded Sagebrush",
                                       "Cheatgrass-dominated", 
                                       "Cheatgrass Dieoff"),
                            labels = c("I. Intact Sagebrush", 
                                       "II. Invaded Sagebrush",
                                       "III. Cheatgrass-dominated", 
                                       "IV. Cheatgrass Dieoff")))%>%
  
  group_by(variable, Site.type,Year) %>%
  dplyr::summarise(vjust = boxplot.stats(value)$stats[4]+boxplot.stats(value)$stats[4]*.1,
            value = median(value, na.rm=T),
            gg_group = first(gg_group)
  ) %>%
  ungroup() %>%
  mutate(x_position = rep(c(0.63,1.63, 0.82,1.82, 1.01,2.01, 1.2,2.2), 
                          length.out = nrow(.))) %>%
  group_by(variable) %>%
  mutate(num_grps = length(unique(gg_group))) %>%
  ungroup() %>%
  group_by(variable) %>%
  mutate(gg_group = replace(gg_group, num_grps ==1, "")) %>%
  ungroup()

# soil plot ======================
letter_positions <- soil_groups %>%
  dplyr::select(variable, Year, x_position, vjust,gg_group, Site.type)%>%
  mutate(x_position = replace(x_position, Site.type=="I. Intact Sagebrush" &
                                Year == "2013", 0.66))%>%
  mutate(x_position = replace(x_position, Site.type=="II. Invaded Sagebrush" &
                                Year == "2013", 0.86))%>%
  mutate(x_position = replace(x_position, Site.type=="III. Cheatgrass-dominated" &
                                Year == "2013", 1.035))%>%
  mutate(x_position = replace(x_position, Site.type=="IV. Cheatgrass Dieoff" &
                                Year == "2013", 1.235))%>%
  mutate(x_position = replace(x_position, Site.type=="I. Intact Sagebrush" &
                                Year == "2016", 1.66))%>%
  mutate(x_position = replace(x_position, Site.type=="II. Invaded Sagebrush" &
                                Year == "2016", 1.86))%>%
  mutate(x_position = replace(x_position, Site.type=="III. Cheatgrass-dominated" &
                                Year == "2016", 2.04))%>%
  mutate(x_position = replace(x_position, Site.type=="IV. Cheatgrass Dieoff" &
                                Year == "2016", 2.24))%>%
  mutate(vjust = replace(vjust, Site.type=="I. Intact Sagebrush" &
                                Year == "2013" & 
                           variable == "B. Soil Total N (kg/ha)", 1700))%>%
  mutate(vjust = replace(vjust, Site.type=="II. Invaded Sagebrush" &
                           Year == "2013" & 
                           variable == "B. Soil Total N (kg/ha)", 1650))%>%
  mutate(vjust = replace(vjust, Site.type=="III. Cheatgrass-dominated" &
                           Year == "2013" & 
                           variable == "B. Soil Total N (kg/ha)", 1900))%>%
  mutate(vjust = replace(vjust, Site.type=="IV. Cheatgrass Dieoff" &
                           Year == "2013" & 
                           variable == "B. Soil Total N (kg/ha)", 1500))%>%
  mutate(vjust = replace(vjust, Site.type=="I. Intact Sagebrush" &
                           Year == "2016" & 
                           variable == "B. Soil Total N (kg/ha)", 1900))%>%
  mutate(vjust = replace(vjust, Site.type=="II. Invaded Sagebrush" &
                           Year == "2016" & 
                           variable == "B. Soil Total N (kg/ha)", 1350))%>%
  mutate(vjust = replace(vjust, Site.type=="III. Cheatgrass-dominated" &
                           Year == "2016" & 
                           variable == "B. Soil Total N (kg/ha)", 1500))%>%
  mutate(vjust = replace(vjust, Site.type=="IV. Cheatgrass Dieoff" &
                           Year == "2016" & 
                           variable == "B. Soil Total N (kg/ha)", 1400));letter_positions[9:22,]

vjusts<- c(1700,1650,1900,1450,
           1900,1250,1500,1290, 
           38.7,   36.5,   37.5,   40.5,
           34,   34, 39, 40,
           1.13,    0.83,    1.1036667,    0.8891667,
           1.25,    1.27,    1.03,    1.2)
letter_positions <- soil_groups %>%
  dplyr::select(variable, Year, x_position, vjust,gg_group, Site.type) %>%
  filter(gg_group != "") %>%
  arrange(variable, Year) %>%
  mutate(x_position = replace(x_position, Site.type=="I. Intact Sagebrush" &
                                Year == "2013", 0.66))%>%
  mutate(x_position = replace(x_position, Site.type=="II. Invaded Sagebrush" &
                                Year == "2013", 0.86))%>%
  mutate(x_position = replace(x_position, Site.type=="III. Cheatgrass-dominated" &
                                Year == "2013", 1.035))%>%
  mutate(x_position = replace(x_position, Site.type=="IV. Cheatgrass Dieoff" &
                                Year == "2013", 1.235))%>%
  mutate(x_position = replace(x_position, Site.type=="I. Intact Sagebrush" &
                                Year == "2016", 1.67))%>%
  mutate(x_position = replace(x_position, Site.type=="II. Invaded Sagebrush" &
                                Year == "2016", 1.86))%>%
  mutate(x_position = replace(x_position, Site.type=="III. Cheatgrass-dominated" &
                                Year == "2016", 2.04))%>%
  mutate(x_position = replace(x_position, Site.type=="IV. Cheatgrass Dieoff" &
                                Year == "2016", 2.24))%>%
  mutate(vjust = vjusts)
letter_positions[10, "x_position"] <- .82
letter_positions[11, "x_position"] <- 1

dummy <- letter_positions%>%
  filter(variable == "(b) Soil Total N (kg/ha)" | variable == "(g) Litter C (%)")
dummy[1,"vjust"] <- 2200
dummy[9, "vjust"] <- 45

all_p %>%
  dplyr::select(Site_number, Year, Site.type,
                soil_n_kg_ha, soil_c_kg_ha, SOIL_CN,
                total_mineral_n, SOIL_Mg_kg_ha, SOIL_Ca_kg_ha,
                Litter_TN_pct, Litter_TC_pct, Litter_CN) %>%
  gather(key =variable,value = value, -Site.type, -Year, -Site_number) %>%
  mutate(variable = factor(lut_variables[variable],
                           levels= c( "(a) Soil Total C (kg/ha)", 
                                      "(b) Soil Total N (kg/ha)",
                                      "(c) Soil C:N",
                                      "(d) Soil Mineral N (kg/ha)",
                                      "(e) Soil Mg (kg/ha)",
                                      "(f) Soil Ca (kg/ha)",
                                      "(g) Litter C (%)",  
                                      "(h) Litter N (%)",
                                      "(i) Litter C:N")),
         Site.type = factor(Site.type, 
                            levels = c("Intact Sagebrush", "Invaded Sagebrush",
                                       "Cheatgrass-dominated", "Cheatgrass Dieoff"),
                            labels = c("I. Intact Sagebrush", 
                                       "II. Invaded Sagebrush",
                                       "III. Cheatgrass-dominated", 
                                       "IV. Cheatgrass Dieoff")))%>%
  ggplot(aes(y = letter_positions$vjust)) +
  geom_boxplot(aes(x=Year, y = value, fill = Site.type),
               outlier.shape = NA) +
  geom_blank(data=dummy, aes(y=vjust))+
  scale_fill_manual(name="Invasion Stage", values = lut_col)+
  facet_wrap(~variable, scales = "free") +
  geom_text(data = letter_positions, vjust = 0,hjust = "left",
            aes(x=x_position, y = vjust, 
                group = Site.type, label = gg_group)) +
  scale_y_continuous(labels = scales::label_number_si(accuracy=.1))+
  theme_pubr()+
  theme(legend.position = "bottom",
        legend.title = element_blank()) +
  xlab(NULL)+
  ylab(NULL)+
  guides(fill=guide_legend(nrow=2,byrow=TRUE)) +
  ggsave("figures/figure_3_soil_raw.pdf", width = 8.25, height = 8.25) 

# alternate figure 3 ===========================================================

gps<- soil_groups %>% rbind(plant_groups) %>%
  mutate(variable = str_sub(variable, 5,length(variable))) %>%
  filter(variable %in% c("Soil Total C (kg/ha)","Soil Total N (kg/ha)",
                         "Litter C:N" ,"italic(B.~tectorum)~C:N" ,
                         "italic(P.~secnuda)~C:N", "Other~C:N" ))

all_p %>%
  dplyr::select(Site_number, Year, Site.type,Bromus_CN, Litter_CN,
                Poa_CN,Other_CN, soil_n_kg_ha, soil_c_kg_ha) %>%
  gather(key =variable,value = value, -Site.type, -Year, -Site_number) %>%
  mutate(variable = factor(lut_variables[variable],
                           levels = c("(a) Soil Total C (kg/ha)",
                                      "(b) Soil Total N (kg/ha)", 
                                        "(i) Litter C:N",
                                        "(c) Bromus C:N",
                                      "(f) Poa C:N",
                                      "(i) Other C:N"),
                           labels=c( "(a)~Soil~Total~C~(kg/ha)",
                                     "(b)~Soil~Total~N~(kg/ha)", 
                                     "(c)~Litter~C:N",
                                     "(d)~italic(B.~tectorum)~C:N",
                                    "(e)~italic(P.~secnuda)~C:N",
                                    "(f)~Other~C:N")),
         Site.type = factor(Site.type, 
                            levels = c("Intact Sagebrush", 
                                       "Invaded Sagebrush",
                                       "Cheatgrass-dominated", 
                                       "Cheatgrass Dieoff"),
                            labels = c("I. Intact Sagebrush", 
                                       "II. Invaded Sagebrush",
                                       "III. Cheatgrass-dominated", 
                                       "IV. Cheatgrass Dieoff")))%>%
  ggplot() +
  geom_boxplot(aes(x=Year, y = value, fill = Site.type),
               outlier.shape = NA) +
  scale_fill_manual(name="Invasion Stage", values = lut_col)+
  facet_wrap(~variable, scales = "free", labeller=label_parsed) +
  theme_pubr()+
  theme(legend.position = "bottom",
        legend.title = element_blank()) +
  xlab(NULL)+
  ylab(NULL)+
  guides(fill=guide_legend(nrow=2,byrow=TRUE)) +
  ggsave("figures/figure_3_new.png", width = 7.5, height = 6) 

# Figure supplement: bivariate predictors of soil =======================================
lut_col <- c("springgreen4", "deepskyblue2", "darkgoldenrod2", "chocolate4")
lut_st <- c("M" = "Invaded Sagebrush", "I" = "Intact Sagebrush",
            "C" = "Cheatgrass-dominated", "D" = "Cheatgrass Dieoff")

lut_variables <- c("Bromus_TN_pct" = "Bromus N (%)", "Bromus_TC_pct" = "Bromus C (%)",
                   "Bromus_CN" = "Bromus C:N", "Other_TN_pct" = "Other N (%)",
                   "Other_TC_pct" = "Other C (%)", "Other_CN" = "Other C:N", 
                   "Poa_TN_pct" = "Poa N (%)", "Poa_TC_pct" = "Poa C (%)",
                   "Poa_CN" = "Poa C:N", "Litter_TN_pct" = "Litter N (%)",
                   "Litter_TC_pct" = "Litter C (%)", "Litter_CN" = "Litter C:N",
                   "SOIL_SurSo4_kg_ha" = "Soil SurSo4 (kg/ha)",
                   "SOIL_Ca_kg_ha" = "Soil Ca (kg/ha)", 
                   "SOIL_Mg_kg_ha" = "Soil Mg (kg/ha)",
                   "SOIL_CN" = "Soil C:N", "soil_n_kg_ha" = "Soil Total N (kg/ha)",
                   "soil_c_kg_ha" = "Soil Total C (kg/ha)",
                   "total_mineral_n" = "Soil Mineral N (kg/ha)",
                   "NO3_kg_ha" = "Soil Nitrate (kg/ha)",
                   "NH4_kg_ha" = "Soil Ammonium (kg/ha)",
                   "ja_ju_tmin_min" = "Minimum Temperature",
                   "ja_ju_aet" = "Actual Evapotranspiration",
                   "AIG" = "Ann. Inv. Grass Cover (%)",
                   "AIF" = "Ann. Inv. Forb Cover (%)",
                   "Bare" = "Bare Ground Cover (%)")

fig5_df<-all_p %>%
  dplyr::select(AIG,AIF, Bare, ja_ju_aet,
                ja_ju_tmin_min, Site.type, Year, Site_number) %>%
  gather(key=variable, value=value,-Site.type,-Year,-Site_number) %>%
  mutate(variable = factor(lut_variables[variable], 
                           levels = c("Actual Evapotranspiration",
                                      "Minimum Temperature",
                                      "Bare Ground Cover (%)",
                                      "Ann. Inv. Grass Cover (%)",
                                      "Ann. Inv. Forb Cover (%)")),
         Site.type = factor(Site.type, levels = c("Intact Sagebrush",
                                                  "Invaded Sagebrush",
                                                  "Cheatgrass-dominated", 
                                                  "Cheatgrass Dieoff"))) %>%
  left_join(dplyr::select(all_p, soil_n_kg_ha, soil_c_kg_ha,Site_number, Year), 
            by = c("Site_number","Year")) %>%
  left_join(all_se %>% dplyr::select(n_se = soil_n_kg_ha, c_se = soil_c_kg_ha,
                Year, Site_number), by = c("Site_number","Year")) 

n_mods <- fig5_df %>%
  nest(-variable) %>%
  mutate(fit = map(data, ~ lmer(soil_n_kg_ha ~ value*Site.type +(value|Site.type), data=.)),
         results = map(fit, augment)) %>%
  unnest(results) %>%
  dplyr::rename(n_fit = .fitted) %>%
  dplyr::select(variable, soil_n_kg_ha, value, n_fit)

c_mods <- fig5_df %>%
  nest(-variable) %>%
  mutate(fit = map(data, ~ lmer(soil_c_kg_ha ~ value*Site.type +(value|Site.type), data=.)),
         results = map(fit, augment)) %>%
  unnest(results) %>%
  dplyr::rename(c_fit = .fitted)%>%
  dplyr::select(variable, soil_c_kg_ha, value, c_fit)

fig5_df1 <- left_join(fig5_df, n_mods, 
                      by = c("variable", "soil_n_kg_ha", "value")) %>%
left_join(c_mods, by = c("variable", "soil_c_kg_ha", "value"))

# making more facet-able
f5_df_n <- dplyr::select(fig5_df1, Site.type, Year, Site_number, variable, value,
                         s_value = soil_n_kg_ha, se = n_se, fit = n_fit) %>%
  mutate(s_variable = "Soil Total Nitrogen (kg/ha)")
f5_df_c <- dplyr::select(fig5_df1, Site.type, Year, Site_number, variable, value,
                         s_value = soil_c_kg_ha, se = c_se, fit = c_fit) %>%
  mutate(s_variable = "Soil Total C (kg/ha)")

f5_both <- rbind(f5_df_c, f5_df_n)

f5 <- ggplot(f5_both, aes(x = value, y = s_value, color = Site.type,
                          shape = Site.type)) +
  geom_line(aes(y=fit), show.legend = FALSE) +
  geom_point() +
  geom_errorbar(aes(x = value,
                    ymin = s_value - se,
                    ymax = s_value + se,
                    color=Site.type),
                alpha = 0.5, show.legend = F) +
  scale_color_manual(name="Invasion Stage", values = lut_col) +
  scale_shape_manual(name = "Invasion Stage", values = c(15,16,17,18)) +
  facet_grid(s_variable~variable, scales="free") +
  theme_bw() +
  theme(axis.title = element_blank())+
  theme(legend.position = "bottom")+
  ggsave("figures/fig_5_soil_univariate.png", dpi = 600, width = 10, height =6.7)

# maybe make a table too


# fig supplement plant tissue scatter plots =============================================

fig6_df<-all_p %>%
  dplyr::select(AIG,AIF, ja_ju_aet,
                ja_ju_tmin_min, Site.type, Year, Site_number) %>%
  gather(key=variable, value=value,-Site.type,-Year,-Site_number) %>%
  mutate(variable = factor(lut_variables[variable], 
                           levels = c("Actual Evapotranspiration",
                                      "Minimum Temperature",
                                      "Ann. Inv. Grass Cover (%)",
                                      "Ann. Inv. Forb Cover (%)")),
         Site.type = factor(Site.type, levels = c("Intact Sagebrush",
                                                  "Invaded Sagebrush",
                                                  "Cheatgrass-dominated", 
                                                  "Cheatgrass Dieoff"))) %>%
  left_join(dplyr::select(all_p, Bromus_CN, Poa_CN, Other_CN,Site_number, Year), 
            by = c("Site_number","Year")) %>%
  left_join(all_se %>% dplyr::select(p_se = Poa_CN, b_se = Bromus_CN, o_se = Other_CN,
                                     Year, Site_number), by = c("Site_number","Year")) 

b_mods <- fig6_df %>%
  nest(-variable) %>%
  mutate(fit = map(data, ~ lmer(Bromus_CN ~ value*Site.type +(value|Site.type), data=.)),
         results = map(fit, augment)) %>%
  unnest(results) %>%
  dplyr::rename(b_fit = .fitted) %>%
  dplyr::select(variable, Bromus_CN, value, b_fit)

b_p <- fig6_df %>%
  nest(-variable) %>%
  mutate(fit = map(data, ~ lmer(Bromus_CN ~ value*Site.type +(value|Site.type), data=.)),
         results = map(fit, summary)) %>%
  unnest(results)

p_mods <- fig6_df %>%
  filter(is.nan(Poa_CN) == FALSE) %>%
  nest(-variable) %>%
  mutate(fit = map(data, ~ lmer(Poa_CN ~ value*Site.type +(value|Site.type), data=.)),
         results = map(fit, augment)) %>%
  unnest(results) %>%
  dplyr::rename(p_fit = .fitted)%>%
  dplyr::select(variable, Poa_CN, value, p_fit)

o_mods <- fig6_df %>%
  filter(is.nan(Other_CN) == FALSE) %>%
  nest(-variable) %>%
  mutate(fit = map(data, ~ lmer(Other_CN ~ value*Site.type +(value|Site.type), data=.)),
         results = map(fit, augment)) %>%
  unnest(results) %>%
  dplyr::rename(o_fit = .fitted)%>%
  dplyr::select(variable, Other_CN, value, o_fit)

fig6_df1 <- left_join(fig6_df, b_mods, 
                      by = c("variable", "Bromus_CN", "value")) %>%
  left_join(p_mods, by = c("variable", "Poa_CN", "value")) %>%
  left_join(o_mods, by = c("variable", "Other_CN", "value"))

# makeing more facetable
f6_b <- dplyr::select(fig6_df1, Site.type, Year, Site_number, variable, value,
                       s_value = Bromus_CN, se = b_se, fit = b_fit) %>%
  mutate(s_variable = "Bromus C:N")
f6_p <-  dplyr::select(fig6_df1, Site.type, Year, Site_number, variable, value,
                       s_value = Poa_CN, se = p_se, fit = p_fit) %>%
  mutate(s_variable = "Poa C:N")
f6_o <-  dplyr::select(fig6_df1, Site.type, Year, Site_number, variable, value,
                       s_value = Other_CN, se = o_se, fit = o_fit) %>%
  mutate(s_variable = "Other C:N")

f6_<-rbind(f6_b, f6_p, f6_o) %>% 
  na.omit() %>%
  mutate(Site.type = as.character(Site.type),
         shrub_b = ifelse(substr(Site.type,1,1) == "I", "Shrub", "Herb"))

# plot

f6 <- ggplot(f6_, aes(x = value, y = s_value, color =Site.type,
                          shape = Site.type)) +
  geom_line(aes(y=fit), show.legend = FALSE, alpha = 0.5) +
  # geom_line(aes(y=fit+se), lty = 3)+
  # geom_line(aes(y=fit-se), lty = 3)+
  geom_point() +
  geom_errorbar(aes(x = value,
                    ymin = s_value - se,
                    ymax = s_value + se,
                    color=Site.type),
                alpha = 0.5, show.legend = F) +
  scale_color_manual(name="Invasion Stage", values = lut_col) +
  scale_shape_manual(name = "Invasion Stage", values = c(15,16,17,18)) +
  facet_grid(s_variable~variable, scales="free") +
  theme_bw() +
  theme(axis.title = element_blank())+
  theme(legend.position = "bottom")+
  ggsave("figures/fig_6_plant_univariate.png", dpi = 600, width = 10, height =10)

f6_sh <- ggplot(f6_, aes(x = value, y = s_value, color =shrub_b,
                      shape = shrub_b)) +
  #geom_line(aes(y=fit), show.legend = FALSE, alpha = 0.5) +
  # geom_line(aes(y=fit+se), lty = 3)+
  # geom_line(aes(y=fit-se), lty = 3)+
  geom_point() +
  geom_errorbar(aes(x = value,
                    ymin = s_value - se,
                    ymax = s_value + se,
                    color=shrub_b),
                alpha = 0.5, show.legend = F) +
  scale_color_manual(name="Invasion Stage", values = lut_col) +
  scale_shape_manual(name = "Invasion Stage", values = c(15,16,17,18)) +
  facet_grid(s_variable~variable, scales="free") +
  theme_bw() +
  theme(axis.title = element_blank())+
  theme(legend.position = "bottom")+
  ggsave("figures/fig_6_plant_univariate_sh.png", dpi = 600, width = 10, height =10)


## c vs n fig (supplement?)==================

cn_mod<- summary(lm(soil_n_kg_ha~soil_c_kg_ha*Year, data = all_p))

ggplot(all_p, aes(x=soil_c_kg_ha, y=soil_n_kg_ha)) +
  geom_smooth(method = "lm", color = "black", show.legend = F) +
  geom_point(aes(color = Site.type, shape = Year)) +
  #facet_wrap(~Year) +
  # theme(legend.position = c(1,0.1),
  #       legend.justification = c(1,0.1)) +
  annotate("text", 
           label = paste("r2 =", round(cn_mod$r.squared,2)),
           x=15000, y=1800) +
  ggsave("figures/supfig_1_soil_n_vs_c.png", width =7, height=4)
