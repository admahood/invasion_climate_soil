# tidymodel syntax approach to diffs
source("R/a_prep_mjcb.R")
library(lme4)
library(ggtext)
library(emmeans)
library(tidymodels)
library(lmerTest)
resp <- c( "Bromus_TC_pct", "Bromus_CN","Bromus_TN_pct",  
           "Poa_TC_pct", "Poa_CN", "Poa_TN_pct",
          "Other_TC_pct",  "Other_CN", "Other_TN_pct",    
          "soil_n_kg_ha", "soil_c_kg_ha","SOIL_CN",
          "SOIL_SurSo4_kg_ha", "SOIL_Ca_kg_ha", "SOIL_Mg_kg_ha" ,
          "Litter_TN_pct", "Litter_TC_pct", "Litter_CN", 
            "total_mineral_n",
          "NO3_kg_ha","NH4_kg_ha")
lut_inv<- c("Intact Sagebrush" = "I",
            "Invaded Sagebrush" = "II",
            "Cheatgrass-dominated" = "III",
            "Cheatgrass Dieoff" = "IV")
plot_labs <- c("(a)","(b)","(c)","(d)","(e)","(f)","(g)","(h)","(i)")

models <- 
  all %>%
  dplyr::select(site_number = Site_number, 
                invasion_stage = Site.type,
                transect = Transect,
                year = Year,
                all_of(resp)) %>%
  mutate(invasion_stage = lut_inv[invasion_stage]) %>%
  pivot_longer(cols = resp, names_to = "response", values_to = "value") %>%
  na.omit() %>%
  group_nest(response) %>%
  mutate(model = map(data, ~lmer(value ~ invasion_stage*year + (1|site_number/transect), 
                                 data=.x)))

post_mod <- 
  models %>%
  mutate(aov = map(model, anova, ddf="Kenward-Roger"),
         stage_by_year = map(model, emmeans, contr = "tukey",
                             specs = ~invasion_stage|year),
         stage = map(model, emmeans, contr = "tukey",
                             specs = ~invasion_stage),
         year_by_stage = map(model, emmeans, contr = "tukey",
                             specs = ~year|invasion_stage)) 

plot_the_stuff<- function(x,  post_mod) {
  gc()
  
  testmod<-post_mod %>%
    filter(response == x) 
  dat<- pluck(testmod$data,1)%>%
    mutate(response = testmod$response)
  aov <- pluck(testmod$aov,1)
  int_p <- aov$`Pr(>F)`[3]
  inv_p <- aov$`Pr(>F)`[1]
  
  # filtering out some "outliers" for better visibility (they're not excluded from 
  # calculating the confidence intervals, p-values, etc.. just a visualization aid)
  if(x=="Poa_TC_pct") dat<- filter(dat, value>35)
  if(x =="SOIL_SurSo4_kg_ha") dat<- filter(dat, value<20)
  if(x == "Other_TC_pct") dat <- filter(dat, value>30)
  
  nudge <- (range(dat$value, na.rm = T) %>% diff)*0.1

  
 
    emmeans_y<-   pluck( testmod$year_by_stage, 1)$emmeans%>%
      multcomp::cld(Letters=letters, decreasing=FALSE) %>% 
      as.data.frame %>%
      mutate(.group = trimws(.group))
    
    ucls<- emmeans_y%>%
      group_by(invasion_stage) %>%
      dplyr::summarise(UCL = max(upper.CL)) %>%
      ungroup
    
    emmeans_s<-
      pluck(testmod$stage, 1)$emmeans%>%
      multcomp::cld(Letters=letters, decreasing=FALSE) %>%
      as.data.frame %>%
      left_join(ucls) %>%
      mutate(.group = str_to_upper(.group) %>% trimws)
    
     p<- ggplot(dat %>%
                      mutate(year = ifelse(year == "2016", "'16", "'13"))) +
        geom_jitter(aes(x=year, y=value, color = invasion_stage), 
                    alpha=.5, width = 0.2, shape=1) +
        geom_errorbar(data=emmeans_y%>%
                        mutate(year = ifelse(year == "2016", "'16", "'13")),
                      width = 0.4, aes(x=year,
                      ymin=lower.CL, 
                      ymax=upper.CL,
                      color=invasion_stage))+
        geom_point(data = emmeans_y%>%
                     mutate(year = ifelse(year == "2016", "'16", "'13")), 
                   aes(x=year, y=emmean, color = invasion_stage),
                   size=3, shape=18)+
        
        scale_color_manual(name = "Invasion Stage", 
                           values=site_colors)+
        scale_y_continuous(labels = scales::label_number())+
        ylab(lut_variables[dat$response])+
        xlab("Year") +
        facet_wrap(~invasion_stage, nrow=1)+
        ggthemes::theme_clean() +
        theme(axis.title.y = element_markdown(),
              legend.position = "none",
              panel.spacing = unit(0, "lines"),
              panel.border = element_rect(fill=NA, size=0))
     if(int_p < 0.05) {p<-p +
       geom_label(data = emmeans_y%>%
                  mutate(year = ifelse(year == "2016", "'16", "'13")),
                nudge_y =nudge,
                size=4,fontface="bold",
                label.size=unit(0.00, "lines"),
                label.padding = unit(0.01,"lines"),
                aes(label=emmeans_y$.group,color=invasion_stage,
                    x=year, y=upper.CL))
     }
     if(int_p > 0.05 & inv_p < 0.05){
       p<- p+ geom_label(data = emmeans_s,
                  nudge_y =nudge,alpha=0.35,
                  size=4,fontface="bold", color = "grey40",
                  label.size=unit(0.00, "lines"),
                  label.padding = unit(0.01,"lines"),
                  aes(label=emmeans_s$.group,
                      x=1.5, y=UCL))
     }
                    
     return(p)
    }


#test plots
# plot_the_stuff(x="soil_n_kg_ha",  post_mod=post_mod)
# plot_the_stuff(x="Litter_TC_pct",  post_mod=post_mod)

plant_vars_y <- resp[1:9] %>%
  lapply(plot_the_stuff, post_mod = post_mod)%>%
  ggarrange(plotlist = ., nrow=3, ncol=3, labels = plot_labs) 

ggsave(plant_vars_y, filename = "figures/tidy_plantvars_year.png",height=8, width=10)

soil_vars_y <- resp[c(10:12, 16:21)] %>%
  lapply(plot_the_stuff, post_mod = post_mod)%>%
  ggarrange(plotlist = ., nrow=3, ncol=3, labels = plot_labs)

ggsave(soil_vars_y, filename="figures/tidy_soilvars_year.png",height=8, width=10)
