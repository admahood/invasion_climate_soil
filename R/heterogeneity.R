source("R/a_prep_mjcb.R")

all_sd <- all_sd %>%
  mutate(`Invasion Stage` = factor(Site.type, 
                                   levels = c("Intact Sagebrush",
                                              "Invaded Sagebrush",
                                              "Cheatgrass-dominated",
                                              "Cheatgrass Dieoff")))

heterogeneity_table <- data.frame(Variable = NA, "Intact_Sagebrush"=NA, 
                                  "Invaded_Sagebrush" = NA, 
                                  "Cheatgrass"=NA, "Dieoff"=NA, 
                                  "Shrub" =NA, "Herbaceous" =NA, 
                                  "y2013"=NA, "y2016"=NA)
het <- all_sd%>%
  dplyr::select(Year,Litter_TN_pct,Litter_TC_pct,Litter_CN,SOIL_SurSo4_kg_ha,
                SOIL_Ca_kg_ha,SOIL_Mg_kg_ha,SOIL_CN,soil_n_kg_ha,soil_c_kg_ha,
                Site.type,shrub_b)
for (i in 1:9){
  heterogeneity_table[i,1] <- names(het)[i+1]
  is <- agricolae::kruskal(het[,i+1], 
                           het$Site.type, 
                           p.adj="bonferroni")$groups %>%
    as.data.frame()%>%
    tibble::rownames_to_column("trt") %>%
    left_join(
      agricolae::kruskal(het[,i+1], 
                         het$Site.type, 
                         p.adj="bonferroni")$means %>%
        as.data.frame()%>%
        tibble::rownames_to_column("trt"))
  is<-dplyr::rename(is,mean = `het...i...1.`)
  
  if(length(unique(is$groups)) == 1 ) is$groups <- "   "
  
  heterogeneity_table[i,2] <- paste(signif(is[is$trt =="Intact Sagebrush",]$mean, 1),
                                    is[is$trt =="Intact Sagebrush",]$groups %>% str_pad(3,"right"))
  heterogeneity_table[i,3] <- paste(signif(is[is$trt =="Invaded Sagebrush",]$mean, 1),
                                    is[is$trt =="Invaded Sagebrush",]$groups%>% str_pad(3,"right"))
  heterogeneity_table[i,4] <- paste(signif(is[is$trt =="Cheatgrass-dominated",]$mean, 1),
                                    is[is$trt =="Cheatgrass-dominated",]$groups%>% str_pad(3,"right"))
  heterogeneity_table[i,5] <- paste(signif(is[is$trt =="Cheatgrass Dieoff",]$mean, 1),
                                    is[is$trt =="Cheatgrass Dieoff",]$groups%>% str_pad(3,"right"))
  year <- agricolae::kruskal(het[,i+1], 
                             het$Year, 
                             p.adj="bonferroni")$means %>%
    as.data.frame()%>%
    tibble::rownames_to_column("trt")
  year<-dplyr::rename(year,mean = `het...i...1.`)
  
  phys <- agricolae::kruskal(het[,i+1], 
                             het$shrub_b, 
                             p.adj="bonferroni")$means %>%
    as.data.frame()%>%
    tibble::rownames_to_column("trt")%>%
    left_join(agricolae::kruskal(het[,i+1], 
                                 het$shrub_b, 
                                 p.adj="bonferroni")$groups %>%
                as.data.frame()%>%
                tibble::rownames_to_column("trt"))
  phys<-dplyr::rename(phys,mean = `het...i...1.`)
  
  if(length(unique(phys$groups)) == 1 ){phys$groups <- "   "}
  
  heterogeneity_table[i,6] <- paste(signif(phys[phys$trt =="Shrub",]$mean, 1),
                                    phys[phys$trt =="Shrub",]$groups%>% str_pad(3,"right"))
  heterogeneity_table[i,7] <- paste(signif(phys[phys$trt =="Grass",]$mean, 1),
                                    phys[phys$trt =="Grass",]$groups%>% str_pad(3,"right"))
  
  year <- agricolae::kruskal(het[,i+1], 
                             het$Year, 
                             p.adj="bonferroni")$means %>%
    as.data.frame()%>%
    tibble::rownames_to_column("trt")%>%
    left_join(agricolae::kruskal(het[,i+1], 
                                 het$Year, 
                                 p.adj="bonferroni")$groups %>%
                as.data.frame()%>%
                tibble::rownames_to_column("trt"))
  year<-dplyr::rename(year,mean = `het...i...1.`)
  
  if(length(unique(year$groups)) == 1 ){year$groups <- "   "}
  
  
  heterogeneity_table[i,8] <- paste(signif(year[year$trt =="2013",]$mean, 1),
                                    year[year$trt =="2013",]$groups%>% str_pad(3,"right"))
  heterogeneity_table[i,9] <- paste(signif(year[year$trt =="2016",]$mean, 1),
                                    year[year$trt =="2016",]$groups%>% str_pad(3,"right"))
}

write_csv(heterogeneity_table, "figures/heterogeneity_table.csv")
