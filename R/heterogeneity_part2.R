source("R/a_prep_mjcb.R")

all_sd <- all_sd %>%
  mutate(`Invasion Stage` = factor(Site.type, 
                                   levels = c("Intact Sagebrush",
                                              "Invaded Sagebrush",
                                              "Cheatgrass-dominated",
                                              "Cheatgrass Dieoff")))

heterogeneity_table_long <- data.frame(Variable = NA,
                                       mean = NA,
                                       group = NA)
het <- all_sd%>%
  dplyr::select(Year,Litter_TN_pct,Litter_TC_pct,Litter_CN,SOIL_SurSo4_kg_ha,
                SOIL_Ca_kg_ha,SOIL_Mg_kg_ha,SOIL_CN,soil_n_kg_ha,soil_c_kg_ha,
                Site.type,shrub_b)
res<- list()
for (i in 1:9){
  # heterogeneity_table[i,1] <- names(het)[i+1]
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
  res[[i]]<-dplyr::select(is, mean = `het...i...1.`, groups, trt) %>%
    mutate(variable = names(het)[i+1])
  if(length(unique(res[[i]]$groups))==1) res[[i]]$groups<-""
}
heterogeneity_table_long <- bind_rows(res) %>%
  pivot_wider(id_cols = variable, names_from = "trt", values_from = c("mean", "groups")) %>%
  dplyr::select(variable,
                "mean_Intact Sagebrush", "groups_Intact Sagebrush",
                "mean_Invaded Sagebrush","groups_Invaded Sagebrush",
                "mean_Cheatgrass-dominated","groups_Cheatgrass-dominated",
                "mean_Cheatgrass Dieoff" ,"groups_Cheatgrass Dieoff")

write_csv(heterogeneity_table_long, "figures/heterogeneity_table_long.csv")
