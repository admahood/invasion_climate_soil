# creating sem tables
library(tidyverse)

lut_et<- c("c"="Contrast", "i"= "Indirect (single pathway)","d" = "Direct",
           "I"= "Indirect (multiple pathways)","T"="Total")

# scn_boot<-readRDS("data/scn_boot.RDS")
load("data/bootstrapped_sems.Rda")

# Tables
# soil C and N, stages 1 & 2 ===================================================
scn_1_cis <- scn_1_boot %>%
  as.data.frame() %>%
  pivot_longer(cols = names(.),names_to = "var", values_to = "vals") %>%
  group_by(var) %>%
  dplyr::summarise(ci_025 = HDInterval::hdi(vals)%>% pluck(1),
                   ci_975 = HDInterval::hdi(vals) %>% pluck(2),
                   median = quantile(vals, probs = 0.5)) %>%
  ungroup() %>%
  mutate(type = ifelse(str_detect(var, "~"), "regression", "user"),
         sig = ifelse(ci_025 *ci_975 > 0, "*", ""))

# whatever direct effects are needed in the table we can add if necessary (pre-bootstrap)
scn_c_1_tab <- scn_1_cis %>%
  filter(type == "user") %>%
  mutate(effect_type = lut_et[str_sub(var, 1,1)],
         var = str_extract(var,"^(.*?):"),
         var = str_sub(var, 3, nchar(var)-1),
         exogenous_var = str_extract(var,"^(.*?)_")%>%
           str_sub(1,nchar(.)-1),
         resp = str_sub(var, nchar(var), nchar(var)),
         mediator = str_match(var, "_(.*?)_c")[,2] %>% str_replace("_"," ")) %>%
  filter(effect_type != "Contrast") %>%
  filter(resp == "c") %>%
  dplyr::select(exogenous_var, effect_type, "mediator(s)" = mediator, median, ci_025, ci_975, sig) %>%
  arrange(exogenous_var,effect_type);scn_c_1_tab

scn_c_1_tab %>%
  write_csv("data/booted_sc_1_table.csv")

scn_n_1_tab <- scn_1_cis %>%
  filter(type == "user") %>%
  mutate(effect_type = lut_et[str_sub(var, 1,1)],
         var = str_extract(var,"^(.*?):"),
         var = str_sub(var, 3, nchar(var)-1),
         exogenous_var = str_extract(var,"^(.*?)_")%>%
           str_sub(1,nchar(.)-1),
         resp = str_sub(var, nchar(var), nchar(var)),
         mediator = str_match(var, "_(.*?)_n")[,2] %>% str_replace("_"," ")) %>%
  filter(effect_type != "Contrast") %>%
  filter(resp == "n") %>%
  dplyr::select(exogenous_var, effect_type, "mediator(s)" = mediator, median, ci_025, ci_975, sig) %>%
  arrange(exogenous_var,effect_type);scn_n_1_tab

scn_n_1_tab %>%
  write_csv("data/booted_sn_1_table.csv")

# soil C and N, stages 3 & 4 ===================================================
scn_3_cis <- scn_3_boot %>%
  as.data.frame() %>%
  pivot_longer(cols = names(.),names_to = "var", values_to = "vals") %>%
  group_by(var) %>%
  dplyr::summarise(ci_025 = HDInterval::hdi(vals)%>% pluck(1),
                   ci_975 = HDInterval::hdi(vals) %>% pluck(2),
                   median = quantile(vals, probs = 0.5)) %>%
  ungroup() %>%
  mutate(type = ifelse(str_detect(var, "~"), "regression", "user"),
         sig = ifelse(ci_025 *ci_975 > 0, "*", ""))

# whatever direct effects are needed in the table we can add if necessary (pre-bootstrap)
scn_c_3_tab <- scn_3_cis %>%
  filter(type == "user") %>%
  mutate(effect_type = lut_et[str_sub(var, 1,1)],
         var = str_extract(var,"^(.*?):"),
         var = str_sub(var, 3, nchar(var)-1),
         exogenous_var = str_extract(var,"^(.*?)_")%>%
           str_sub(1,nchar(.)-1),
         resp = str_sub(var, nchar(var), nchar(var)),
         mediator = str_match(var, "_(.*?)_c")[,2] %>% str_replace("_"," ")) %>%
  filter(effect_type != "Contrast") %>%
  filter(resp == "c") %>%
  dplyr::select(exogenous_var, effect_type, 
                "mediator(s)" = mediator, median, ci_025, ci_975, sig) %>%
  arrange(exogenous_var,effect_type);scn_c_3_tab

scn_c_3_tab %>%
  write_csv("data/booted_sc_3_table.csv")

scn_n_3_tab <- scn_3_cis %>%
  filter(type == "user") %>%
  mutate(effect_type = lut_et[str_sub(var, 1,1)],
         var = str_extract(var,"^(.*?):"),
         var = str_sub(var, 3, nchar(var)-1),
         exogenous_var = str_extract(var,"^(.*?)_")%>%
           str_sub(1,nchar(.)-1),
         resp = str_sub(var, nchar(var), nchar(var)),
         mediator = str_match(var, "_(.*?)_n")[,2] %>% str_replace("_"," ")) %>%
  filter(effect_type != "Contrast") %>%
  filter(resp == "n") %>%
  dplyr::select(exogenous_var, effect_type, "mediator(s)" = mediator, median, ci_025, ci_975, sig) %>%
  arrange(exogenous_var,effect_type);scn_n_3_tab

scn_n_3_tab %>%
  write_csv("data/booted_sn_3_table.csv")
# bromus c:n ===================================================================
bcn_cis <- bcn_boot %>%
  as.data.frame() %>%
  pivot_longer(cols = names(.),names_to = "var", values_to = "vals") %>%
  group_by(var) %>%
  dplyr::summarise(ci_025 = quantile(vals, probs = 0.025),
            ci_975 = quantile(vals, probs = 0.975),
            median = quantile(vals, probs = 0.5)) %>%
  ungroup() %>%
  mutate(type = ifelse(str_detect(var, "~"), "regression", "user"),
         sig = ifelse(ci_025 *ci_975 > 0, "*", ""))

bcn_tab<-bcn_cis %>%
  filter(type == "user") %>%
  mutate(effect_type = lut_et[str_sub(var, 1,1)],
         var = str_extract(var,"^(.*?):"),
         var = str_sub(var, 3, nchar(var)-1),
         exogenous_var = str_extract(var,"^(.*?)_")%>%
           str_sub(1,nchar(.)-1),
         mediator = str_match(var, "_(.*?)_bcn")[,2] %>% str_replace_all("_"," ")) %>%
  filter(effect_type != "Contrast") %>%
  dplyr::select(exogenous_var, effect_type, "mediator(s)" = mediator, median, ci_025, ci_975, sig) %>%
  arrange(exogenous_var,effect_type);bcn_tab

bcn_tab%>%write_csv("data/booted_bcn_table.csv")

# poa c:n ======================================================================
pcn_cis <- pcn_boot %>%
  as.data.frame() %>%
  pivot_longer(cols = names(.),names_to = "var", values_to = "vals") %>%
  group_by(var) %>%
  dplyr::summarise(ci_025 = quantile(vals, probs = 0.025),
            ci_975 = quantile(vals, probs = 0.975),
            median = quantile(vals, probs = 0.5)) %>%
  ungroup() %>%
  mutate(type = ifelse(str_detect(var, "~"), "regression", "user"),
         sig = ifelse(ci_025 *ci_975 > 0, "*", ""))

pcn_tab<- pcn_cis%>%
  filter(type == "user") %>%
  mutate(effect_type = lut_et[str_sub(var, 1,1)],
         var = str_extract(var,"^(.*?):"),
         var = str_sub(var, 3, nchar(var)-1),
         exogenous_var = str_extract(var,"^(.*?)_")%>%
           str_sub(1,nchar(.)-1),
         mediator = str_match(var, "_(.*?)_pcn")[,2] %>% str_replace_all("_"," ")) %>%
  filter(effect_type != "Contrast") %>%
  dplyr::select(exogenous_var, effect_type, "mediator(s)" = mediator, median, ci_025, ci_975, sig) %>%
  arrange(exogenous_var,effect_type);pcn_tab

pcn_tab%>%write_csv("data/booted_pcn_table.csv")

# other c:n ====================================================================
ocn_cis <- ocn_boot %>%
  as.data.frame() %>%
  pivot_longer(cols = names(.),names_to = "var", values_to = "vals") %>%
  group_by(var) %>%
  dplyr::summarise(ci_025 = quantile(vals, probs = 0.025),
            ci_975 = quantile(vals, probs = 0.975),
            median = quantile(vals, probs = 0.5)) %>%
  ungroup() %>%
  mutate(type = ifelse(str_detect(var, "~"), "regression", "user"),
         sig = ifelse(ci_025 *ci_975 > 0, "*", ""))

ocn_tab<-ocn_cis %>%
  filter(type == "user") %>%
  mutate(effect_type = lut_et[str_sub(var, 1,1)],
         var = str_extract(var,"^(.*?):"),
         var = str_sub(var, 3, nchar(var)-1),
         exogenous_var = str_extract(var,"^(.*?)_")%>%
           str_sub(1,nchar(.)-1),
         mediator = str_match(var, "_(.*?)_ocn")[,2] %>% str_replace_all("_"," ")) %>%
  filter(effect_type != "Contrast") %>%
  dplyr::select(exogenous_var, effect_type, "mediator(s)" = mediator, median, ci_025, ci_975, sig) %>%
  arrange(exogenous_var,effect_type);ocn_tab

ocn_tab%>%write_csv("data/booted_ocn_table.csv")



# figures ======================================================================
# scn_boot %>%
#   as.data.frame() %>%
#   pivot_longer(cols = names(.),names_to = "var", values_to = "vals") %>%
#   mutate(type = ifelse(str_detect(var, "~"), "regression", "user"))%>%
#   filter(type == "user")%>%
#   mutate(var = str_extract(var,"^(.*?):"))%>%
#   ggplot(aes(x=vals)) +
#   geom_vline(xintercept = 0, lty=2)+
#   geom_density() +
#   facet_wrap(~var, scales="free") +
#   theme(axis.title = element_blank(),
#         axis.ticks = element_blank(),
#         axis.text = element_blank()) 

# covariance matrices =====================
load("data/sem_fits.Rda")
library(lavaan)
# fitted(scn_fit)$cov %>%
#   rbind(fitted(scn_fit)$mean)%>%
#   round(3) %>%
#   as_tibble(rownames = "x") %>%
#   mutate(x = replace(x, x == "", "mean")) %>%
#   write_csv("data/cm_scn.csv")

fitted(scn_12_fit)$cov %>%
  rbind(fitted(scn_12_fit)$mean)%>%
  round(3) %>%
  as_tibble(rownames = "x") %>%
  mutate(x = replace(x, x == "", "mean")) %>%
  write_csv("data/cm_scn12.csv")

fitted(scn_34_fit)$cov %>%
  rbind(fitted(scn_34_fit)$mean)%>%
  round(3) %>%
  as_tibble(rownames = "x") %>%
  mutate(x = replace(x, x == "", "mean")) %>%
  write_csv("data/cm_scn34.csv")

fitted(bcn_fit)$cov %>%
  rbind(fitted(bcn_fit)$mean)%>%
  round(3) %>%
  as_tibble(rownames = "x") %>%
  mutate(x = replace(x, x == "", "mean"))%>%
  write_csv("data/cm_bcn.csv")

fitted(ocn_fit)$cov %>%
  rbind(fitted(ocn_fit)$mean)%>%
  round(3) %>%
  as_tibble(rownames = "x") %>%
  mutate(x = replace(x, x == "", "mean"))%>%
  write_csv("data/cm_ocn.csv")

fitted(pcn_fit)$cov %>%
  rbind(fitted(pcn_fit)$mean)%>%
  round(3) %>%
  as_tibble(rownames = "x") %>%
  mutate(x = replace(x, x == "", "mean"))%>%
  write_csv("data/cm_pcn.csv")

