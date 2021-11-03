# final sem models
# setup ------------------------------------------------------------------------
set.seed(1312)
source("R/a_prep_mjcb.R")

libs<-c("lavaan",  "kableExtra")
lapply(libs, library, character.only = TRUE)
source("R/functions.R")
source("R/ggplot_sem.R")

if(file.exists("data/sem_fits.Rda")) load("data/sem_fits.Rda")
if(file.exists("data/bootstrapped_sems.Rda")) load("data/bootstrapped_sems.Rda")
#functions =====================================================================

# function for boostraping the standardised model coeficients
boot_std<-function(x) {
  ss<-standardizedsolution(x) %>%
    filter(op != "~~")
  est_stds <- ss$est.std
  names_est_stds <- paste0(ss$lhs,ss$op, ss$rhs)
  names(est_stds) <- names_est_stds
  
  return(est_stds)
}

# function for bootstrapping r2 values
boot_r2 <- function(x) {
  sm <- summary(x, rsquare=T)
  r2s<- filter(sm$PE, op=="r2") %>%
    dplyr::select(est) %>%
    pull()
  names(r2s) <-filter(sm$PE, op=="r2") %>%
    dplyr::select(lhs) %>%
    pull()
  return(r2s)
}

# data import ==================================================================

all_p_sc <- all_p %>%
  left_join(read_csv("data/climate_anomalies_summary.csv") %>%
              mutate(Year = as.character(Year), Site_number = as.character(Site_number))) %>%
  dplyr::rename(sN = soil_n_kg_ha,
                sC = soil_c_kg_ha,
                sCN = SOIL_CN,
                sMg = SOIL_Mg_kg_ha,
                sSO = SOIL_SurSo4_kg_ha,
                sCa = SOIL_Ca_kg_ha,
                OCN = Other_CN,
                BCN = Bromus_CN,
                PCN = Poa_CN,
                ltN = Litter_TN_pct,
                ltC = Litter_TC_pct,
                lCN = Litter_CN,
                aet = mean_aet,
                aez = ja_ju_aet_z,
                a2z = twoyr_aet_z,
                ae2 = twoyr_aet,
                def = ja_ju_def,
                dez = ja_ju_def_z,
                d2z = twoyr_def_z,
                de2 = twoyr_def,
                tmn = mean_tmin,
                tmz = ja_ju_tmin_z,
                p2t = ppt_2yrtot,
                p2 = ppt_2yr,
                p1 = ppt_1yr,
                sd_cwd = sd_def,
                BtC = Bromus_TC_pct,
                PsC = Poa_TC_pct,
                O_C = Other_TC_pct,
                BtN = Bromus_TN_pct,
                PsN = Poa_TN_pct,
                O_N = Other_TN_pct,
                BSC = Crypto) %>%
  mutate(NF = ANF + PNF,
         THC = NF+AIF+AIG+PNG,
         AIGr = AIG/THC,
         AIFr = AIF/THC,
         PNGr = PNG/THC,
         NFr = NF/THC,
         GtF = (AIG+PNG)/(AIF+NF))

# rescaling variables ==========================================================
for(i in 1:ncol(all_p_sc)){
  if(is.numeric(as.data.frame(all_p_sc)[,i])==TRUE){
    all_p_sc[,i] <- as.numeric(scale(all_p_sc[,i],center = F))
  }
}

# exploring correlation of climate variables ===================================
# cor(all_p_sc[,c(20:47,50:73)])
# cor(all_p[,50:60])

# #soil c and n ===========

soil_cn_model<-'
sN ~ aig_n*AIG + lit_n*Litter + cwd_n*sd_cwd + aet_n*aet +p2_n*p2+ bsc_n*BSC+ png_n*PNG + tmn_n*tmn
sC ~ AIF + aig_c*AIG +shr_c*Shrubs + cwd_c*sd_cwd +p2_c*p2 + lit_c*Litter + aet_c*aet+bsc_c*BSC +tmn_c*tmn+ png_c*PNG

#lCN ~ AIG + AIF +sd_cwd+aet+tmn+p2+Litter + PNG +BSC
AIG ~ tmn_aig*tmn + aif_aig*AIF + shr_aig*Shrubs+p2_aig*p2+ aet_aig*aet +cwd_aig*sd_cwd + png_aig*PNG
AIF ~ aig_aif*AIG+ shr_aif*Shrubs+tmn+aet
Litter ~ aig_lit*AIG + shr_lit*Shrubs + bsc_lit*BSC + sd_cwd + PNG
PNG ~ p2_png*p2 + aet_png*aet + tmn_png*tmn + cwd_png*sd_cwd + lit_png*Litter + Shrubs

# sN~~sC

## indirect effects
#C

# aet through aig -----------
i_aet_aig_c:= aet_aig * aig_c
i_aet_png_aig_c:= aet_png * png_aig * aig_c
I_aet_png_aig_c:= aig_c * ((aet_png * png_aig) + aet_aig)
# T_aet_aig_c:= (aet_aig * aig_c)+aet_c
T_aet_aig_png_c:= (aig_c * ((aet_png * png_aig) + aet_aig)) + aet_c
# tmin --------------------------
i_tmn_aig_c:= tmn_aig * aig_c
i_tmn_png_aig_c:= tmn_png * png_aig * aig_c
I_tmn_png_aig_c:= ((tmn_png * png_aig) + tmn_aig) * aig_c
# T_tmn_aig_c:= (tmn_aig * aig_c) + tmn_c
T_tmn_aig_png_c:= (((tmn_png * png_aig) + tmn_aig) * aig_c) + tmn_c
# p2 through aig -----------------
i_p2_aig_c:= p2_aig * aig_c
i_p2_png_aig_c:= p2_png*png_aig*aig_c
I_p2_png_aig_c:= ((p2_png * png_aig) + p2_aig) * aig_c
# T_p2_aig_c:= (p2_aig * aig_c)+ p2_c
T_p2_aig_png_c:= (((p2_png * png_aig) + p2_aig) * aig_c) + p2_c
# sd_cwd through aig
i_cwd_aig_c:= cwd_aig * aig_c
i_cwd_png_aig_c:= cwd_png * png_aig * aig_c
I_cwd_png_aig_c:= ((cwd_png * png_aig) + cwd_aig) * aig_c
# T_cwd_aig_c:= (cwd_aig * aig_c)+ cwd_c
T_cwd_aig_png_c:= (((cwd_png * png_aig) + cwd_aig) * aig_c) + cwd_c

# #shrubs through AIG
i_shr_aig_c:= shr_aig * aig_c
T_shr_aig_c:= (shr_aig * aig_c) + shr_c
i_shr_lit_c:= shr_lit * lit_c
I_shr_lit_aig_c:=i_shr_lit_c + i_shr_aig_c

# i_shr_bsc_lit_c:= (shr_lit + bsc_lit) * lit_c
# T_shr_lit_c:= i_shr_lit_c + shr_c
# T_shr_c:=i_shr_lit_c + shr_c + i_shr_aig_c
# c_shr_aig_vs_lit_c:= i_shr_aig_c - i_shr_lit_c
i_lit_png_aig_c:= lit_png * png_aig * aig_c
T_lit_png_aig_c:= (lit_png * png_aig * aig_c) + lit_c

# N
i_aet_aig_n:= aet_aig * aig_n
i_aet_png_aig_n:= aet_png * png_aig * aig_n
# T_aet_aig_n:= (aet_aig * aig_n)+aet_n
I_aet_png_aig_n:= aig_n * ((aet_png * png_aig) + aet_aig)
T_aet_aig_png_n:= (aig_n * ((aet_png * png_aig) + aet_aig)) + aet_n 

i_tmn_aig_n:= tmn_aig * aig_n
i_tmn_png_aig_n:= tmn_png * png_aig * aig_n
I_tmn_png_aig_n:= ((tmn_png * png_aig) + tmn_aig) * aig_n
T_tmn_aig_png_n:= (((tmn_png * png_aig) + tmn_aig) * aig_n) + tmn_n

i_cwd_aig_n:= cwd_aig * aig_n
i_cwd_png_aig_n:= cwd_png * png_aig * aig_n
I_cwd_png_aig_n:= ((cwd_png * png_aig) + cwd_aig) * aig_n
T_cwd_aig_png_n:= (((cwd_png * png_aig) + cwd_aig) * aig_n) + cwd_n
# T_cwd_aig_n:= (cwd_aig * aig_n) + cwd_n

# T_p2_aig_n:= (p2_aig * aig_n) + p2_n
i_p2_aig_n:= p2_aig * aig_n
i_p2_png_aig_n:= p2_png*png_aig*aig_n
I_p2_png_aig_n:= ((p2_png * png_aig) + p2_aig) * aig_n
T_p2_aig_png_n:= (((p2_png * png_aig) + p2_aig) * aig_n) + p2_n


i_shr_bsc_lit_c:= (shr_lit + bsc_lit) * lit_c
i_shr_bsc_lit_n:= (shr_lit + bsc_lit) * lit_n
i_shr_aig_n:= shr_aig * aig_n
# T_shr_aig_n:= (shr_aig * aig_n) + shr_n

i_shr_lit_n:= shr_lit * lit_n
I_shr_lit_aig_n:=i_shr_lit_n + i_shr_aig_n
# c_shr_aig_vs_lit_n:= i_shr_aig_n - i_shr_lit_n

i_lit_png_aig_n:= lit_png*png_aig*aig_n
i_lit_png_aig_c:= lit_png*png_aig*aig_c
T_lit_png_aig_n:= i_lit_png_aig_n + lit_n
'

scn_fit <- sem(soil_cn_model, all_p_sc)

# scn 1 and 2 ==================================================================
soil_cn_1_model <-'
sN ~  png_n*PNG + aet_n*aet +p2_n*p2 + shr_n*Shrubs + bsc_n*BSC
sC ~  png_c*PNG + aet_c*aet +p2_c*p2 + srh_c*Shrubs + lcn_c*lCN + cwd_c*sd_cwd
lCN ~ cwd_lcn*sd_cwd + aet_lcn*aet + p2_lcn*p2 + bsc_lcn*BSC
PNG ~ p2_png*p2 + aet_png*aet +  shr_png*Shrubs + nf_png*NF
NF  ~ bsc_nf*BSC + shr_nf*Shrubs + cwd_nf*sd_cwd + lcn_nf*lCN


## indirect effects

# aet through png -----------
i_aet_png_c:= aet_png * png_c
i_aet_png_n:= aet_png * png_n
T_aet_png_c:= (aet_png * png_c) + aet_c
T_aet_png_n:= (aet_png * png_n) + aet_n

# p2 through png -----------
i_p2_png_c:= p2_png * png_c
i_p2_png_n:= p2_png * png_n
T_p2_png_c:= (p2_png * png_c) + p2_c
T_p2_png_n:= (p2_png * png_n) + p2_n

# p2 through lcn
i_p2_lcn_c:= p2_lcn * lcn_c
T_p2_lcn_c:= (p2_lcn * lcn_c) + p2_c
'


scn_12_fit <- sem(soil_cn_1_model, all_p_sc %>% 
                   filter(Site.type == "Intact Sagebrush"|
                            Site.type == "Invaded Sagebrush"));scn_12_fit

summary(scn_12_fit, rsquare=T, fit.measures=T, standardized=T)
ggsem(scn_12_fit, 
      variable = "Soil Total C & N, Stages I & II",
      filename = "figures/newsems/scn_no_gps_stage1.png")
modificationindices(scn_12_fit, sort. = TRUE)[1:10,]
parameterestimates(scn_12_fit) %>%
  filter(op=="~")%>%
  group_by(lhs, op, rhs) %>%
  dplyr::summarise(max = max(pvalue, na.rm=T),
                   min=min(pvalue, na.rm=T),
                   mean=mean(pvalue, na.rm=T)) %>%
  ungroup()%>%
  arrange(desc(mean)) %>%
  filter(max >0.05, min >0.05)
resid(scn_12_fit, "cor")$cov

# soil stages 3&4 ==============================================================

soil_cn_3_model<-'
sN ~  aig_n*AIG + sd_cwd + p2_n*p2 + aif_n*AIF + aet_n*aet + tmn_n*tmn
sC ~  aig_c*AIG + sd_cwd + p2_c*p2 + aif_c*AIF + aet_c*aet
AIG ~ tmn_aig*tmn  +p2_aig*p2 + aif_aig*AIF + aet_aig*aet+ sd_cwd
AIF ~ aig_aif*AIG +  aet_aif*aet+sd_cwd

# indirect effects
i_p2_aig_c:= p2_aig * aig_c
i_p2_aig_n:= p2_aig * aig_n
T_p2_aig_c:= (p2_aig * aig_c) + p2_c
T_p2_aig_n:= (p2_aig * aig_n) + p2_n

i_tmn_aig_c:= tmn_aig * aig_c
i_tmn_aig_n:= tmn_aig * aig_n
# T_tmn_aig_c:= (tmn_aig * aig_c) + tmn_c
T_tmn_aig_n:= (tmn_aig * aig_n) + tmn_n

i_aet_aig_c:= aet_aig * aig_c
i_aet_aig_n:= aet_aig * aig_n
T_aet_aig_c:= (aet_aig * aig_c) + aet_c
T_aet_aig_n:= (aet_aig * aig_n) + aet_n

I_aet_aif_aig_n:= ((aet_aif * aif_aig) + aet_aig) * aig_n
I_aet_aif_aig_c:= ((aet_aif * aif_aig) + aet_aig) * aig_c

T_aet_aif_aig_n:= ((aet_aif * aif_aig * aig_aif) + aet_aig) * aig_n
T_aet_aif_aig_c:= ((aet_aif * aif_aig * aig_aif) + aet_aig) * aig_c

'

scn_34_fit <- sem(soil_cn_3_model, all_p_sc %>% 
                   filter(Site.type == "Cheatgrass-dominated"|
                            Site.type == "Cheatgrass Dieoff"));scn_34_fit
summary(scn_34_fit, rsquare=T, fit.measures=T, standardized=T)
ggsem(scn_34_fit,layout_df = layout_df,
      variable = "Soil Total C & N, Stages III & IV",
      filename = "figures/newsems/scn_no_gps_stage3.png")
modificationindices(scn_34_fit, sort. = TRUE)[1:10,]
parameterestimates(scn_34_fit) %>%
  filter(op=="~")%>%
  group_by(lhs, op, rhs) %>%
  dplyr::summarise(max = max(pvalue, na.rm=T),
                   min=min(pvalue, na.rm=T),
                   mean=mean(pvalue, na.rm=T)) %>%
  ungroup()%>%
  arrange(desc(mean)) %>%
  filter(max >0.05, min >0.05)
resid(scn_34_fit, "cor") 

summary(scn_3_fit, rsquare=T, fit.measures=T, standardized=T)

# 
ggsem(scn_34_fit,
      variable = "Soil Total C & N, Cheatgrass Dominated Sites",
      filename = "figures/newsems/scn_no_gps_stage3.png",
      groups = "All Sites")


# bromus ======================
bcn_model <- '
BCN ~ aig_bcn*AIG  + p2_bcn*p2 + shr_bcn*Shrubs + PNG + aif_bcn*AIF + aet_bcn*aet + tmn_bcn*tmn

AIG ~ cwd_aig*sd_cwd + aet_aig*aet  + png_aig*PNG + tmn_aig*tmn +p2_aig*p2+ shr_aig*Shrubs # + Litter
AIF ~ aet_aif*aet + cwd_aif*sd_cwd + tmn+ p2_aif*p2 +PNG + AIG+ shr_aif*Shrubs
# Litter ~ AIG +p2+ Shrubs
PNG ~ shr_png*Shrubs+AIG #+Litter

#indirect effects - climate----
i_p2_aig_bcn:= p2_aig * aig_bcn 
i_p2_aif_bcn:= p2_aif * aif_bcn
I_p2_aig_aif_bcn:= i_p2_aig_bcn + i_p2_aif_bcn
T_p2_aig_aif_bcn:= i_p2_aig_bcn + p2_bcn + i_p2_aif_bcn

i_tmn_aig_bcn:= tmn_aig*aig_bcn
T_tmn_aig_bcn:= i_tmn_aig_bcn + tmn_bcn

i_aet_aig_bcn:= aet_aig * aig_bcn
i_aet_aif_bcn:= aet_aif * aif_bcn
I_aet_aig_aif_bcn:=i_aet_aig_bcn + i_aet_aif_bcn
T_aet_aig_aif_bcn:= i_aet_aig_bcn + aet_bcn + i_aet_aif_bcn

i_cwd_aig_bcn:= cwd_aig * aig_bcn
i_cwd_aif_bcn:= cwd_aif * aif_bcn
I_cwd_aig_aif_bcn:= i_cwd_aig_bcn + i_cwd_aif_bcn

i_shr_aig_bcn:= shr_aig*aig_bcn
i_shr_png_aig_bcn:= shr_png*png_aig*aig_bcn
i_shr_aif_bcn:= shr_aif*aif_bcn
I_shr_png_aig_aif_bcn:= i_shr_aig_bcn + i_shr_aif_bcn + i_shr_png_aig_bcn
'

bcn_fit <- sem(bcn_model, all_p_sc);bcn_fit
# bromus shrub=============
bcn_12_model <- '
BCN ~ p2_bcn*p2 + shr_bcn*Shrubs +aif_bcn*AIF + aet_bcn*aet + tmn_bcn*tmn +Litter+lCN

# AIG ~ cwd_aig*sd_cwd + aet_aig*aet + tmn_aig*tmn +p2_aig*p2+ shr_aig*Shrubs # + Litter
AIF ~ cwd_aif*sd_cwd + tmn+ p2_aif*p2 +shr_aif*Shrubs+Litter+sd_cwd#aet_aif*aet + 
Litter ~ p2+ Shrubs+tmn + sd_cwd
# PNG ~ shr_png*Shrubs+AIG + Litter +sd_cwd+AIF
lCN ~ aet + tmn + Litter + p2 + Shrubs +sd_cwd# + AIG

# #indirect effects - climate----
# i_p2_aig_bcn:= p2_aig * aig_bcn 
# i_p2_aif_bcn:= p2_aif * aif_bcn
# I_p2_aig_aif_bcn:= i_p2_aig_bcn + i_p2_aif_bcn
# T_p2_aig_aif_bcn:= i_p2_aig_bcn + p2_bcn + i_p2_aif_bcn
# 
# i_tmn_aig_bcn:= tmn_aig*aig_bcn
# T_tmn_aig_bcn:= i_tmn_aig_bcn + tmn_bcn
# 
# i_aet_aig_bcn:= aet_aig * aig_bcn
# i_aet_aif_bcn:= aet_aif * aif_bcn
# I_aet_aig_aif_bcn:=i_aet_aig_bcn + i_aet_aif_bcn
# T_aet_aig_aif_bcn:= i_aet_aig_bcn + aet_bcn + i_aet_aif_bcn
# 
# i_cwd_aig_bcn:= cwd_aig * aig_bcn
# i_cwd_aif_bcn:= cwd_aif * aif_bcn
# I_cwd_aig_aif_bcn:= i_cwd_aig_bcn + i_cwd_aif_bcn
# 
# i_shr_aig_bcn:= shr_aig*aig_bcn
# i_shr_png_aig_bcn:= shr_png*png_aig*aig_bcn
# i_shr_aif_bcn:= shr_aif*aif_bcn
# I_shr_png_aig_aif_bcn:= i_shr_aig_bcn + i_shr_aif_bcn + i_shr_png_aig_bcn
'



bcn_12_fit <- sem(bcn_12_model, all_p_sc%>% 
  filter(Site.type == "Intact Sagebrush"|
           Site.type == "Invaded Sagebrush"));bcn_12_fit

summary(bcn_12_fit, rsquare=T, fit.measures=T, standardized=T)
ggsem(bcn_12_fit, layout_df = layout_df_bcn12,
      variable = "Soil Total C & N, Stages I & II",
      filename = "figures/newsems/bcn_no_gps_stage12.png")
modificationindices(bcn_12_fit, sort. = TRUE)[1:10,]
parameterestimates(bcn_12_fit) %>%
  filter(op=="~")%>%
  group_by(lhs, op, rhs) %>%
  dplyr::summarise(max = max(pvalue, na.rm=T),
                   min=min(pvalue, na.rm=T),
                   mean=mean(pvalue, na.rm=T)) %>%
  ungroup()%>%
  arrange(desc(mean)) %>%
  filter(max >0.05, min >0.05)
resid(bcn_12_fit, "cor")$cov

bcn_34_model <- '
BCN ~ p2_bcn*p2 + aif_bcn*AIF+AIG+aet+BSC+PNG+Litter
AIG ~ sd_cwd + aet_aig*aet +AIF + tmn_aig*tmn +p2_aig*p2
AIF ~ aet_aif*aet + cwd_aif*sd_cwd + tmn+ p2_aif*p2+Litter+BSC
Litter ~ AIG +p2+ BSC+sd_cwd
PNG ~ p2+tmn+BSC+AIG+Litter

# #indirect effects - climate----
# i_p2_aig_bcn:= p2_aig * aig_bcn 
# i_p2_aif_bcn:= p2_aif * aif_bcn
# I_p2_aig_aif_bcn:= i_p2_aig_bcn + i_p2_aif_bcn
# T_p2_aig_aif_bcn:= i_p2_aig_bcn + p2_bcn + i_p2_aif_bcn
# 
# i_tmn_aig_bcn:= tmn_aig*aig_bcn
# T_tmn_aig_bcn:= i_tmn_aig_bcn + tmn_bcn
# 
# i_aet_aig_bcn:= aet_aig * aig_bcn
# i_aet_aif_bcn:= aet_aif * aif_bcn
# I_aet_aig_aif_bcn:=i_aet_aig_bcn + i_aet_aif_bcn
# T_aet_aig_aif_bcn:= i_aet_aig_bcn + aet_bcn + i_aet_aif_bcn
# 
# i_cwd_aig_bcn:= cwd_aig * aig_bcn
# i_cwd_aif_bcn:= cwd_aif * aif_bcn
# I_cwd_aig_aif_bcn:= i_cwd_aig_bcn + i_cwd_aif_bcn
# 
# i_shr_aig_bcn:= shr_aig*aig_bcn
# i_shr_png_aig_bcn:= shr_png*png_aig*aig_bcn
# i_shr_aif_bcn:= shr_aif*aif_bcn
# I_shr_png_aig_aif_bcn:= i_shr_aig_bcn + i_shr_aif_bcn + i_shr_png_aig_bcn
'
layout_df_bcn34 <- read_csv("data/sem_layout.csv") %>%
  mutate(x=replace(x, metric == "sd_cwd", 0.35),
         x = replace(x, metric == "p2", -.075),
         x = replace(x, metric == "Litter", -0.5),
         y = replace(y, metric == "Litter", -0.4),
         x = replace(x, metric == "Shrubs", -0.2),
         y = replace(y, metric == "Shrubs", -0.95),
         y = replace(y, metric == "BSC", -0.95),
         x = replace(x, metric == "BSC", -0.95),
         y = replace(y, metric == "NF", 0.11),
         x = replace(x, metric == "NF", -0.9),
         y = replace(y, metric == "lCN", 0.6),
         x = replace(x, metric == "lCN", 0.95))
bcn_34_fit <- sem(bcn_34_model, all_p_sc%>% 
                    filter(Site.type == "Cheatgrass-dominated"|
                             Site.type == "Cheatgrass Dieoff"));bcn_34_fit
summary(bcn_34_fit, rsquare=T, fit.measures=T, standardized=T)

ggsem(bcn_34_fit, layout_df = layout_df_bcn34,
      variable = "Soil Total C & N, Stages I & II",
      filename = "figures/newsems/bcn_no_gps_stage34.png")
modificationindices(bcn_34_fit, sort. = TRUE)[1:10,]
parameterestimates(bcn_34_fit) %>%
  filter(op=="~")%>%
  group_by(lhs, op, rhs) %>%
  dplyr::summarise(max = max(pvalue, na.rm=T),
                   min=min(pvalue, na.rm=T),
                   mean=mean(pvalue, na.rm=T)) %>%
  ungroup()%>%
  arrange(desc(mean)) %>%
  filter(max >0.05, min >0.05)
resid(bcn_34_fit, "cor")$cov

# poa ================
pcn_model <- '
PCN ~ tmn_pcn*tmn + aet_pcn*aet + png_pcn*PNG + p2_pcn*p2 +cwd_pcn*sd_cwd

AIG ~ cwd_aig*sd_cwd + aet_aig*aet + tmn + p2_aig*p2+ shr_aig*Shrubs
Litter ~ aig_lit*AIG + shr_lit*Shrubs+p2
PNG ~ lit_png*Litter+tmn_png*tmn+p2+sd_cwd

#indirect effects - climate
i_tmn_png_pcn:= tmn_png*png_pcn
T_tmn_png_pcn:= i_tmn_png_pcn + tmn_pcn
i_lit_png_pcn:= lit_png*png_pcn

i_aet_aig_lit_png_pcn:= aet_aig*aig_lit*lit_png*png_pcn
i_p2_aig_lit_png_pcn:= p2_aig*aig_lit*lit_png*png_pcn
i_cwd_aig_lit_png_pcn:= cwd_aig*aig_lit*lit_png*png_pcn
T_aet_aig_lit_png_pcn:= i_aet_aig_lit_png_pcn + aet_pcn
T_p2_aig_lit_png_pcn:= i_p2_aig_lit_png_pcn + p2_pcn
T_cwd_aig_lit_png_pcn:= i_cwd_aig_lit_png_pcn + cwd_pcn

I_shr_aig_lit_png_pcn:= ((shr_aig*aig_lit) +shr_lit)*lit_png*png_pcn

'

pcn_fit <- sem(pcn_model, all_p_sc);pcn_fit

# poa shrubs ==================================================================
pcn_12_model <- '
PCN ~ tmn_pcn*tmn + aet_pcn*aet + png_pcn*PNG + p2_pcn*p2 +cwd_pcn*sd_cwd+lCN+AIG
lCN ~ aet + tmn + Litter + p2 + Shrubs +sd_cwd+ AIG+PNG

AIG ~ cwd_aig*sd_cwd + aet_aig*aet + tmn + p2_aig*p2+ shr_aig*Shrubs
Litter ~ aig_lit*AIG + shr_lit*Shrubs+p2+PNG+aet+sd_cwd
PNG ~ lit_png*Litter+tmn_png*tmn+p2+sd_cwd+Shrubs+aet

# #indirect effects - climate
# i_tmn_png_pcn:= tmn_png*png_pcn
# T_tmn_png_pcn:= i_tmn_png_pcn + tmn_pcn
# i_lit_png_pcn:= lit_png*png_pcn
# 
# i_aet_aig_lit_png_pcn:= aet_aig*aig_lit*lit_png*png_pcn
# i_p2_aig_lit_png_pcn:= p2_aig*aig_lit*lit_png*png_pcn
# i_cwd_aig_lit_png_pcn:= cwd_aig*aig_lit*lit_png*png_pcn
# T_aet_aig_lit_png_pcn:= i_aet_aig_lit_png_pcn + aet_pcn
# T_p2_aig_lit_png_pcn:= i_p2_aig_lit_png_pcn + p2_pcn
# T_cwd_aig_lit_png_pcn:= i_cwd_aig_lit_png_pcn + cwd_pcn
# 
# I_shr_aig_lit_png_pcn:= ((shr_aig*aig_lit) +shr_lit)*lit_png*png_pcn

'
pcn_12_fit <- sem(pcn_12_model, all_p_sc%>% 
                    filter(Site.type == "Intact Sagebrush"|
                             Site.type == "Invaded Sagebrush"));pcn_12_fit

summary(pcn_12_fit, rsquare=T, fit.measures=T, standardized=T)
ggsem(pcn_12_fit, layout_df = layout_df_bcn12,
      variable = "Soil Total C & N, Stages I & II",
      filename = "figures/newsems/pcn_no_gps_stage12.png")
modificationindices(pcn_12_fit, sort. = TRUE)[1:10,]
parameterestimates(pcn_12_fit) %>%
  filter(op=="~")%>%
  group_by(lhs, op, rhs) %>%
  dplyr::summarise(max = max(pvalue, na.rm=T),
                   min=min(pvalue, na.rm=T),
                   mean=mean(pvalue, na.rm=T)) %>%
  ungroup()%>%
  arrange(desc(mean)) %>%
  filter(max >0.05, min >0.05)
resid(pcn_12_fit, "cor")$cov
# poa grass ==================================================================
pcn_34_model <- '
PCN ~ tmn_pcn*tmn + aet_pcn*aet + png_pcn*PNG + p2_pcn*p2 +cwd_pcn*sd_cwd+Litter#+AIG#+AIF

# AIG ~ cwd_aig*sd_cwd + aet_aig*aet + tmn + p2_aig*p2
# AIF ~ cwd_aig*sd_cwd + aet_aig*aet + tmn + p2_aig*p2

Litter ~ p2+tmn+sd_cwd
PNG ~ lit_png*Litter+aet

# #indirect effects - climate
# i_tmn_png_pcn:= tmn_png*png_pcn
# T_tmn_png_pcn:= i_tmn_png_pcn + tmn_pcn
# i_lit_png_pcn:= lit_png*png_pcn
# 
# i_aet_aig_lit_png_pcn:= aet_aig*aig_lit*lit_png*png_pcn
# i_p2_aig_lit_png_pcn:= p2_aig*aig_lit*lit_png*png_pcn
# i_cwd_aig_lit_png_pcn:= cwd_aig*aig_lit*lit_png*png_pcn
# T_aet_aig_lit_png_pcn:= i_aet_aig_lit_png_pcn + aet_pcn
# T_p2_aig_lit_png_pcn:= i_p2_aig_lit_png_pcn + p2_pcn
# T_cwd_aig_lit_png_pcn:= i_cwd_aig_lit_png_pcn + cwd_pcn
# 
# I_shr_aig_lit_png_pcn:= ((shr_aig*aig_lit) +shr_lit)*lit_png*png_pcn

'
pcn_34_fit <- sem(pcn_34_model, all_p_sc%>% 
                    filter(Site.type == "Cheatgrass-dominated"|
                             Site.type == "Cheatgrass Dieoff"));pcn_34_fit

summary(pcn_34_fit, rsquare=T, fit.measures=T, standardized=T)
ggsem(pcn_34_fit, layout_df = layout_df_bcn34,
      variable = "Poa C:N, Stages III & IV",
      filename = "figures/newsems/pcn_stage34.png")
modificationindices(pcn_34_fit, sort. = TRUE)[1:10,]
parameterestimates(pcn_34_fit) %>%
  filter(op=="~")%>%
  group_by(lhs, op, rhs) %>%
  dplyr::summarise(max = max(pvalue, na.rm=T),
                   min=min(pvalue, na.rm=T),
                   mean=mean(pvalue, na.rm=T)) %>%
  ungroup()%>%
  arrange(desc(mean)) %>%
  filter(max >0.05, min >0.05)
resid(pcn_34_fit, "cor")$cov
# other plants ============

ocn_model <- '
OCN ~ aig_ocn*AIG + AIF + tmn_ocn*tmn + aet_ocn*aet  + p2_ocn*p2 +cwd_ocn*sd_cwd

AIG ~ cwd_aig*sd_cwd + aet_aig*aet + tmn_aig*tmn +p2_aig*p2+shr_aig*Shrubs+AIF
AIF ~ PNG + AIG
Litter ~ AIG + Shrubs+BSC
PNG ~ Litter+tmn+p2+sd_cwd+aet+AIG+BSC

#indirect effects - climate
i_tmn_aig_ocn:= tmn_aig*aig_ocn
T_tmn_aig_ocn:= i_tmn_aig_ocn + tmn_ocn

i_p2_aig_ocn:= p2_aig*aig_ocn
T_p2_aig_ocn:= i_p2_aig_ocn + p2_ocn

i_aet_aig_ocn:= aet_aig*aig_ocn
T_aet_aig_ocn:= i_aet_aig_ocn + aet_ocn

i_cwd_aig_ocn:= cwd_aig*aig_ocn
T_cwd_aig_ocn:= i_cwd_aig_ocn + cwd_ocn

i_shr_aig_ocn:= shr_aig*aig_ocn
'

ocn_fit <- sem(ocn_model, all_p_sc);ocn_fit

# other shrub ==================================================================
ocn_12_model <- '
OCN ~ AIF + aet_ocn*aet  + p2_ocn*p2+Shrubs
AIG ~ cwd_aig*sd_cwd + aet_aig*aet + tmn_aig*tmn +p2_aig*p2+shr_aig*Shrubs+AIF
AIF ~ Shrubs+AIG+aet
# Litter ~ AIG + Shrubs+BSC+aet+tmn
# PNG ~ tmn+p2+sd_cwd+aet+AIG+BSC+Shrubs

#indirect effects - climate
# i_tmn_aig_ocn:= tmn_aig*aig_ocn
# T_tmn_aig_ocn:= i_tmn_aig_ocn + tmn_ocn
# 
# i_p2_aig_ocn:= p2_aig*aig_ocn
# T_p2_aig_ocn:= i_p2_aig_ocn + p2_ocn
# 
# i_aet_aig_ocn:= aet_aig*aig_ocn
# T_aet_aig_ocn:= i_aet_aig_ocn + aet_ocn
# 
# i_cwd_aig_ocn:= cwd_aig*aig_ocn
# T_cwd_aig_ocn:= i_cwd_aig_ocn + cwd_ocn
# 
# i_shr_aig_ocn:= shr_aig*aig_ocn
'

ocn_12_fit <- sem(ocn_12_model, all_p_sc%>% 
                    filter(Site.type == "Intact Sagebrush"|
                             Site.type == "Invaded Sagebrush"));ocn_12_fit

summary(ocn_12_fit, rsquare=T, fit.measures=T, standardized=T)
ggsem(ocn_12_fit, layout_df = layout_df,
      variable = "Other plants C:N, Stages I & II",
      filename = "figures/newsems/ocn_stage12.png")
modificationindices(ocn_12_fit, sort. = TRUE)[1:10,]
parameterestimates(ocn_12_fit) %>%
  filter(op=="~")%>%
  group_by(lhs, op, rhs) %>%
  dplyr::summarise(max = max(pvalue, na.rm=T),
                   min=min(pvalue, na.rm=T),
                   mean=mean(pvalue, na.rm=T)) %>%
  ungroup()%>%
  arrange(desc(mean)) %>%
  filter(max >0.05, min >0.05)
resid(pcn_12_fit, "cor")$cov

# other grass ===========================
ocn_34_model <- '
OCN ~ AIG+AIF + cwd_ocn*sd_cwd

AIG ~ cwd_aig*sd_cwd + aet_aig*aet + tmn_aig*tmn +p2_aig*p2+AIF
AIF ~ AIG +aet
# Litter ~ BSC+sd_cwd+tmn
# PNG ~ Litter+p2+AIG+BSC

#indirect effects - climate
# i_tmn_aig_ocn:= tmn_aig*aig_ocn
# T_tmn_aig_ocn:= i_tmn_aig_ocn + tmn_ocn
# 
# i_p2_aig_ocn:= p2_aig*aig_ocn
# T_p2_aig_ocn:= i_p2_aig_ocn + p2_ocn
# 
# i_aet_aig_ocn:= aet_aig*aig_ocn
# T_aet_aig_ocn:= i_aet_aig_ocn + aet_ocn
# 
# i_cwd_aig_ocn:= cwd_aig*aig_ocn
# T_cwd_aig_ocn:= i_cwd_aig_ocn + cwd_ocn
# 
# i_shr_aig_ocn:= shr_aig*aig_ocn
'

ocn_34_fit <- sem(ocn_34_model, all_p_sc%>% 
                    filter(Site.type == "Cheatgrass-dominated"|
                             Site.type == "Cheatgrass Dieoff"));ocn_34_fit

summary(ocn_34_fit, rsquare=T, fit.measures=T, standardized=T)
ggsem(ocn_34_fit, layout_df = layout_df_bcn34,
      variable = "Other plants C:N, Stages III & IV",
      filename = "figures/newsems/ocn_stage12.png")
modificationindices(ocn_34_fit, sort. = TRUE)[1:10,]
parameterestimates(ocn_34_fit) %>%
  filter(op=="~")%>%
  group_by(lhs, op, rhs) %>%
  dplyr::summarise(max = max(pvalue, na.rm=T),
                   min=min(pvalue, na.rm=T),
                   mean=mean(pvalue, na.rm=T)) %>%
  ungroup()%>%
  arrange(desc(mean)) %>%
  filter(max >0.05, min >0.05)
resid(pcn_12_fit, "cor")$cov
# Litter CN =================

# better to just do an lm
l1<-lm(lCN~NF*BSC*p2, data = all_p_sc)
l2<-lm(lCN~NF*BSC +p2, data = all_p_sc)
l3<-lm(lCN~NF+BSC*p2, data = all_p_sc)
l4<-lm(lCN~NF+BSC +p2, data = all_p_sc)
AIC(l1,l2,l3,l4)
Anova(l4)
vif(l4)
summary(l4)
lcn_boot<- simpleboot::lm.boot(l4, 5000)

# saving all the fits ==========================================================
save(scn_fit, bcn_fit, ocn_fit, pcn_fit, 
     scn_12_fit,scn_34_fit,bcn_12_fit,bcn_34_fit,
     pcn_12_fit,pcn_34_fit,ocn_12_fit,ocn_34_fit,file = "data/sem_fits.Rda")

# doing the bootstrapping ======================================================

# this takes a long time per line
if(!file.exists("data/bootstrapped_sems.Rda")){
  scn_3_boot <- bootstrapLavaan(scn_34_fit,FUN = boot_std,R = 5000, verbose =T)
  scn_3r2_boot <- bootstrapLavaan(scn_34_fit,FUN = boot_r2,R = 5000, verbose =T)
  
  scn_1_boot <- bootstrapLavaan(scn_12_fit,FUN = boot_std,R = 5000, verbose =T)
  scn_1r2_boot <- bootstrapLavaan(scn_12_fit,FUN = boot_r2,R = 5000, verbose =T)

  scn_boot <- bootstrapLavaan(scn_fit,FUN = boot_std,R = 5000, verbose =T)
  scn_r2_boot <- bootstrapLavaan(scn_fit,FUN = boot_r2,R = 5000, verbose =T)
 
  bcn_boot <- bootstrapLavaan(bcn_fit,FUN = boot_std,R = 5000, verbose =T)
  bcn_r2_boot <- bootstrapLavaan(bcn_fit,FUN = boot_r2,R = 5000, verbose =T)
 
  pcn_boot <- bootstrapLavaan(pcn_fit,FUN = boot_std,R = 5000, verbose =T)
  pcn_r2_boot <- bootstrapLavaan(pcn_fit,FUN = boot_r2,R = 5000, verbose =T)
  
  ocn_boot <- bootstrapLavaan(ocn_fit,FUN = boot_std,R = 5000, verbose =T)
  ocn_r2_boot <- bootstrapLavaan(ocn_fit,FUN = boot_r2,R = 5000, verbose =T)
  
  save(scn_3_boot,  scn_1_boot,   scn_boot,    bcn_boot,    pcn_boot,    ocn_boot,
       scn_3r2_boot,scn_1r2_boot, scn_r2_boot, bcn_r2_boot, pcn_r2_boot, ocn_r2_boot,
       scn_1_boot, scn_1r2_boot, 
       file = "data/bootstrapped_sems.Rda")
}else{
load("data/bootstrapped_sems.Rda")
  }

# make a df for a table ========================================================

sem_df <- data.frame(Model = NA,
                     df = NA,
                     p = NA,
                     Chisq = NA, 
                     CFI = NA,
                     TLI = NA,
                     RMSEA = NA,
                     SRMR = NA)
semlist <- list(scn_12_fit, scn_34_fit, bcn_fit, pcn_fit,ocn_fit)
names<- c("Soil Total C & N, Stages I & II", "Soil Total C & N, Stages III & IV",
          "Bromus C:N", "Poa C:N", "Other C:N")

for(i in 1:length(semlist)){
  sem_df[i,1] <- names[i]
  sem_df[i,2] <- semlist[[i]]@Fit@test[[1]]$df
  sem_df[i,3] <- round(semlist[[i]]@Fit@test[[1]]$pvalue,2)
  sem_df[i,4] <- round(semlist[[i]]@Fit@test[[1]]$stat,2)
  sem_df[i,5] <- round(fitmeasures(semlist[[i]], "cfi") %>% as.numeric,2)
  sem_df[i,6] <- round(fitmeasures(semlist[[i]], "tli") %>% as.numeric,2)
  sem_df[i,7] <- round(fitmeasures(semlist[[i]], "rmsea") %>% as.numeric,2)
  sem_df[i,8] <- round(fitmeasures(semlist[[i]], "srmr") %>% as.numeric,2)
}
write_csv(sem_df, "data/sem_df.csv")

# plots =================================================================

layout_df <- read_csv("data/sem_layout.csv") %>%
  mutate(x=replace(x, metric == "sd_cwd", 0.35),
         x = replace(x, metric == "p2", -.075),
         x = replace(x, metric == "Litter", -0.5),
         y = replace(y, metric == "Litter", -0.4),
         x = replace(x, metric == "Shrubs", -0.2),
         y = replace(y, metric == "Shrubs", -0.95),
         y = replace(y, metric == "BSC", -0.95),
         x = replace(x, metric == "BSC", -0.95),
         y = replace(y, metric == "NF", 0.11),
         x = replace(x, metric == "NF", -0.9),
         y = replace(y, metric == "lCN", 0.2),
         x = replace(x, metric == "lCN", -0.1))

ps <- ggsem_boot(boot_fit=scn_boot,
                 fit=scn_fit, r2 =TRUE, r2_boot = scn_r2_boot,
                 layout_df=layout_df,
                 variable = "(a) Soil Total C and N, All Sites", 
                 legend=TRUE,
                 filename="figures/newsems/scn_no_gps.png")

leg<- ggsem(scn_fit, variable = "Soil Total C and N", layout_df = layout_df,
            legend=TRUE, just_leg = TRUE,
            filename="figures/newsems/legend.png")
pb <- ggsem_boot(bcn_boot,bcn_fit, variable = "(a) *B. tectorum* C:N", 
                 legend=FALSE,r2 =TRUE, r2_boot = bcn_r2_boot,
           filename="figures/newsems/bcn.png", layout_df = layout_df %>%
             mutate(x=replace(x, metric == "AIF", 0.5),
                    x=replace(x, metric == "AIG", -0.05)));pb
pp <- ggsem_boot(pcn_boot,pcn_fit, variable = "(b) *P. secunda* C:N", 
                 legend=FALSE, r2_boot = pcn_r2_boot,r2=TRUE,
           filename="figures/newsems/pcn.png", layout_df = layout_df)
po <- ggsem_boot(ocn_boot,ocn_fit, variable = "(c) Other Plants C:N", 
                 legend=FALSE, r2_boot = ocn_r2_boot,r2=TRUE,
           filename="figures/newsems/ocn.png", layout_df = layout_df)
ggarrange(ggarrange(pb, pp, nrow=1,ncol=2),
          ggarrange(po, leg, nrow=1, ncol=2),
          nrow=2)+
  ggsave("figures/newsems/plant_cn_square.png", height = 14, width=14)

ggarrange(pb, pp, po, leg, nrow=1, ncol=4, widths=c(1,1,1,.3)) +
  ggsave("figures/newsems/plant_cn_wide.png", height=8, width=23)



ps3 <- ggsem_boot(boot_fit=scn_3_boot, 
                  layout_df = layout_df%>%
                    mutate(y=replace(y, metric == "AIF", -0.955)),
                  fit=scn_34_fit, r2 =TRUE, r2_boot = scn_3r2_boot,
                  variable = "(b) Soil Total N & C: Stages III & IV", 
                  legend=FALSE,
                  filename="figures/newsems/scn_3_no_gps.png");ps3


layout_df_ps1 <- layout_df %>%
  mutate(x = replace(x, metric == "PNG", -0.1),
         y = replace(y, metric == "PNG", -0.25),
         x = replace(x, metric == "lCN", -0.35),
         y = replace(y, metric == "lCN", 0.35),
         x = replace(x, metric == "aet", 0.3),
         y = replace(y, metric == "aet", -0.95))
ps1 <- ggsem_boot(boot_fit=scn_1_boot,
                  fit=scn_12_fit, r2 =TRUE, r2_boot = scn_1r2_boot,
                  layout_df = layout_df_ps1,
                  variable = "(a) Soil Total N & C: Stages I & II", 
                  legend=FALSE,
                  filename="figures/newsems/scn_1_no_aig.png");ps1

soil_panel<- ggarrange(ps1, ps3, leg, nrow=1, 
          widths = c(1,1, 0.5))

ggsave("figures/newsems/soilcn_panel.png", soil_panel, height=8, width = 20, 
       dpi=600, bg="white")
