# graphical abstract fig

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


leg<- ggsem_boot_ga(boot_fit = scn_boot,fit=scn_fit, variable = "Soil Total C and N", layout_df = layout_df,
            legend=TRUE, just_leg = TRUE)

# leg_df<- data.frame(vals = c("< -0.75",
#                              "-0.5 - -0.75",
#                              "0 - -0.5",
#                              "0 - 0.5",
#                              "0.5 - 0.75",
#                              "> 0.75"),
#                     cols = c("red",
#                       "orange",
#                       "peachpuff",
#                       "grey80",
#                      "grey40",
#                      "black"),
#       lty=c(6,6,6,1,1,1))
# lut_cols<-c("< -0.75" = "red",
#             "-0.5 - -0.75" = "orange",
#             "0 - -0.5" = "peachpuff",
#             "0 - 0.5" = "grey80",
#             "0.5 - 0.75" = "grey40",
#             "> 0.75"= "black")
# lut_lty<-c("< -0.75" = 6,
#            "-0.5 - -0.75" = 6,
#            "0 - -0.5" = 6,
#            "0 - 0.5" = 1,
#            "0.5 - 0.75" = 1,
#            "> 0.75"= 1)
# leg <- get_legend(
# 
#   ggplot(leg_df, aes(x=lty, y=lty)) +
#     geom_line(aes(color = vals, lty=vals), lwd=3) +
#     theme_classic()+
#     theme(legend.key.size = unit(1.5,"cm"),
#           legend.text = element_text(size=15),
#           legend.title = element_text(size=15),
#           legend.background = element_rect(fill="transparent"))+
#     guides(color=guide_legend(ncol=1))+
#     scale_linetype_manual(values = lut_lty,
#                           name = "Standardized\nEstimates")+
#     scale_color_manual(values = lut_cols,
#                        name = "Standardized\nEstimates")
# )


ps3 <- ggsem_boot_ga(boot_fit=scn_3_boot, 
                  layout_df = layout_df%>%
                    mutate(y=replace(y, metric == "AIF", -0.955)),
                  fit=scn_34_fit, r2 =FALSE, r2_boot = scn_3r2_boot,
                  variable = "(b) Burned Sites, Annual-Dominated", 
                  legend=FALSE,
                  filename="figures/newsems/scn_3_graphical_abstract.png")


layout_df_ps1 <- layout_df %>%
  mutate(x = replace(x, metric == "PNG", -0.1),
         y = replace(y, metric == "PNG", -0.25),
         x = replace(x, metric == "lCN", -0.25),
         y = replace(y, metric == "lCN", 0.35))

ps1 <- ggsem_boot_ga(boot_fit=scn_1_boot,
                  fit=scn_12_fit, r2 =FALSE, r2_boot = scn_1r2_boot,
                  layout_df = layout_df_ps1,
                  variable = "(a) Unburned Sites, Shrubs Intact", 
                  legend=FALSE,
                  filename="figures/newsems/scn_1_graphical_abstract.png");ps1

ga_plot<- ggarrange(ps1, ps3, leg, nrow=1, 
          widths = c(1,1, 0.5))
ggsave(ga_plot, filename="figures/newsems/soilcn_graphical_abstract.png", height=8, width = 20, bg="white")

ggsave(cowplot::ggdraw(xlim = c(0,18), ylim = c(0,8)) +
  draw_plot(ps1, height= 8, width = 9, x = 0, y=0) +
  draw_plot(ps3, height = 8, width = 9, x=9, y=0) +
  draw_plot(leg, x=9.5, y=2 ),
  filename = "figures/newsems/graphical_abstract_cowplot.png", height = 8, width = 16, bg="white")
# Graphical abstract text
# 
# The sensitivity of herbaceous cover to inter-annual climatic variability mediates the effects 
# of invasion of exotic annual grasses and changes in plant functional group composition 
# on soil total C and N

# Different variables determined soil total C and soil total N in different invasion stages.
# Climate-driven differences in herbaceous plant cover drove soil total C and soil 
# total N in different directions depending on the composition of the herbaceous 
# vegetation, with annual invaders resulting in losses of soil total C and soil total N, 
# and native perennials resulting in increases.
