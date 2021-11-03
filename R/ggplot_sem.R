# setup ========================================================================
libs <- c("tidyverse", "tidygraph","cowplot", "ggraph", "lavaan", "ggpubr", 
          "classInt", "ggtext","HDInterval")
# install.packages("ggraph")
lapply(libs, library, character.only = TRUE)

# if there's a font error:
# extrafont::font_import()
# then restart rstudio
# extrafont::loadfonts()

# plotting functions adapted, heavily modified, from: ==========================
# https://drsimonj.svbtle.com/ggsem-plot-sem-models-with-ggplot2 
# also this:
# https://cran.r-project.org/web/packages/ggraph/vignettes/Edges.html

ggsem <- function(fit, filename, variable,layout_df = layout_df, alpha = 0.05,
                  legend=TRUE, just_leg = FALSE) {

     # Extract standardized parameters
  params <- lavaan::standardizedSolution(fit) 
  # Edge properties
  
  param_edges <- params %>% 
    filter(op %in% c("=~", "~", "~~"), lhs != rhs) %>% #, pvalue < .10) %>%
    mutate(sig = replace(pvalue, pvalue > alpha, 2)) %>%
    mutate(sig = replace(sig, pvalue < alpha, 1)) %>%
    filter(sig != 2) %>%
    transmute(to = lhs,
              from = rhs,
              sig=sig,
              val = est.std,
              type = dplyr::case_when(
              op == "~"  ~ "regression",
              op == "~~" ~ "correlation",
              TRUE ~ NA_character_))
  
  lut_cols<-c("< -0.75" = "red",
              "-0.5 - -0.75" = "orange",
              "0 - -0.5" = "peachpuff", 
              "0 - 0.5" = "grey80", 
              "0.5 - 0.75" = "grey40", 
              "> 0.75"= "black")
  lut_lty<-c("< -0.75" = 6,
             "-0.5 - -0.75" = 6,
             "0 - -0.5" = 6, 
             "0 - 0.5" = 1, 
             "0.5 - 0.75" = 1, 
             "> 0.75"= 1)
  
  param_edges <- param_edges %>%
    mutate(class = val)%>%
    mutate(class = replace(class, val >= 0.75, "> 0.75")) %>%
    mutate(class = replace(class, val >= 0.5 & val< 0.75, "0.5 - 0.75")) %>%
    mutate(class = replace(class, val >= 0 & val< 0.5, "0 - 0.5")) %>%
    mutate(class = replace(class, val >= -0.5 & val< 0, "0 - -0.5")) %>%
    mutate(class = replace(class, val >= -0.75 & val < -0.5, "-0.5 - -0.75")) %>%
    mutate(class = replace(class, val <= -0.75, "< -0.75")) %>%
    mutate(class = factor(class, levels = c("> 0.75",
                                            "0.5 - 0.75",
                                            "0 - 0.5",
                                            "0 - -0.5",
                                            "-0.5 - -0.75",
                                            "< -0.75"))) %>%
    mutate(sign = ifelse(val>0, 1,2))
  
  # Node properties
  param_nodes <- params %>% 
    filter(lhs == rhs) %>% 
    transmute(metric = lhs, e = est.std)
  
  # Complete Graph Object
  param_graph1 <- tidygraph::tbl_graph(param_nodes, 
                                      param_edges)
  
  # setting up the manual layout
  lut_x <- layout_df$x; names(lut_x) <- layout_df$metric
  lut_y <- layout_df$y; names(lut_y) <- layout_df$metric
  
  layout_man <- create_layout(param_graph1, layout = "linear") %>%
    mutate(x = lut_x[metric],
           y=lut_y[metric]) %>%
    dplyr::select(x,y) %>%
    as.data.frame()
  
  # applying the manual layout to the graph objects, one for each group
  layout1 <- create_layout(param_graph1, layout = layout_man)
  dummy <- data.frame(val = c(0.9, 0.55, 0.25, -0.25,-0.55,-0.8),
                      class = names(lut_cols))
  leg_plot <-get_legend(
    ggplot(dummy, aes(x=val, y=val)) +
      geom_line(aes(color = class, lty=class), lwd=3) +
      theme_classic()+
      theme(legend.key.size = unit(1.5,"cm"),
            legend.text = element_text(size=15),
            legend.title = element_text(size=15),
            legend.background = element_rect(fill="transparent"))+
      guides(color=guide_legend(ncol=1))+
      scale_linetype_manual(values = lut_lty,
                            name = "Standardized\nEstimates")+
      scale_color_manual(values = lut_cols, 
                            name = "Standardized\nEstimates")
    )
  
  if(just_leg==TRUE){return(leg_plot)}
  
  p1_title <- variable[1]
  
  # Plot
  p1 <- ggraph(layout1) +
    geom_edge_arc(aes(color=class, #width = abs(val),
                      linetype = class), width=1,
                  strength = 0.1,
                  angle_calc = "along", vjust = -.5,family = 'Times',
                  check_overlap = TRUE,
                  arrow = arrow(25, length = unit(0.3, "inches"), type = "open"),
                  label_colour = "grey20",
                  end_cap = circle(0.5, "inches"),
                  start_cap = circle(0.5, "inches")
                  )+
    geom_node_text(aes(label = metric),family = 'Times', fontface = "bold", #node names
                   nudge_y = 0.05, size = 10) +
    scale_edge_color_manual(values = lut_cols, 
                            name = "Standardized\nEstimates") +
    scale_edge_linetype_manual(values = lut_lty, 
                               name = "Standardized\nEstimates") +
    scale_alpha_manual(values = c(1,0.5,0.5,0.5,0.5,1)) +
    scale_edge_width_continuous(guide = FALSE, range = c(.5,4))+
    scale_size(guide = FALSE) +
    xlim(c(-1,1))+
    theme_graph(fg_text_colour = 'white', 
                base_family = 'Times')+
    ggtitle(p1_title) +
    theme(plot.title = element_text(size = 30),
          legend.position = "none")
  
  if(legend==TRUE){plot_main <- ggarrange(p1,leg_plot, nrow=1,widths = c(3,1)) +
      ggsave(filename = filename, width = 12, height =8)
  }else{
        plot_main<-p1+
          ggsave(filename = filename, width = 12, height =8)
      }
 
  
  return(plot_main)
}

ggsem_boot <- function(boot_fit,fit, filename, variable, layout_df, r2 = FALSE,
                       r2_boot = NULL,
                       booted=TRUE, alpha = 0.05,legend=TRUE, just_leg = FALSE, 
                       values = FALSE) {

    # Extract standardized parameters
  params <- lavaan::standardizedSolution(fit)
  
  param_edges_boot <- boot_fit %>%
    as.data.frame() %>%
    pivot_longer(cols = names(.),names_to = "var", values_to = "vals") %>%
    group_by(var) %>%
    dplyr::summarise(ci_025 = HDInterval::hdi(vals)%>% pluck(1), # default is 95
                     ci_975 = HDInterval::hdi(vals)%>% pluck(2),
                     val = quantile(vals, probs = 0.5)) %>%
    ungroup() %>%
    mutate(type = ifelse(str_detect(var, "~"), "regression", "user"),
           sig = ifelse(ci_025 *ci_975 > 0, "*", "")) %>%
    filter(type == "regression" & sig == "*") %>%
    tidyr::separate(var, into = c("to", "from"), sep = "~")
  
  # Edge properties
  
  param_edges <- params %>% 
    filter(op %in% c("=~", "~", "~~"), lhs != rhs) %>% #, pvalue < .10) %>%
    mutate(sig = replace(pvalue, pvalue > alpha, 2)) %>%
    mutate(sig = replace(sig, pvalue < alpha, 1)) %>%
    filter(sig != 2) %>%
    transmute(to = lhs,
              from = rhs,
              sig=sig,
              val = est.std,
              type = dplyr::case_when(
                op == "~"  ~ "regression",
                op == "~~" ~ "correlation",
                TRUE ~ NA_character_))
  
  lut_cols<-c("< -0.75" = "red",
              "-0.5 - -0.75" = "orange",
              "0 - -0.5" = "peachpuff", 
              "0 - 0.5" = "grey80", 
              "0.5 - 0.75" = "grey50", 
              "> 0.75"= "black")
  lut_lty<-c("< -0.75" = 6,
             "-0.5 - -0.75" = 6,
             "0 - -0.5" = 6, 
             "0 - 0.5" = 1, 
             "0.5 - 0.75" = 1, 
             "> 0.75"= 1)
  
  param_edges <- param_edges %>%
    mutate(class = val)%>%
    mutate(class = replace(class, val >= 0.75, "> 0.75")) %>%
    mutate(class = replace(class, val >= 0.5 & val< 0.75, "0.5 - 0.75")) %>%
    mutate(class = replace(class, val >= 0 & val< 0.5, "0 - 0.5")) %>%
    mutate(class = replace(class, val >= -0.5 & val< 0, "0 - -0.5")) %>%
    mutate(class = replace(class, val >= -0.75 & val < -0.5, "-0.5 - -0.75")) %>%
    mutate(class = replace(class, val <= -0.75, "< -0.75")) %>%
    mutate(class = factor(class, levels = c("> 0.75",
                                            "0.5 - 0.75",
                                            "0 - 0.5",
                                            "0 - -0.5",
                                            "-0.5 - -0.75",
                                            "< -0.75"))) %>%
    mutate(sign = ifelse(val>0, 1,2))
  param_edges_boot <- param_edges_boot %>%
    mutate(class = val)%>%
    mutate(class = replace(class, val >= 0.75, "> 0.75")) %>%
    mutate(class = replace(class, val >= 0.5 & val< 0.75, "0.5 - 0.75")) %>%
    mutate(class = replace(class, val >= 0 & val< 0.5, "0 - 0.5")) %>%
    mutate(class = replace(class, val >= -0.5 & val< 0, "0 - -0.5")) %>%
    mutate(class = replace(class, val >= -0.75 & val < -0.5, "-0.5 - -0.75")) %>%
    mutate(class = replace(class, val <= -0.75, "< -0.75")) %>%
    mutate(class = factor(class, levels = c("> 0.75",
                                            "0.5 - 0.75",
                                            "0 - 0.5",
                                            "0 - -0.5",
                                            "-0.5 - -0.75",
                                            "< -0.75"))) %>%
    mutate(sign = ifelse(val>0, 1,2))
  
  sigvars <- c((param_edges_boot$to), (param_edges_boot$from))
  
  # Node properties
  param_nodes <- params %>% 
    filter(lhs == rhs) %>% 
    transmute(metric = lhs, e = est.std) %>%
    filter(metric %in% sigvars)
  
  if(booted == TRUE) {param_edges <- param_edges_boot}
  
  # Complete Graph Object
  param_graph1 <- tidygraph::tbl_graph(param_nodes, 
                                       param_edges)
  
  lut_x <- layout_df$x; names(lut_x) <- layout_df$metric
  lut_y <- layout_df$y; names(lut_y) <- layout_df$metric
    
  layout_man <- create_layout(param_graph1, layout = "linear") %>%
    mutate(x = lut_x[metric],
           y=lut_y[metric]) %>%
    dplyr::select(x,y) %>%
    as.data.frame()
  
  # applying the manual layout to the graph objects, one for each group
  layout1 <- create_layout(param_graph1, layout = layout_man) %>%
    mutate(metric = str_replace_all(metric, c("aet" = "AET")) %>% as.factor) %>%
    mutate(metric = recode_factor(metric,`p2` = 'P[ant]', 
                                  `sd_cwd` = 'sigma[CWD]',
                                  `sC` = 'C[soil]',
                                  `sN` = 'N[soil]', 
                                  `tmn` = 'T[min]',
                                  `lCN` = 'C:N[litter]'))
    
 
  
  leg_plot <-get_legend(ggplot(param_edges, aes(x=val, y=sig)) +
                          geom_line(aes(color = class, lty=class), lwd=3) +
                          scale_linetype_manual(values = c(1,1,1,6,6,6),
                                                name = "Standardized\nEstimates")+
                          theme_classic()+
                          theme(legend.key.size = unit(1.5,"cm"),
                                legend.text = element_text(size=15),
                                legend.title = element_text(size=15),
                                legend.background = element_rect(fill="transparent"))+
                          guides(color=guide_legend(ncol=1))+
                          scale_color_manual(values = lut_cols, 
                                             name = "Standardized\nEstimates"))
  
  if(just_leg==TRUE){return(leg_plot)}
  
  p1_title <- variable
  
  # Plot
  p1 <- ggraph(layout1) +
    geom_edge_arc(aes(color=class, width = abs(val),
                      linetype = class, #label=round(val,2)
                      ), #width=1,
                  strength = 0.1,
                  angle_calc = "along", vjust = -.5,family = 'Times',
                  check_overlap = TRUE,label_parse = TRUE,
                  arrow = arrow(25, length = unit(0.3, "inches"), type = "open"),
                  label_colour = "grey20", label_size = 6,
                  end_cap = circle(0.5, "inches"),
                  start_cap = circle(0.5, "inches")
    )+
    geom_node_text(aes(label = metric),family = 'Times', fontface = "bold", #node names
                   nudge_y = 0.05, size = 10,parse = T) +

    scale_edge_color_manual(values = lut_cols, 
                            name = "Standardized\nEstimates") +
    scale_edge_linetype_manual(values = lut_lty, 
                               name = "Standardized\nEstimates") +
    scale_alpha_manual(values = c(1,0.5,0.5,0.5,0.5,1)) +
    scale_edge_width_continuous(guide = FALSE, range = c(1,3))+
    scale_size(guide = FALSE) +
    xlim(c(-1,1))+
    theme_graph(fg_text_colour = 'white', 
                base_family = 'Times')+
    ggtitle(p1_title) +
    theme(plot.title = ggtext::element_markdown(size = 30),
          legend.position = "none")
  
  if(r2 == TRUE){
    r2_tab0 <- r2_boot %>%
      as_tibble() %>%
      pivot_longer(cols = names(.),names_to = "metric", values_to = "r2") %>%
      group_by(metric) %>%
      dplyr::summarise(r2=median(r2) %>% round(2)) %>%
      ungroup()%>%
      mutate(metric = recode_factor(metric,`p2` = 'P[ant]', 
                                    `sd_cwd` = 'sigma[CWD]',
                                    `sC` = 'C[soil]',
                                    `sN` = 'N[soil]', 
                                    `tmn` = 'T[min]',
                                    `lCN` = 'C:N[litter]'))
    test <- data.frame(x = layout1$x, y = layout1$y, metric = layout1$metric)
    r2_tab<- left_join(x=test,y=r2_tab0, by="metric") %>%
      mutate(y=y-0.15)
    
    p1 <- p1 +
      geom_node_text(data = r2_tab, aes(label = r2),
                     family = 'Times', fontface = "plain",
                     nudge_y = 0.05, size = 6,parse = T)
  }

  if(legend==TRUE){
    plot_main <- ggarrange(p1,leg_plot, nrow=1,widths = c(3,1)) 
    }else{
    plot_main<-p1
    }
  ggsave(plot = plot_main,filename = filename, width = 12, height =8)
  
  
  return(plot_main)
}

ggsem_boot_ga <- function(boot_fit,fit, filename, variable, layout_df, r2 = FALSE,
                       r2_boot = NULL,
                       booted=TRUE, alpha = 0.05,legend=TRUE, just_leg = FALSE, 
                       values = FALSE) {
  
  # Extract standardized parameters
  params <- lavaan::standardizedSolution(fit)
  
  param_edges_boot <- boot_fit %>%
    as.data.frame() %>%
    pivot_longer(cols = names(.),names_to = "var", values_to = "vals") %>%
    group_by(var) %>%
    dplyr::summarise(ci_025 = HDInterval::hdi(vals)%>% pluck(1), # default is 95
                     ci_975 = HDInterval::hdi(vals)%>% pluck(2),
                     val = quantile(vals, probs = 0.5)) %>%
    ungroup() %>%
    mutate(type = ifelse(str_detect(var, "~"), "regression", "user"),
           sig = ifelse(ci_025 *ci_975 > 0, "*", "")) %>%
    filter(type == "regression" & sig == "*") %>%
    tidyr::separate(var, into = c("to", "from"), sep = "~")
  
  # Edge properties
  
  param_edges <- params %>% 
    filter(op %in% c("=~", "~", "~~"), lhs != rhs) %>% #, pvalue < .10) %>%
    mutate(sig = replace(pvalue, pvalue > alpha, 2)) %>%
    mutate(sig = replace(sig, pvalue < alpha, 1)) %>%
    filter(sig != 2) %>%
    transmute(to = lhs,
              from = rhs,
              sig=sig,
              val = est.std,
              type = dplyr::case_when(
                op == "~"  ~ "regression",
                op == "~~" ~ "correlation",
                TRUE ~ NA_character_))
  
  lut_cols<-c("Positive" = "black",
              "Negative" = "red")
  
  param_edges <- param_edges %>%
    mutate(class = val)%>%
    mutate(class = replace(class, val >= 0, "Positive")) %>%
    mutate(class = replace(class, val < 0, "Negative")) %>%
    mutate(class = factor(class, levels = c("Positive", "Negative"))) %>%
    mutate(sign = ifelse(val>0, 1,2))
  
  param_edges_boot <- param_edges_boot %>%
    mutate(class = val)%>%
    mutate(class = replace(class, val >= 0, "Positive")) %>%
    mutate(class = replace(class, val < 0, "Negative")) %>%
    mutate(class = factor(class, levels = c("Positive", "Negative"))) %>%
    mutate(sign = ifelse(val>0, 1,2))
  
  sigvars <- c((param_edges_boot$to), (param_edges_boot$from))
  
  # Node properties
  param_nodes <- params %>% 
    filter(lhs == rhs) %>% 
    transmute(metric = lhs, e = est.std) %>%
    filter(metric %in% sigvars)
  
  if(booted == TRUE) {param_edges <- param_edges_boot}
  
  # Complete Graph Object
  param_graph1 <- tidygraph::tbl_graph(param_nodes, 
                                       param_edges)
  
  lut_x <- layout_df$x; names(lut_x) <- layout_df$metric
  lut_y <- layout_df$y; names(lut_y) <- layout_df$metric
  
  layout_man <- create_layout(param_graph1, layout = "linear") %>%
    mutate(x = lut_x[metric],
           y=lut_y[metric]) %>%
    dplyr::select(x,y) %>%
    as.data.frame()
  
  # applying the manual layout to the graph objects, one for each group
  layout1 <- create_layout(param_graph1, layout = layout_man) %>%
    mutate(metric = str_replace_all(metric, c("aet" = "AET")) %>% as.factor) %>%
    mutate(metric = recode_factor(metric,`p2` = 'P[ant]', 
                                  `sd_cwd` = 'sigma[CWD]',
                                  `sC` = 'C[soil]',
                                  `sN` = 'N[soil]', 
                                  `tmn` = 'T[min]',
                                  `lCN` = 'LCN',
                                  `NF` = 'NF',
                                  `BSC` = 'Crust'))
  
  
  
  leg_plot <-get_legend(ggplot(param_edges, aes(x=val, y=sig)) +
                          geom_line(aes(color = class, lty=class), lwd=3) +
                          scale_linetype_manual(values = c(1,1,1,6,6,6),
                                                name = "Standardized\nEstimates")+
                          theme_classic()+
                          theme(legend.key.size = unit(1.5,"cm"),
                                legend.text = element_text(size=15),
                                legend.title = element_text(size=15),
                                legend.background = element_rect(fill="transparent"))+
                          guides(color=guide_legend(ncol=1))+
                          scale_color_manual(values = lut_cols, 
                                             name = "Standardized\nEstimates"))
  
  if(just_leg==TRUE){return(leg_plot)}
  
  p1_title <- variable
  
  # Plot
  p1 <- ggraph(layout1) +
    geom_edge_arc(aes(color=class, #width = abs(val),
                      #linetype = class, #label=round(val,2)
    ), width=3,
    strength = 0.1,
    angle_calc = "along", vjust = -.5,family = 'Times',
    check_overlap = TRUE,label_parse = TRUE,
    arrow = arrow(25, length = unit(0.3, "inches"), type = "open"),
    label_colour = "grey20", label_size = 6,
    end_cap = circle(0.5, "inches"),
    start_cap = circle(0.5, "inches")
    )+
    geom_node_text(aes(label = metric),family = 'Times', fontface = "bold", #node names
                   nudge_y = 0.05, size = 10,parse = T) +
    
    scale_edge_color_manual(values = lut_cols, 
                            name = "Standardized\nEstimates") +
    scale_edge_linetype_manual(values = lut_lty, 
                               name = "Standardized\nEstimates") +
    scale_alpha_manual(values = c(1,0.5,0.5,0.5,0.5,1)) +
    scale_edge_width_continuous(guide = FALSE, range = c(1,3))+
    scale_size(guide = FALSE) +
    xlim(c(-1,1))+
    theme_graph(fg_text_colour = 'white', 
                base_family = 'Times')+
    ggtitle(p1_title) +
    theme(plot.title = ggtext::element_markdown(size = 30),
          legend.position = "none")
  
  if(legend==TRUE){
    plot_main <- ggarrange(p1,leg_plot, nrow=1,widths = c(3,1)) 
  }else{
    plot_main<-p1
  }
  ggsave(plot = plot_main,filename = filename, width = 12, height =8)
  
  
  return(plot_main)
}
