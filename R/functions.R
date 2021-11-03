library(ggpubr)
library(sjPlot)
## functions -------------------------------------------------------------------

mass_krus = function(x, y){
  # This function inputs a data frame, and outputs another data frame
  # of Kruskall-Wallis test results on that data. x is the many response
  # variables, y is one column - the treatment. It runs the K-W test 
  # column by column.
  
  krus = list()
  for(i in 1:(length(x))){
    means = as.data.frame(kruskal(x[,i], y)$means)
    groups = as.data.frame(kruskal(x[,i], y)$groups)
    means$Treatment = rownames(means)
    groups$Treatment = rownames(groups)
    krus[[i]] = plyr::join(means, groups, by = "Treatment")
    #colnames(krus[[i]]) = c("variable", "something","std", "r", "Min", "Max", "trt", "means", "Groups")
    names(krus[[i]])[1] = "means"
    krus[[i]]$variable = as.factor(colnames(x[i]))
  }
  tuks_df = as.data.frame(do.call("rbind", krus))
  tuks_df$medians = round(tuks_df$Q50, 2)
  tuks_df$Treatment = as.factor(tuks_df$Treatment)
  tuks_df$Treatment = factor(tuks_df$Treatment, levels(tuks_df$Treatment)[c(3,4,1,2)])
  tuks_df$median_std = paste0(tuks_df$medians, " (", round(tuks_df$std,2), ")")
  tuks_df$groups <- factor(tuks_df$groups, levels = c("a", "ab", "b", "bc", "c"))
  #colnames(tuks_df) = c("Variable", "something","std", "r", "min", "max", "Treatment", "means", "Group", "mean_std")
  return(tuks_df)
}

plotit = function(kdf, title){
  # this function plots the results of the kw tests
  
  only_sigs = kdf %>%
    #mutate(M_id = class_(groups)) %>%
    group_by(variable) # %>%
  #summarise(min_ID = sum(M_id)) %>%
  #filter(min_ID >= 1 ) %>%
  #left_join(.,kdf, by="variable")
  
  ggplot(na.omit(only_sigs), aes(x=Treatment, y = variable)) +
    geom_tile(aes(fill = as.factor(groups)), color = "grey10") +
    theme_bw() + 
    scale_fill_grey(start = 1, end = 0.5) +
    geom_text(aes(label=median_std)) +
    ggtitle(paste(title))
}

# these are from the internet... for the correlation matrices
reorder_cormat <- function(cormat){
  # Use correlation between variables as distance
  dd <- as.dist((1-cormat)/2)
  hc <- hclust(dd)
  cormat <-cormat[hc$order, hc$order]
}
get_lower_tri<-function(cormat){
  cormat[upper.tri(cormat)] <- NA
  return(cormat)
}
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}

# these are to make the correlation matrices, they output
# a ggplot and also save as a .png and depend on the preceeding functions
tri_cor <-function(resp, meth= "pearson"){
  
  cormat <- round(cor(resp, use="pairwise.complete.obs", method = meth),2)
  
  # Get upper triangle of the correlation matrix
  # Reorder the correlation matrix
  cormat <- reorder_cormat(cormat)
  upper_tri <- get_upper_tri(cormat)
  melted_cormat <- melt(upper_tri, na.rm = TRUE)
  
  # take out low correlation values in the interest of visualization
  melted_cormat$vis = ifelse(abs(melted_cormat$value) >0.3, melted_cormat$value, NA ) 
  
  # Create a ggheatmap
  g <- ggplot(melted_cormat, aes(Var2, Var1, fill = value))+
    geom_tile(color = "white")+
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                         midpoint = 0, limit = c(-1,1), space = "Lab", 
                         name=paste0(meth,"\nCorrelation")) +
    theme_minimal() + # minimal theme
    theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                     size = 12, hjust = 1)) +
    coord_fixed() + 
    geom_text(aes(Var2, Var1, label = value), color = "grey60", size = 3) +
    geom_text(aes(Var2, Var1, label = vis), color = "black", size = 3) +
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      panel.grid.major = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank(),
      axis.ticks = element_blank(),
      legend.justification = c(1, 0),
      legend.position = c(0.6, 0.7),
      legend.direction = "horizontal")+
    scale_y_discrete(position = "right") +
    guides(fill = guide_colorbar(barwidth = 7, 
                                 barheight = 1,
                                 title.position = "top", 
                                 title.hjust = 0.5))
  ggsave("Images/cormat_soil_plant.png", scale=1.45,limitsize = FALSE)
  return(g)
}
square_cor <- function(pred, resp){
  cormat <- round(cor(pred, resp, use="pairwise.complete.obs",method="kendall"),2)
  melted_cormat <- melt(cormat)
  # take out low correlation values in the interest of visualization
  melted_cormat$vis = ifelse(abs(melted_cormat$value) >0.3, melted_cormat$value, NA ) 
  
  # Create a ggheatmap
  g <- ggplot(melted_cormat, aes(Var2, Var1, fill = value)) +
    geom_tile(color = "white") +
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                         midpoint = 0, limit = c(-1,1), space = "Lab", 
                         name="Pearson\nCorrelation") +
    theme_minimal()+ # minimal theme
    theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                     size = 12, hjust = 1))+
    coord_fixed()+ 
    geom_text(aes(Var2, Var1, label = value), color = "grey60", size = 4) +
    geom_text(aes(Var2, Var1, label = vis), color = "black", size = 4) +
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      panel.grid.major = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank(),
      axis.ticks = element_blank())
  ggsave("Images/cor_veg_soil_plant.png", scale=1.2)
  return(g)
}

# functions --------------------------------------------------------------------
best_random <- function(df){
  # this tests lm vs lme with random slope vs lme with random slope and intercept
  # f is a formula using the formula() function NEEDS TO BE DEFINED IN THE 
  # GLOBAL ENVIRONMENT 
  # df is the data frame
  ctrl <- lmeControl(opt = "optim", maxIter = 1000, msMaxIter = 1000)
  # gls <- gls(model = f, data = df,na.action = na.pass) # straight up lm - REML 
  # is default
  ri <- lme(fixed = f, 
            random = ~1|Site_number,
            na.action = na.pass,
            control = ctrl,
            method = "REML", data = df) # lme with random intercept
  ris <- lme(fixed = f, 
             random = ~1 + Year|Site_number, 
             na.action = na.pass,
             method = "REML",
             control = ctrl,
             data = df) # random slope and intercept
  anova(ri,ris)
}

plot_res_col <- function(mod, df, pred){
  # this function plots the residuals vs fitted, colored by all the variables
  Residuals <- resid(mod, type = "normalized") 
  Fitted <- fitted(mod)
  resp <- dplyr::select(df, Rock, Crypto, Shrubs, PNF, PNG, ANF, 
                        AIG, AIF, Litter, ppt)
  coefs2 <- c("Rock","Crypto","Shrubs","PNF","PNG",
              "ANF","AIG","AIF","Litter","ppt")
  dir.create(paste0("Images/model_viz_jan18/",pred))
  
  for(i in 1:length(coefs2)){
    ggplot(data=df, aes(x=Fitted,y=Residuals,colour=resp[i])) + 
      geom_point() +
      scale_color_continuous()+
      ggtitle(coefs2[i])
    
    ggsave(paste0("Images/model_viz_jan18/", pred,"/", coefs2[i], ".png"), limitsize = TRUE)
    
  }}

plot_res_size <- function(mod, df){
  Residuals <- resid(mod, type = "normalized") 
  Fitted <- fitted(mod)
  resp <- dplyr::select(df, Rock, Crypto, Shrubs, PNF, PNG, ANF, 
                        AIG, AIF, Litter, ppt)
  coefs2 <- c("Rock","Crypto","Shrubs","PNF","PNG",
              "ANF","AIG","AIF","Litter","ppt")
  
  op <- par(mfrow = c(2, 2), mar = c(4, 4, 3, 2)) 
  MyYlab <- "Residuals"
  for(i in 1:length(resp)){
    plot(x = Fitted, y = Residuals, 
         xlab = "Fitted values", 
         ylab = MyYlab, 
         cex = scale(resp[,i]),
         main = coefs2[i])
  }
  par(op)
}

plot_res_vs_all <- function(df, mod){
  # plots model residuals vs all predictors
  E2 <- resid(mod, type = "normalized") 
  F2 <- fitted(mod)
  op <- par(mfrow = c(2, 2), mar = c(4, 4, 3, 2)) 
  MyYlab <- "Residuals"
  plot(x = F2, y = E2, xlab = "Fitted values", ylab = MyYlab)
  plot(y = E2, x= df$AIG, main = "AIG", ylab = MyYlab)
  plot(y = E2, x= df$Litter, main = "Litter", ylab = MyYlab)
  plot(y = E2, x= df$AIF, main = "AIF", ylab = MyYlab)
  plot(y = E2, x= df$PNF, main = "PNF", ylab = MyYlab)
  plot(y = E2, x= df$Crypto, main = "Crypto", ylab = MyYlab)
  plot(y = E2, x= df$PNG, main = "PNG", ylab = MyYlab)
  plot(y = E2, x= df$ANF, main = "ANF", ylab = MyYlab)
  plot(y = E2, x= df$Shrubs, main = "Shrubs", ylab = MyYlab)
  plot(y = E2, x= df$Rock, main = "Rock", ylab = MyYlab)
  plot(y = E2, x= df$ppt, main = "precipitation", ylab = MyYlab)
  par(op)
  dev.off()
}

droptest <- function(mod, coefs){
  # this function runs an anova with every variable dropped
  for(i in 1:length(coefs)){
    x <- update(mod, as.formula(paste0(".~. -",coefs[i])))
    print(coefs[i])            
    print(anova(mod,x))
  }}

viz_mod <- function(mod, title, path, file){
  # plots the model coefficients and saves to png
  return(plot_model(mod, 
                    type = "est", 
                    title = title,
                    colors = "bw",
                    sort.est = TRUE,
                    show.values = TRUE,
                    show.p = TRUE))
  
  ggsave(paste0(path,file,".png"), limitsize = FALSE)
}

viz_int <- function(mod){
  return(
    sjp.int(mod, 
            swap.pred = T, 
            #int.term = "Shrubs:ppt", 
            mdrt.values = "quart",
            legend.title = "Precipitation",
            title = "",
            axis.title.y = "",
            show.ci = T,
            facet.grid = F,
            geom.colors = c("black","grey20", "grey40", "grey60", "grey80"))
  )
}


v_struct <- function(v = formula(~ppt),df,r){
  # tests different weight structures
  # f is formula for fixed effects - DEFINE IN GLOBAL ENV
  # v is formula for variable to test e.g. formula(~ppt)
  # df is data
  # r is formula for random structure i.e. formula(~1 + Year|Site_number)
  vcp <- varConstPower(form = v)
  ve <- varExp(form = v)
  vp <- varPower(form = v)
  vf <- varFixed(v)
  n <- lme(f, random = r, 
           method = "REML", 
           control = ctrl,
           data = df) 
  lvcp <- lme(f, random = r, 
              method = "REML", 
              control = ctrl,
              weights = vcp,
              data = df)
  lve <- lme(f, random = r, 
             method = "REML", 
             control = ctrl,
             weights = ve,
             data = df) 
  lvp <- lme(f, random = r, 
             method = "REML", 
             control = ctrl,
             weights = vp,
             data = df) 
  lvf <- lme(f, random = r, 
             method = "REML", 
             control = ctrl,
             weights = vf,
             data = df)
  print(anova(n,lvcp,lve,lvp,lvf))
}

set_theme(theme_pubr())
model_list = list()

pct_to_kg_ha <- function(pct, bulk_density, depth){
  #1 % = 1000 mg.kg-1 = 10 g.kg-1= 0.001 kg.kg-1
  ha= 10^8 #cm2
  soil_kg_per_ha = depth * ha * bulk_density / 1000 
  return(soil_kg_per_ha * pct/100)
}
  
