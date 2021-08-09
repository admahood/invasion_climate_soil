# extract climate z-scores
library(foreach)
library(doParallel)
library(tidyverse)
library(readxl)
library(sf)
# z_path <- "/home/a/data/climate/climate_zscores"
z_path <- "data/climate_zscores"
system(paste("aws s3 sync",
             "s3://earthlab-amahood/climate/climate_zscores",
             z_path))
latlong<- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs" 

plots<-readxl::read_xlsx("data/plot locations.xlsx") %>%
  st_as_sf(coords=c("Longitude","Latitude"), crs = latlong) %>%
  dplyr::select(site_type = `Site type`, plot = `Site number`,
                elevation_m = `Elevation (m)`)
years<- 2011:2016

t0 = Sys.time()
registerDoParallel(detectCores()-1)

system(paste("echo", plots$plot[p]))

vars <- c("def", "tmn", "aet")

results <- foreach(v = vars, .combine = rbind)%dopar%{
  
  zfiles<- Sys.glob(paste0(z_path,"/", v, "_201*" ))
  r <- raster::stack(zfiles)

  output_data <- raster::extract(r, plots, df=T)%>%
    pivot_longer(cols = -ID,names_to = "variable",
                 values_to = "value") %>%
    mutate(month = str_sub(variable, 12,13) %>% str_pad(2,side="left",pad="0"),
           year = str_sub(variable, 5,8),
           date = str_c(year,"-", month,"-01") %>% as.Date(),
           variable = str_sub(variable,1,3)) %>%
    dplyr::select(plot=ID, variable, value, date)
  
  return(output_data)
}
write.csv(results, "data/js_climate_z.csv")
print(Sys.time()-t0)
