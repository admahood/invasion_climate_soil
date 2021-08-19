# invasion_climate_soil



[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5176664.svg)](https://doi.org/10.5281/zenodo.5176664)


This is the repository that houses the data and code to reproduce the analysis of our upcoming paper "Interannual climate variability mediates changes in carbon and nitrogen pools caused by annual grass invasion in a semi-arid shrubland"

The script to create figures 1 and 2 is `R/c_js_figs.R`, the script for figures 3 and 4 is `R/tidymods.R`.

`R/ggplot_sem.R` contains functions to create figures from lavaan objects, and these are used to create figure 5 in `R/d_sem_mods.R`. An older version of the manuscript had path models for plant C:N ratios as well, and the code for those is still in `R/d_sem_mods.R`.

Things might be a little less than straightforward, but the data are all there. Feel free to reach out with any questions!
