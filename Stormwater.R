# Tacoma Stormwater Analysis, December 2022


# SETUP ----

library(tidyverse)
library(readxl)
library(lubridate)
library(R.utils)
library(coin)
library(rlang)
library(purrr)
library(broom)
library(stats)
library(lmPerm)
library(boot)
library(magrittr)
library(RColorBrewer)
library(raster)


## R.utils::sourceDirectory(path='/Users/kcmacstat/Stats-Research/MacStat R Code',modifiedOnly=F)

source('scripts_original/eplot.R')
source('scripts_original/trend_map_compute.R')
source('scripts_original/trend_map_plot.R')
source('scripts_original/oneway_mc.R')
source('scripts_original/ci_band_lm_kirkedit3.R')
source('scripts_original/ci_band_plot_kirkedit3.R')
source('scripts_original/trend_map.R')
source('scripts_original/pus.R')
## source('add_gis.R')
## source('convex_hull.R')
source('scripts_original/oneway_mc_plot.R')
source('scripts_original/ci_compute.R')
source('scripts_original/mc_vector.R')
source('scripts_original/ci_boot_kirkedit2.R')
source('scripts_original/whdquantile.R')
source('scripts_original/wquantile_generic.R')
source('scripts_original/impute_km_kirkedit3.R')
source('scripts_original/pp_km.R')
##  source('trans.R')
##  source('backtrans.R')
##  source('pp_tricube.R')
source('scripts_original/pctdiff_boot_kirkedit.R')
source('scripts_original/Mode.R')

site <- 'Tacoma'; site.tag <- 'Tacoma'; site.filetag <- 'Tacoma_Stormwater'
user.cut <- '2022-01-01'
new.yr <- '2022'
yr.cut <- '2022-01-01'

fdate <- format.Date(Sys.Date(),format='%y%m%d')

# WRANGLE ----
a0 <- read_excel('data_raw/WY2021_lab_data_kc.xlsx',col_types = c(rep('text',3),rep('numeric',2),rep('text',3),rep('numeric',3),'date',rep('text',9),rep('numeric',9),rep('text',2)))

## TODO: Script will break when the troublesome compound is included, 
## but runs when it is excluded.
a0 <- a0 %>% filter(CAS_RN %in% c("68334-30-5")) # Breaks at location 230
# a0 <- a0 %>% filter(CAS_RN %in% c("117-81-7")) # Worksfor a single cas_rn
#a0 <- a0 %>% filter(CAS_RN != "68334-30-5", locid != "230") # Works for all cas_rns with troublesome one removed



## Continue with wrangling
names(a0) <- tolower(names(a0))

a00 <- a0 %>% dplyr::select(-c(sys_loc_code,loc_number,`_237a_comp`,task_code,sys_sample_code,cipp_lining,halfnd,loghalfnd,totphthsum,number,forstats,forstats_5yr,forstats_2yr,validator_qualifiers))

a00 <- a00 %>% mutate(id=1:n(),wt=1)
a00 <- a00 %>% mutate(nd=1-detect) %>% rename(conc=result_numeric,units=result_unit,qual=lab_qualifiers,cas=cas_rn) %>% dplyr::select(-c(detect,result_text,fraction_label))

a00 <- a00 %>% mutate(date=as.Date(logdate)) %>% dplyr::select(-c(logdate)) %>% rename(matrix=matrix_code)

a01 <- a00 %>% mutate(locid=case_when(locid == '237A New'~'237A',TRUE~locid))

a01 <- a01 %>% mutate(units=case_when(units == '#/100mL'~'#/100ml',units %in% c('CFU/100 ml','CFU/100 mL','CFU/100mL')~'cfu/100ml',units=='mg/l'~'mg/L',units=='ug/l'~'ug/L',units=='pH Units'~'SU',TRUE~units))

a01 <- a01 %>% group_by(coc) %>% mutate(units=Mode(units)) %>% ungroup() 

tst.coc <- a01 %>% group_by(coc) %>% summarise(min.d=min(date),max.d=max(date),.groups = 'drop')
coc.elim <- tst.coc %>% filter(max.d < ymd('2021-01-01') | min.d > ymd('2019-01-01')) %>% pull(coc)

a02 <- a01 %>% filter(!coc %in% coc.elim)
a02 <- a02 %>% mutate(locid=factor(locid,ordered = T))


f <- a02


# CI BANDS, PCT REDUCTIONS ----
outfall.pair.cband <- f %>% group_by(coc,locid) %>% do({
  tcoc <- .$coc[1]; tloc <- .$locid[1]
  print(paste(tcoc,tloc))
  ci_band_lm(.,clev=0.99,side='both')
}) %>% ungroup()


## plot confidence bands
fn <- paste0(site.filetag,'_cband_',fdate,'.pdf')
pdf(file=fn,w=11,h=8.5)
outfall.cband.plots <- outfall.pair.cband %>% group_by(coc) %>% do(plot={
  tcoc= .$coc[1]
  print(tcoc)
  hdr <- paste0('Confidence Bands by Outfall for ',tcoc,' With 99% Confidence')
  td <- f %>% filter(coc==tcoc)
  ci_band_plot(td,.,hdr,show.lim=F,pt.size = 2,lab.size = 3)
})
dev.off()


## Last 5 year bands and plots
outfall.pair.cband.5y <- f %>% filter(date >= ymd('2017-01-01')) %>% group_by(coc,locid) %>% do({
  tcoc <- .$coc[1]; tloc <- .$locid[1]
  print(paste(tcoc,tloc))
  ci_band_lm(.,clev=0.99,side='both')
}) %>% ungroup()



fn <- paste0(site.filetag,'_cband5y_',fdate,'.pdf')
pdf(file=fn,w=11,h=8.5)
outfall.cband.plots.5y <- outfall.pair.cband.5y %>% group_by(coc) %>% do(plot={
  tcoc= .$coc[1]
  print(tcoc)
  hdr <- paste0('Confidence Bands by Outfall Since 2017 for ',tcoc,' With 99% Confidence')
  td <- f %>% filter(coc==tcoc,date >= ymd('2017-01-01'))
  ci_band_plot(td,.,hdr,show.lim=F,pt.size = 2,lab.size = 3)
})
dev.off()

## pct changes across range of each coc-outfall pair
outfall.pctdiff <- f %>% group_by(coc,locid) %>% do({
  tcoc <- .$coc[1]; tloc <- .$locid[1]
  print(paste(tcoc,tloc))
  pctdiff_boot(.,side = 'both',impute = 'MC',ns=20,nb=100)
})


outfall.pctdiff.5y <- f %>% filter(date >= ymd('2017-01-01')) %>% group_by(coc,locid) %>% do({
  tcoc <- .$coc[1]; tloc <- .$locid[1]
  print(paste(tcoc,tloc))
  pctdiff_boot(.,side = 'both',impute = 'MC',ns=20,nb=100)
})


fn <- paste0(site.filetag,'_outfall_pctdiff_',fdate,'.csv')
write_excel_csv(outfall.pctdiff,file=fn)

fn <- paste0(site.filetag,'_outfall_pctdiff_5y_',fdate,'.csv')
write_excel_csv(outfall.pctdiff.5y,file=fn)


