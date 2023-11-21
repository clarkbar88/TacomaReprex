## Debug Reprex

# STORMWATER SETUP ----

library(tidyverse)
library(readxl)

source('scripts_original/ci_band_lm_kirkedit2.R')
source('scripts_original/ci_band_plot.R')
source('scripts_original/pus.R')
source('scripts_original/ci_boot_kirkedit2.R')
source('scripts_original/impute_km.R')

site <- 'Tacoma'; site.tag <- 'Tacoma'; site.filetag <- 'Tacoma_Stormwater'
user.cut <- '2022-01-01'
new.yr <- '2022'
yr.cut <- '2022-01-01'
fdate <- format.Date(Sys.Date(),format='%y%m%d')

# STORMWATER WRANGLE ----

a0 <- read_excel('data_raw/WY2021_lab_data_kc.xlsx',col_types = c(rep('text',3),rep('numeric',2),rep('text',3),rep('numeric',3),'date',rep('text',9),rep('numeric',9),rep('text',2)))
#TODO: Added a line to print the unique names of CAS_RN.  
print(unique(a0$CAS_RN))

a0 <- a0 %>% filter(CAS_RN %in% c("68334-30-5"))
#TODO: Filtering for TSS will not produce an error in the stormwater log, and the code executes normally.
#TODO: The following CAS_RN will produce the error: "68334-30-5"
# The error also occurs when there is no filtering (breaks when the code reaches that CAS_RN.)
a0 <- a0 



#Unchanged
names(a0) <- tolower(names(a0))
a00 <- a0 %>% dplyr::select(-c(sys_loc_code,loc_number,`_237a_comp`,task_code,sys_sample_code,cipp_lining,halfnd,loghalfnd,totphthsum,number,forstats,forstats_5yr,forstats_2yr,validator_qualifiers))
a00 <- a00 %>% mutate(id=1:n(),wt=1)
a00 <- a00 %>% mutate(nd=1-detect) %>% rename(conc=result_numeric,units=result_unit,qual=lab_qualifiers,cas=cas_rn) %>% dplyr::select(-c(detect,result_text,fraction_label))
a00 <- a00 %>% mutate(date=as.Date(logdate)) %>% dplyr::select(-c(logdate)) %>% rename(matrix=matrix_code)
a01 <- a00 %>% mutate(locid=case_when(locid == '237A New'~'237A',TRUE~locid))
a01 <- a01 %>% mutate(units=case_when(units == '#/100mL'~'#/100ml',units %in% c('CFU/100 ml','CFU/100 mL','CFU/100mL')~'cfu/100ml',units=='mg/l'~'mg/L',units=='ug/l'~'ug/L',units=='pH Units'~'SU',TRUE~units))
a01 <- a01 %>% group_by(coc) %>% mutate(units=mode(units)) %>% ungroup() #TODO: changed Mode to mode, unclear on functionality
tst.coc <- a01 %>% group_by(coc) %>% summarise(min.d=min(date),max.d=max(date),.groups = 'drop')
coc.elim <- tst.coc %>% filter(max.d < ymd('2021-01-01') | min.d > ymd('2019-01-01')) %>% pull(coc)
a02 <- a01 %>% filter(!coc %in% coc.elim)
a02 <- a02 %>% mutate(locid=factor(locid,ordered = T))

f <- a02

### Error reproduction
# STORMWATER CI BANDS, PCT REDUCTIONS ----
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


