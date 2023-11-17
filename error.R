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

a0 <- a0 %>% filter(CAS_RN %in% c("TSS"))
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



# SEDIMENT SETUP ----
library(tidyverse)
library(readxl)

source('scripts_original/ci_band_lm_kirkedit2.R')
source('scripts_original/ci_band_plot.R')
source('scripts_original/pus.R')
source('scripts_original/ci_boot_kirkedit2.R')
source('scripts_original/impute_km.R')

site <- 'Tacoma'; site.tag <- 'Tacoma'; site.filetag <- 'Tacoma_Sediment'
user.cut <- '2022-01-01'
new.yr <- '2022'
yr.cut <- '2022-01-01'

fdate <- format.Date(Sys.Date(),format='%y%m%d')

# SEDIMENT WRANGLE ----

a0 <- read_excel('data_raw/AllWYSedTrap2021_kc.xlsx',col_types = c(rep('text',10),'date',rep('text',6),'numeric','text',rep('numeric',4),rep('text',3),rep('numeric',2),'text'))
a00 <- a0 %>% dplyr::select(-c(data_provider,matrix_desc,matrix_code,sample_name,sample_class,fraction_desc,MTCASTAT_input,half_ND,log_halfND,pthal_sum,for_stats_5yr,`qual code 2013`))
a00 <- a00 %>% mutate(id=1:n())
a00 <- a00 %>% dplyr::mutate(nd=1-detect_flag) %>% rename(conc=result_numeric,coc=chemical_name,units=result_unit,qual=lab_qualifiers) %>% dplyr::select(-detect_flag)
a00 <- a00 %>% mutate(tmp=sub('Water Year ','',task_code),tst.date=as.Date(paste0(tmp,'-07-01')))
a00 <- a00 %>% mutate(date=as.Date(logdate),date=case_when(is.na(date)~tst.date,TRUE~date)) %>% dplyr::select(-c(logdate,tmp,tst.date))
a00 <- a00 %>% rename(matrix=matrix_class_desc)
a00 <- a00 %>% mutate(nd=case_when(is.na(nd)~1,TRUE~nd))
a01 <- a00 %>% mutate(locid=case_when(locid %in% c('FD3-New','FD3-NEW')~'FD-3NEW',TRUE~locid),units=case_when(units %in% c('ug/kg dry','ug/Kg','ug/Kg dry','ug/L')~'ug/kg',units %in% c('mg/kg dry','mg/Kg','mg/Kg dry')~'mg/kg',units=='ng/Kg'~'ng/kg',TRUE~units))
a01 <- a01 %>% mutate(coc=case_when(coc=='Diethyl phthalate'~'Diethylphthalate',coc=='Endosulfan Sulfate'~'Endosulfan sulfate',coc=='Endrin Aldehyde'~'Endrin aldehyde',coc=='Endrin Ketone'~'Endrin ketone',coc=='Indeno(1,2,3-cd)pyrene'~'Indeno(1,2,3-c,d)pyrene',coc=='Solids-Total Volatile'~'Total Volatile Solids',coc=='Benzo(b,k)fluoranthenes'~'Benzo(b,k)fluoranthene',coc=="4,4'-DDD"~'4,4-DDD',coc=="4,4'-DDE"~'4,4-DDE',coc=="4,4'-DDT"~'4,4-DDT',coc=='Di-n-butyl phthalate'~'Di-n-butylphthalate',coc=='Di-n-octyl phthalate'~'Di-n-octylphthalate',TRUE~coc))
a01 <- a01 %>% mutate(units=case_when(coc=='Clay/Silt'~'%',coc=='Diesel'~'ug/kg',coc=='Gravel'~'%',coc=='Heavy Oil Range Hydrocarbons'~'ug/kg',coc=='Lead'~'ug/kg',coc=='Mercury'~'ug/kg',coc=='Sand'~'%',coc=='Total Organic Carbon'~'%',coc=='Total Solids'~'%',coc=='Zinc'~'ug/kg',coc=='Total Volatile Solids'~'%',TRUE~units))
a02 <- a01 %>% mutate(conc=case_when(coc=='2-Fluorobiphenyl' & units == 'mg/kg'~conc*1000,coc=='Terphenyl-d14' & units=='mg/kg'~conc*1000,TRUE~conc),units=case_when(coc=='2-Fluorobiphenyl'~'ug/kg',coc=='Terphenyl-d14'~'ug/kg',TRUE~units))

# link sediment traps to specific outfalls; remove <DA-1 Line>, harmonize location names
a03 <- a02 %>% filter(locid != 'DA-1 Line') %>% mutate(locid=case_when(locid == 'FD3-A'~'FD-3A',locid == 'FD6-B'~'FD-6B',locid == 'FD1'~'FD-1',locid == 'FD2'~'FD-2',locid == 'FD22'~'FD-22',locid == 'FD23'~'FD-23',locid == 'FD16'~'FD-16',locid == 'FD18'~'FD-18',locid == 'FD10-C'~'FD-10C',locid == 'FD13'~'FD-13',locid == 'FD3-New'~'FD-3NEW',locid == 'FD13-B New'~'FD-13BNEW',locid == 'FD30'~'FD-30',locid == 'FD32'~'FD-32',locid == 'FD33'~'FD-33',locid == 'FD36'~'FD-36',locid == 'FD37'~'FD-37',locid == 'FD38'~'FD-38',TRUE~locid))
a03 <- a03 %>% mutate(locid=case_when(locid %in% c('FD-3','FD-3NEW')~'FD-3NEW',locid %in% c('FD-13B','FD-13B New','FD-13BNEW')~'FD-13BNEW',TRUE~locid))
loc.vec <- a03 %>% pus(locid)
of.tab <- tibble(locid=loc.vec,outfall=c('OF237B',rep('OF237A',5),rep('OF230',4),'OF237A','OF245','OF248','OF243','OF237A',rep('OF237B',9),rep('OF230',3),'OF237A',rep('OF235',3),'OF245'))
of.tab <- of.tab %>% arrange(outfall,locid)

# arrange sediment traps by outfall grouping
a03 <- a03 %>% left_join(.,of.tab)
a03 <- a03 %>% mutate(locid=factor(locid,levels=of.tab$locid,ordered = T))

# eliminate cocs with insufficient data
tmp <- a03 %>% group_by(coc,units) %>% summarise(mx.d=max(date),N=n(),Nloc=length(unique(locid)),n.per.loc=round(N/Nloc,1),ave=mean(conc),med=median(conc),min=min(conc),max=max(conc),.groups = 'drop')
coc.elim <- tmp %>% filter(n.per.loc < 4.5 | mx.d < ymd('2021-01-01')) %>% pus(coc)
a04 <- a03 %>% filter(!coc %in% coc.elim)

f <- a04 %>%
  filter(coc == "Anthracene")
#TODO: In the sediment runs, filtering for some contaminants will produce results, 
# but the entire run will not go. 


## SEDIMENT Linear Confidence Bands ----
outfall.pair.cband <- f %>% group_by(coc,outfall) %>% do({
  td <- .
  tcoc <- td$coc[1]; tout <- td$outfall[1]
  print(paste(tcoc,tout))
  ci_band_lm(td,clev=0.99,side='both')
}) %>% ungroup()

## plot confidence bands
fn <- paste0(site.filetag,'_cband_',fdate,'.pdf')
pdf(file=fn,w=11,h=8.5)
outfall.cband.plots <- outfall.pair.cband %>% group_by(coc) %>% do(plot={
  tcoc= .$coc[1]
  print(tcoc)
  hdr <- paste0('99% Confidence Bands for Traps Grouped by Outfall for ',tcoc)
  td <- f %>% filter(coc==tcoc)
  ci_band_plot(td,.,hdr,vloc=outfall,show.lim=F,pt.size = 2,lab.size = 3)
})
dev.off()

## Last 5 year bands and plots
outfall.pair.cband.5y <- f %>% filter(date >= ymd('2017-01-01')) %>% group_by(coc,outfall) %>% do({
  td <- .
  tcoc <- td$coc[1]; tout <- td$outfall[1]
  print(paste(tcoc,tout))
  ci_band_lm(td,clev=0.99,side='both')
}) %>% ungroup()

fn <- paste0(site.filetag,'_cband5y_',fdate,'.pdf')
pdf(file=fn,w=11,h=8.5)
outfall.cband.plots.5y <- outfall.pair.cband.5y %>% group_by(coc) %>% do(plot={
  tcoc= .$coc[1]
  print(tcoc)
  hdr <- paste0('99% Confidence Bands for Traps Grouped by Outfall Since 2017 for ',tcoc)
  td <- f %>% filter(coc==tcoc,date >= ymd('2017-01-01'))
  ci_band_plot(td,.,hdr,vloc=outfall,show.lim=F,pt.size = 2,lab.size = 3)
})
dev.off()
