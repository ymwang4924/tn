# Validation

rm(list=ls())

library(rstudioapi)
setwd(dirname(getActiveDocumentContext()$path))
library(dplyr)
library(ggplot2)
library(ggsci)
library(ggpmisc)
library(foreign)
library(ggpointdensity)
library(viridis)
library(grid)
library(lubridate)



# ****************************************************************************************
# WRTDS (https://www.sciencebase.gov/catalog/item/60e618acd34e2a7685cf6358)
# ****************************************************************************************
wrtds20201113rn = read.csv('./validation/sparrow/Water_qualityan/yetiRun_20201113ran_data/data/WRTDS_2017data2.csv')
wrtds20210219nw = read.csv('./validation/sparrow/Water_qualityan/yetiRun_20210219new_data/data/WRTDS_2017data1124.csv')
wrtds20210219re = read.csv('./validation/sparrow/Water_qualityan/yetiRun_20210219rerun_data/data/WRTDS_2017datarerun1124.csv')
wrtds = rbind(wrtds20201113rn,wrtds20210219nw,wrtds20210219re)
rm(wrtds20201113rn,wrtds20210219nw,wrtds20210219re)
doc = wrtds %>% filter(parameter=="Dissolved organic carbon") %>% 
  mutate(conc = (ConcLow+ConcHigh)/2, var = 'doc') %>% select(Site_no,date,conc,var)
sediment = wrtds %>% filter(parameter=="Suspended sediment concentration") %>% 
  mutate(conc = (ConcLow+ConcHigh)/2, var = 'sediment') %>% select(Site_no,date,conc,var)
tn = wrtds %>% filter(parameter=="Total Nitrogen") %>% 
  mutate(conc = (ConcLow+ConcHigh)/2, var = 'tn') %>% select(Site_no,date,conc,var)
tp = wrtds %>% filter(parameter=="Total Phosphorus") %>% 
  mutate(conc = (ConcLow+ConcHigh)/2, var = 'tp') %>% select(Site_no,date,conc,var)
po4 = wrtds %>% filter(parameter%in%c("Orthophosphate, filtered","Orthophosphate, unfiltered")) %>% 
  mutate(conc = (ConcLow+ConcHigh)/2, var = 'po4') %>% select(Site_no,date,conc,var)
toc = wrtds %>% filter(parameter=="Total organic carbon") %>% 
  mutate(conc = (ConcLow+ConcHigh)/2, var = 'toc') %>% select(Site_no,date,conc,var)
poc = wrtds %>% filter(parameter=="Particulate organic carbon") %>% 
  mutate(conc = (ConcLow+ConcHigh)/2, var = 'poc') %>% select(Site_no,date,conc,var)
rm(wrtds)
wrtds = rbind(doc,sediment,tn,tp,po4,toc,poc)
para = unique(wrtds$var)

for (i in 1:length(para)) {
  print(para[i])
  
  wrtds_nutr = wrtds %>% filter(var==para[i])
  sites = unique(wrtds_nutr$Site_no)
  for (s in 1:length(sites)) {
    print(paste0(s,'/',length(sites)))
    
    wrtds_nutrSite = wrtds_nutr %>% filter(Site_no==sites[s]) %>% mutate(date=mdy(date))
    # streamflow: fts to l/s
    q = read.csv(paste0('./validation/sparrow/Water_qualityan/flow/Q_',sites[s],'.csv')) %>%
      mutate(flow=scaledQ*28.3169, date=as.Date(start_date)) %>% select(date,flow)
    # load: ton/month
    nutrWq = inner_join(wrtds_nutrSite, q, by = join_by(date)) %>% 
      mutate(load = conc*flow*86400*0.001*0.001*0.001, year = substr(date,1,4), month = substr(date,6,7), dayInMonth = days_in_month(date)) %>%
      group_by(year,month) %>%
      reframe(loadMonth = mean(load)*dayInMonth, n = n()) %>%
      distinct(year, month, .keep_all = T) %>%
      mutate(gauge = sites[s])
    
    if (s==1) nutrWqAll = nutrWq
    if (s>1) nutrWqAll = rbind(nutrWqAll, nutrWq)
  }
  write.dbf(as.data.frame(nutrWqAll), file = paste0('./validation/sparrow/Water_qualityan/',para[i],'.dbf'))
}




# Function for plot regression
plotFun = function(joinDF, pathFileName) {
  lregre = lm(joinDF$load~joinDF$meanLoad)
  r2 = round(summary(lregre)$r.squared, 3)
  rmse = sqrt(mean((joinDF$meanLoad - joinDF$load)^2))
  grobR2 = grobTree(textGrob(paste0('R square = ', round(r2, 3)), x=0.1,  y=0.85, hjust=0,
                             gp=gpar(col="red", fontsize=38, fontface="italic")))
  grobN = grobTree(textGrob(paste0('n = ', nrow(joinDF)), x=0.1,  y=0.9, hjust=0,
                            gp=gpar(col="red", fontsize=38, fontface="italic")))
  grobRmse = grobTree(textGrob(paste0('RMSE = ', round(rmse)), x=0.1,  y=0.8, hjust=0,
                               gp=gpar(col="red", fontsize=38, fontface="italic")))
  joinPlot = ggplot(joinDF, aes(x=meanLoad, y=load)) +
    geom_abline(color = 'black', size = 1) +
    geom_point(size = 8, shape = 21, color=pal_npg("nrc")(2)[2], fill=pal_npg("nrc", alpha=0.8)(2)[2], stroke=1) + 
    labs(x=expression('LOADEST (MgN yr'^-1*')'), y = expression('WRTDS (MgN yr'^-1*')')) +
    scale_color_viridis() +
    coord_cartesian(ylim = c(min(c(joinDF$meanLoad,joinDF$load), na.rm = T),max(c(joinDF$meanLoad,joinDF$load), na.rm = T)),
                    xlim = c(min(c(joinDF$meanLoad,joinDF$load), na.rm = T),max(c(joinDF$meanLoad,joinDF$load), na.rm = T))) +
    geom_smooth(method=lm, color = pal_npg("nrc")(2)[1], se = F, size = 1, linetype = "dashed") +
    annotation_custom(grobR2) +
    annotation_custom(grobN) +
    annotation_custom(grobRmse) +
    theme(legend.position = 'none', legend.justification='left',
          plot.title = element_text(hjust = 0, size = 28, face = "bold"),
          axis.text = element_text(size = 38, colour = 'black'),
          axis.text.x = element_text(hjust = 0.9),
          axis.text.y = element_text(vjust = 0.6),
          axis.ticks.length = unit(-0.25, "cm"),
          axis.line = element_line(colour = "black", size = 1, linetype = "solid"),
          panel.grid.major.y = element_line(colour = "black", linetype = "dashed"),
          panel.grid.major.x = element_line(colour = "black", linetype = "dashed"),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.title = element_text(size=42),
          legend.text = element_text(size=28),
          legend.title = element_text(size=28),
          panel.border = element_rect(colour = "black", fill=NA, size=2))
  
  ggsave(pathFileName, joinPlot, width = 16, height = 14.4)
}



# Monthly average load in each station
level = '1'
for (i in 1:length(para)) {
  print(para[i])
  
  switch (level,
          '1' = {
            # Annual average load of multi-year data from WRTDS
            nutrWqAll = read.dbf(paste0('./validation/sparrow/Water_qualityan/',para[i],'.dbf')) %>%
              group_by(gauge) %>%
              reframe(load = mean(loadMonth)*12)
            
            # Annual average load from LOADEST
            nutrType = switch (para[i],
                               "doc" = "Organic carbon_Dissolved",
                               "sediment" = "sediment (suspended solid)",
                               "tn" = "TN (Water)",
                               "tp" = "TP (Water)",
                               "po4" = "PO4",
                               "toc" = "Organic carbon_Total",
                               "poc" = "Organic carbon_Suspended"
            )
            nutr = read.dbf(paste0('./analysis/table/wqAvgLocal_',nutrType,'.dbf')) %>% 
              mutate(meanLoad = meanLoad*0.001*0.001*365) %>% 
              select(site_no,meanLoad)
            
            # Join by station number
            joinDF = inner_join(nutrWqAll, nutr, by = join_by(gauge==site_no))
          },
          '2' = {
            # Monthly load from WRTDS
            nutrWqAll = read.dbf(paste0('./validation/sparrow/Water_qualityan/',para[i],'.dbf')) %>% 
              mutate(gauge = as.numeric(as.character(gauge)),
                     year = as.numeric(as.character(year)),
                     month = as.numeric(as.character(month))) %>% 
              na.omit()
            
            # Monthly load from LOADEST
            nutrType = switch (para[i],
                               "doc" = "Organic carbon_Dissolved",
                               "sediment" = "sediment (suspended solid)",
                               "tn" = "TN (Water)",
                               "tp" = "TP (Water)",
                               "po4" = "PO4",
                               "toc" = "Organic carbon_Total",
                               "poc" = "Organic carbon_Suspended"
            )
            nutr = read.csv(paste0('./result_dataset/0nutrient complie/',nutrType,'.csv')) %>% 
              mutate(month = case_when(month=='Jan.'~1, month=='Feb.'~2, month=='Mar.'~3, month=='Apr.'~4, month=='May'~5, month=='June'~6,
                                       month=='July'~7, month=='Aug.'~8, month=='Sep.'~9, month=='Oct.'~10, month=='Nov.'~11, month=='Dec.'~12),
                     meanLoad = meanLoad*0.001*0.001*365/12) %>%
              select(station,year,month,meanLoad)
            
            # Join by station number, year, and month
            joinDF = inner_join(nutrWqAll, nutr, by = join_by(gauge==station, year==year, month==month)) %>% rename('load'=3)
          }
  )
  
  plotFun(joinDF, paste0('./validation/sparrow/Water_qualityan/',para[i],'_annual.png'))
}




# ****************************************************************************************
# correct tohuc info and evaluate drainage area
# ****************************************************************************************
load('./shapefiles/huc/HUC12_traceUp_v3.RData')
usgs_raw = read.csv('./validation/USGS_drainage/drainageArea.csv') %>% select(site_no,site_tp_cd)

usgs = read.dbf('./validation/USGS_drainage/usgs_drainage_huc12.dbf') %>%
  select(site_no,drain_area,drainSqKm,huc12) %>%
  rename('HUC12'=4)
usgs = inner_join(usgs, usgs_raw, by = join_by(site_no)) %>%
  filter(site_tp_cd=='ST') %>% select(site_no,drain_area,drainSqKm,HUC12)
evlt = inner_join(traceUpAll, usgs, by=join_by('HUC12')) %>%
  mutate(drainage = drainage/2.59, 
         huc2usgs = drainage/drain_area,
         usgs2huc = drain_area/drainage) %>%
  select(-upHUC12)

# Revise drainage area based on NHDflowline
evlt[which(evlt$site_no=='8355490'),]$HUC12 = '130202031104'
evlt[which(evlt$site_no=='8355490'),]$drainage = traceUpAll[which(traceUpAll$HUC12=='130202031104'),]$drainage/2.59
evlt[which(evlt$site_no=='8355490'),]$upArea = traceUpAll[which(traceUpAll$HUC12=='130202031104'),]$upArea/2.59

evlt[which(evlt$site_no=='8473700'),]$HUC12 = '130900020403'
evlt[which(evlt$site_no=='8473700'),]$drainage = traceUpAll[which(traceUpAll$HUC12=='130900020403'),]$drainage/2.59
evlt[which(evlt$site_no=='8473700'),]$upArea = traceUpAll[which(traceUpAll$HUC12=='130900020403'),]$upArea/2.59

evlt[which(evlt$site_no=='9419000'),]$HUC12 = '150100120904'
evlt[which(evlt$site_no=='9419000'),]$drainage = traceUpAll[which(traceUpAll$HUC12=='150100120904'),]$drainage/2.59
evlt[which(evlt$site_no=='9419000'),]$upArea = traceUpAll[which(traceUpAll$HUC12=='150100120904'),]$upArea/2.59

# ?
evlt[which(evlt$site_no=='2172002'),]$HUC12 = '030502010101'
evlt[which(evlt$site_no=='2172002'),]$drainage = traceUpAll[which(traceUpAll$HUC12=='030502010101'),]$drainage/2.59
evlt[which(evlt$site_no=='2172002'),]$upArea = traceUpAll[which(traceUpAll$HUC12=='030502010101'),]$upArea/2.59

evlt[which(evlt$site_no=='4087170'),]$HUC12 = '040400030606'
evlt[which(evlt$site_no=='4087170'),]$drainage = traceUpAll[which(traceUpAll$HUC12=='040400030606'),]$drainage/2.59
evlt[which(evlt$site_no=='4087170'),]$upArea = traceUpAll[which(traceUpAll$HUC12=='040400030606'),]$upArea/2.59

evlt[which(evlt$site_no=='7141300'),]$HUC12 = '110300040702'
evlt[which(evlt$site_no=='7141300'),]$drainage = traceUpAll[which(traceUpAll$HUC12=='110300040702'),]$drainage/2.59
evlt[which(evlt$site_no=='7141300'),]$upArea = traceUpAll[which(traceUpAll$HUC12=='110300040702'),]$upArea/2.59

evlt[which(evlt$site_no=='8064100'),]$HUC12 = '120301090402'
evlt[which(evlt$site_no=='8064100'),]$drainage = traceUpAll[which(traceUpAll$HUC12=='120301090402'),]$drainage/2.59
evlt[which(evlt$site_no=='8064100'),]$upArea = traceUpAll[which(traceUpAll$HUC12=='120301090402'),]$upArea/2.59

evlt[which(evlt$site_no=='5212700'),]$HUC12 = '070101030206'
evlt[which(evlt$site_no=='5212700'),]$drainage = traceUpAll[which(traceUpAll$HUC12=='070101030206'),]$drainage/2.59
evlt[which(evlt$site_no=='5212700'),]$upArea = traceUpAll[which(traceUpAll$HUC12=='070101030206'),]$upArea/2.59

evlt[which(evlt$site_no=='5050000'),]$HUC12 = '090201010502'
evlt[which(evlt$site_no=='5050000'),]$drainage = traceUpAll[which(traceUpAll$HUC12=='090201010502'),]$drainage/2.59
evlt[which(evlt$site_no=='5050000'),]$upArea = traceUpAll[which(traceUpAll$HUC12=='090201010502'),]$upArea/2.59

evlt[which(evlt$site_no=='2045320'),]$HUC12 = '030102010401'
evlt[which(evlt$site_no=='2045320'),]$drainage = traceUpAll[which(traceUpAll$HUC12=='030102010401'),]$drainage/2.59
evlt[which(evlt$site_no=='2045320'),]$upArea = traceUpAll[which(traceUpAll$HUC12=='030102010401'),]$upArea/2.59

evlt[which(evlt$site_no=='3402900'),]$HUC12 = '051301010303'
evlt[which(evlt$site_no=='3402900'),]$drainage = traceUpAll[which(traceUpAll$HUC12=='051301010303'),]$drainage/2.59
evlt[which(evlt$site_no=='3402900'),]$upArea = traceUpAll[which(traceUpAll$HUC12=='051301010303'),]$upArea/2.59

evlt[which(evlt$site_no=='15478000'),]$HUC12 = '190803070901'
evlt[which(evlt$site_no=='15478000'),]$drainage = traceUpAll[which(traceUpAll$HUC12=='190803070901'),]$drainage/2.59
evlt[which(evlt$site_no=='15478000'),]$upArea = traceUpAll[which(traceUpAll$HUC12=='190803070901'),]$upArea/2.59

evlt[which(evlt$site_no=='10172727'),]$HUC12 = '160203040104'
evlt[which(evlt$site_no=='10172727'),]$drainage = traceUpAll[which(traceUpAll$HUC12=='160203040104'),]$drainage/2.59
evlt[which(evlt$site_no=='10172727'),]$upArea = traceUpAll[which(traceUpAll$HUC12=='160203040104'),]$upArea/2.59

evlt[which(evlt$site_no=='5592900'),]$HUC12 = '071402020508'
evlt[which(evlt$site_no=='5592900'),]$drainage = traceUpAll[which(traceUpAll$HUC12=='071402020508'),]$drainage/2.59
evlt[which(evlt$site_no=='5592900'),]$upArea = traceUpAll[which(traceUpAll$HUC12=='071402020508'),]$upArea/2.59

evlt[which(evlt$site_no=='9517490'),]$HUC12 = '150701040706'
evlt[which(evlt$site_no=='9517490'),]$drainage = traceUpAll[which(traceUpAll$HUC12=='150701040706'),]$drainage/2.59
evlt[which(evlt$site_no=='9517490'),]$upArea = traceUpAll[which(traceUpAll$HUC12=='150701040706'),]$upArea/2.59

evlt[which(evlt$site_no=='3196600'),]$HUC12 = '050500070606'
evlt[which(evlt$site_no=='3196600'),]$drainage = traceUpAll[which(traceUpAll$HUC12=='050500070606'),]$drainage/2.59
evlt[which(evlt$site_no=='3196600'),]$upArea = traceUpAll[which(traceUpAll$HUC12=='050500070606'),]$upArea/2.59

evlt[which(evlt$site_no=='4115000'),]$HUC12 = '040500050208'
evlt[which(evlt$site_no=='4115000'),]$drainage = traceUpAll[which(traceUpAll$HUC12=='040500050208'),]$drainage/2.59
evlt[which(evlt$site_no=='4115000'),]$upArea = traceUpAll[which(traceUpAll$HUC12=='040500050208'),]$upArea/2.59

evlt[which(evlt$site_no=='3150700'),]$HUC12 = '050302011010'
evlt[which(evlt$site_no=='3150700'),]$drainage = traceUpAll[which(traceUpAll$HUC12=='050302011010'),]$drainage/2.59
evlt[which(evlt$site_no=='3150700'),]$upArea = traceUpAll[which(traceUpAll$HUC12=='050302011010'),]$upArea/2.59

evlt[which(evlt$site_no=='12415070'),]$HUC12 = '170103041104'
evlt[which(evlt$site_no=='12415070'),]$drainage = traceUpAll[which(traceUpAll$HUC12=='170103041104'),]$drainage/2.59
evlt[which(evlt$site_no=='12415070'),]$upArea = traceUpAll[which(traceUpAll$HUC12=='170103041104'),]$upArea/2.59

evlt[which(evlt$site_no=='12465400'),]$HUC12 = '170200130608'
evlt[which(evlt$site_no=='12465400'),]$drainage = traceUpAll[which(traceUpAll$HUC12=='170200130608'),]$drainage/2.59
evlt[which(evlt$site_no=='12465400'),]$upArea = traceUpAll[which(traceUpAll$HUC12=='170200130608'),]$upArea/2.59

evlt[which(evlt$site_no=='13306385'),]$HUC12 = '170602031003'
evlt[which(evlt$site_no=='13306385'),]$drainage = traceUpAll[which(traceUpAll$HUC12=='170602031003'),]$drainage/2.59
evlt[which(evlt$site_no=='13306385'),]$upArea = traceUpAll[which(traceUpAll$HUC12=='170602031003'),]$upArea/2.59

evlt[which(evlt$site_no=='2171645'),]$HUC12 = '030501120106'
evlt[which(evlt$site_no=='2171645'),]$drainage = traceUpAll[which(traceUpAll$HUC12=='030501120106'),]$drainage/2.59
evlt[which(evlt$site_no=='2171645'),]$upArea = traceUpAll[which(traceUpAll$HUC12=='030501120106'),]$upArea/2.59

evlt = evlt %>% filter(site_no!='14246900') %>% filter(site_no!='6889000') %>% filter(site_no!='234296910') %>% 
  filter(site_no!='7144301') %>% filter(site_no!='4127885') %>% filter(site_no!='5056500') %>% filter(site_no!='2187010') %>% 
  filter(site_no!='7385500') %>% filter(site_no!='7368000') %>% filter(site_no!='7369649') %>% filter(site_no!='14092050') %>% 
  filter(site_no!='2188100') %>% filter(site_no!='12470000') %>% filter(site_no!='3268100') %>% filter(site_no!='4159130') %>% 
  filter(site_no!='4165710') %>% filter(site_no!='14145100') %>% filter(site_no!='15041200') %>% filter(site_no!='8079700') %>% 
  filter(site_no!='2357500') %>% filter(site_no!='15453500') %>% filter(HUC12 != '190203021704') %>% filter(site_no!='5536191') %>%
  mutate(huc2usgs = drainage/drain_area,
         usgs2huc = drain_area/drainage)

conusHUC12 = read.dbf('./shapefiles/huc/trace check/WBDHU12_conus.dbf')
evlt = evlt %>% group_by(HUC12) %>% mutate(flag1 = abs(1-usgs2huc), flag2 = abs(1-huc2usgs)) %>%
  slice(which.min(flag1)) %>% 
  filter(HUC12%in%conusHUC12$huc12) %>%
  filter(HUC12!='103002000804') %>% filter(HUC12!='071401010401') %>% filter(HUC12!='120902050307') %>%
  filter(HUC12!='071401010507') %>% filter(HUC12!='103001010305') %>% filter(HUC12!='050902030204') %>% 
  filter(HUC12!='050902030201')
rm(conusHUC12,traceUpAll,usgs,usgs_raw)

gages = read.dbf('./validation/GAGES-II/gagesII_9322_sept30_2011.dbf') %>%
  mutate(drain_area = DRAIN_SQKM/2.59) %>%
  mutate(site_no = as.numeric(levels(STAID)[as.numeric(STAID)])) %>%
  select(site_no,drain_area)
usgs = read.dbf('./validation/USGS_drainage/usgs_drainage.dbf') %>%
  inner_join(gages, by=join_by(site_no)) %>%
  rename('USGS'=5, 'GAGES'=8)

evlt = evlt %>% mutate(drain_area=drain_area*2.58998811, drainage=drainage*2.58998811)
lregre = lm(evlt$drain_area~evlt$drainage)
r2 = summary(lregre)$r.squared
rmse = sqrt(mean((evlt$drain_area - evlt$drainage)^2))
grobN = grobTree(textGrob(paste0('n = ', nrow(evlt)), x=0.1,  y=0.9, hjust=0,
                          gp=gpar(col="red", fontsize=38, fontface="italic")))
grobR2 = grobTree(textGrob(paste0('R square = ', round(r2, 3)), x=0.1,  y=0.85, hjust=0,
                           gp=gpar(col="red", fontsize=38, fontface="italic")))
grobRmse = grobTree(textGrob(paste0('RMSE = ', round(rmse)), x=0.1,  y=0.8, hjust=0,
                             gp=gpar(col="red", fontsize=38, fontface="italic")))
ggplot(evlt,aes(x=drain_area, y=drainage)) +
  geom_abline(color = 'black', size = 1) +
  geom_point(size = 8, shape = 21, color=pal_npg("nrc")(2)[2], fill=pal_npg("nrc", alpha=0.8)(2)[2], stroke=2) +
  annotation_custom(grobN) +
  annotation_custom(grobR2) +
  annotation_custom(grobRmse) +
  labs(x=expression('USGS reported drainage area (km'^2*')'), y= expression('Drainage area derived from this study (km'^2*')')) +
  geom_smooth(method=lm,color = pal_npg("nrc")(2)[1], se = FALSE, size = 1.5, linetype = "dashed") +
  coord_cartesian(ylim = c(0,max(c(evlt$drainage,evlt$drain_area), na.rm = T)),
                  xlim = c(0,max(c(evlt$drainage,evlt$drain_area), na.rm = T))) +
  theme(legend.position = 'none', legend.justification='left',
        plot.title = element_text(hjust = 0, size = 28, face = "bold"),
        axis.text = element_text(size = 38, colour = 'black'),
        axis.text.x = element_text(hjust = .8),
        axis.text.y = element_text(vjust = .8),
        axis.ticks.length = unit(-0.25, "cm"),
        axis.line = element_line(colour = "black", size = 1, linetype = "solid"),
        panel.grid.major.y = element_line(colour = "black", linetype = "dashed"),
        panel.grid.major.x = element_line(colour = "black", linetype = "dashed"),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.title = element_text(size=42),
        legend.text = element_text(size=28),
        legend.title = element_text(size=28),
        panel.border = element_rect(colour = "black", fill=NA, size=2))
ggsave(paste0('./manuscript/N P data/figures/drainageArea.png'), width = 16, height = 14.4)