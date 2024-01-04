rm(list=ls())

library(rstudioapi)
library(dplyr)
library(foreign)
library(ggplot2)
library(ggpointdensity)
library(viridis)

setwd(dirname(getActiveDocumentContext()$path))


# --------------------------------------------------------------------------------------------------------------------
# Point-source nutrient
# --------------------------------------------------------------------------------------------------------------------

# Export monthly N with location info, unit: kg/yr
facilities = read.csv('./factors/Point_SourceNut/Facilities.csv')
loadMonth = read.csv('./factors/Point_SourceNut/load_summary_by_discharger.csv') %>%
  rename('NPDES'=1, 'Nitrogen'=5, 'Phosphorous'=6) %>%
  left_join(facilities %>% select(NPDES, Latitude, Longitude), by = join_by(NPDES)) %>%
  select(NPDES, Latitude, Longitude, Nitrogen, Phosphorous) %>%
  mutate(Nitrogen=Nitrogen/12, Phosphorous=Phosphorous/12)
write.csv(loadMonth, file = './factors/Point_SourceNut/load_cal_monthly.csv')

# After adding HUC12 info, read point source load again to calculate at the HUC12 level
loadMonthHUC = read.dbf('./factors/Point_SourceNut/load_pts_shp/load_pts_monthly_huc.dbf') %>%
  select(Nitrogen, Phosphorou, huc12) %>%
  group_by(huc12) %>%
  summarise(Nitrogen = sum(Nitrogen), Phosphorou = sum(Phosphorou))

load(paste0('./analysis/table/gaugeSloc_hucGroup_TN (Water).RData'))
tnHucGroup = gaugeSloc %>% mutate(ptsNutr = NA)
load(paste0('./analysis/table/gaugeSloc_hucGroup_TP (Water).RData'))
tpHucGroup = gaugeSloc %>% mutate(ptsNutr = NA)
rm(gaugeSloc)

tnGauge = tnHucGroup$gauge[!is.na(tnHucGroup$gauge)]
for (i in 1:length(tnGauge)) {
  print(paste0(i, '/', length(tnGauge)))
  
  hucGroup = unlist(strsplit(tnHucGroup %>% filter(gauge==tnGauge[i]) %>% pull(hucGroup), ", "))
  tnHucGroup[which(tnHucGroup$gauge==tnGauge[i]),]$ptsNutr = sum(loadMonthHUC %>% filter(huc12%in%hucGroup) %>% pull(Nitrogen))
}

tpGauge = tpHucGroup$gauge[!is.na(tpHucGroup$gauge)]
for (i in 1:length(tpGauge)) {
  print(paste0(i, '/', length(tpGauge)))
  
  hucGroup = unlist(strsplit(tpHucGroup %>% filter(gauge==tpGauge[i]) %>% pull(hucGroup), ", "))
  tpHucGroup[which(tpHucGroup$gauge==tpGauge[i]),]$ptsNutr = sum(loadMonthHUC %>% filter(huc12%in%hucGroup) %>% pull(Nitrogen))
}

# The unit from LOADEST is gram/day, convert it to kg/month
tn = read.dbf('./analysis/table/wqAvgLocal_TN (Water)_Whucgroup.dbf') %>%
  left_join(tnHucGroup %>% select(huc12,ptsNutr), by = join_by(huc12)) %>%
  mutate(localLoadWoPts = localLoad*0.001*365/12 - ptsNutr)
tp = read.dbf('./analysis/table/wqAvgLocal_TP (Water)_Whucgroup.dbf') %>%
  left_join(tpHucGroup %>% select(huc12,ptsNutr), by = join_by(huc12)) %>%
  mutate(localLoadWoPts = localLoad*0.001*365/12 - ptsNutr)
tnWoNa = tn %>% na.omit() %>% select(hucGroupId,ptsNutr,localLoadWoPts)
tpWoNa = tp %>% na.omit() %>% select(hucGroupId,ptsNutr,localLoadWoPts)
tn = tn  %>% select(-ptsNutr,-localLoadWoPts) %>% left_join(tnWoNa, by = join_by(hucGroupId)) %>% filter(!is.na(hucGroupId))
tp = tp  %>% select(-ptsNutr,-localLoadWoPts) %>% left_join(tpWoNa, by = join_by(hucGroupId)) %>% filter(!is.na(hucGroupId))
write.dbf(tn, file = './factors/Point_SourceNut/wqAvgLocal_TN (Water)_Whucgroup_WoPts.dbf')
write.dbf(tp, file = './factors/Point_SourceNut/wqAvgLocal_TP (Water)_Whucgroup_WoPts.dbf')