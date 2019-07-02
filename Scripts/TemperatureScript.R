### plot temperature data from Panama and Bermuda
## created by Nyssa Silbiger
## edited on 11/30/2018


#############
library(tidyverse)
library(lubridate)
library(sf)
library(gridExtra)
library(grid)

# Panama temperature
PData<-read.csv('../Data/PanamaTemperature.csv')
# remove the missing data
PData<-PData[PData$chk_note=='good',]

# make dates
PData$datetime<-parse_date_time(PData$datetime, c('dmy_hms','dmy_hm')) # they are in multiple formats
PData$date<-dmy(PData$date)

# Bermuda Data
BData2016<-read.table('../Data/bepb6h2016.txt')
colnames(BData2016)<-c('YY',  'MM', 'DD', 'hh', 'mm', 'WDIR', 'WSPD', 'GST',  'WVHT',  
               'DPD',   'APD', 'MWD',   'PRES',  'ATMP',  'WTMP',  'DEWP',  'VIS',  'TIDE')
BData2017<-read.table('../Data/bepb6h2017.txt')
colnames(BData2017)<-c('YY',  'MM', 'DD', 'hh', 'mm', 'WDIR', 'WSPD', 'GST',  'WVHT',  
                       'DPD',   'APD', 'MWD',   'PRES',  'ATMP',  'WTMP',  'DEWP',  'VIS',  'TIDE')

BData<-rbind(BData2016,BData2017)
# make dates
BData$Date.Time<-mdy_hm(paste(BData$MM, BData$DD, BData$YY, BData$hh, ':', BData$mm))

# remove the bad temp data
bad<-which(BData$WTMP==999)
BData<-BData[-bad,]

png('../Output/TempTimeSeries.png', width = 5, height = 4.5, units = 'in', res = 300)
par(mar=c(4.1,4.1,3.1,2.1))
plot(PData$datetime[PData$date>"2016-01-01" & PData$date < "2017-12-31" ],
     PData$wt[PData$date>"2016-01-01"& PData$date < "2017-12-31" ], 
     type = 'l', col = 'tomato', xlab = "", ylab = 'Temperature'~degree~C,
     ylim = c(min(BData$WTMP, na.rm=TRUE), max(PData$wt, na.rm=TRUE)))
lines(BData$Date.Time, BData$WTMP, col = 'skyblue')
legend('bottomright',c('Panama','Bermuda'), 
       lty=1, col = c('tomato','skyblue'), bty = 'n')
dev.off()

PData<-PData %>%
  filter(date>"2016-01-01"& date < "2017-12-31") 

PData$MM<-month(PData$datetime, label = FALSE)
PData$Y<-year(PData$datetime)

# make a plot with boxplots by month for 2017
pdf('../Output/TempTimeSeries_boxplots.pdf', width = 8, height = 5)

boxplot(BData$WTMP~BData$MM, col = 'skyblue',
        ylim = c(min(BData$WTMP, na.rm=TRUE), max(PData$wt, na.rm=TRUE)), xaxt = 'n' )
par(new=T)
boxplot(PData$wt~PData$MM, col = 'tomato',xaxt='n',
        ylim = c(min(BData$WTMP, na.rm=TRUE), max(PData$wt, na.rm=TRUE)) )
axis(1, at = 1:12, labels = c('Jan','Feb','Mar','Apr','May','June','Jul','Aug','Sept','Oct','Nov','Dec'))
dev.off()

map.world <- map_data("world")

Atlantic<-map.world %>%
  filter(region == c('Bermuda','Canada', 'USA', 'Mexico', 
                     'Guatemala','Hondurus', 'El Salvador','Belize',
                     'Panama', 'Columbia'))

# make some maps
# shapefiles from http://datapages.com/gis-map-publishing-program/gis-open-files/global-framework/global-heat-flow-database/shapefiles-list
# http://www.diva-gis.org/gdata
aoi_boundary_Bermuda <- st_read(
  "../Data/BMU_adm/BMU_adm0.shp")

aoi_boundary_Panama <- st_read(
  "../Data/PAN_adm/PAN_adm0.shp")

# aoi_boundary_NAO <- st_read(
#   "../Data/North_Atlantic_Ocean_Shapefile/North_Atlantic_Ocean.shp")

berm<-ggplot() + 
  geom_sf(data = aoi_boundary_Bermuda, size = 2, color = "black", fill = "lightblue") + 
  ggtitle("Bermuda") + 
  coord_sf()+
  theme_light()
ggsave('../Output/berm.png', berm, device = 'png', width = 5, height = 5)

pan<-ggplot() + 
  geom_sf(data = aoi_boundary_Panama, size = 2, color = "black", fill = "tomato") + 
  ggtitle("Panama") + 
  coord_sf()+
  theme_light()
ggsave('../Output/pan.png', pan, device = 'png', width = 5, height = 5)

NAO<-ggplot() +
  geom_polygon(data = map.world, aes(x = long, y = lat, group = group))+
  coord_sf(xlim = c(-100,0), ylim = c(0,50))+
  xlab("")+
  ylab("")+
  theme_light()

ggsave('../Output/NAO.pdf', NAO, device = 'pdf', width = 10, height = 5)

# calculate average yearly max for year site

BData %>%
  group_by(YY) %>%
  summarise(max = max(WTMP, na.rm=T)) %>% # take yearly max
  summarise(mean = mean(max)) # take average max

PData %>%
  group_by(Y) %>%
  summarise(max = max(wt, na.rm=T))%>% # take yearly max
  summarise(mean = mean(max)) # take average max
  
