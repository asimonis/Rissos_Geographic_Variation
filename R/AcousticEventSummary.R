


# 1. Decide which events will be included in analysis
# 2. Read in each acoustic study, clean out the unwanted detections, and then calculate the number of clicks in each event (Parts 1-3 from SpectrogramsAndSpectra code)
# 3. Aggregate all click event data into one data.frame
# 4. Calculate summary statistics (total number of clicks, mean, median, IQR)




#Start by loading the required packages

library("easypackages")
library("lubridate")
library("ggplot2")
library("PAMpal")
library("PAMmisc")
library("dplyr")
library("here")
library("beepr")
here()

#Define details of deployment to process
  
###USER-DEFINED FIELDS#### 
baseDir <- 'E:/Analysis/'
DriftID<-'PASCAL_012'
binFolder <- paste0(baseDir,'Binaries/',DriftID)    #Folder with binaries
# this database should be a COPY of the original because we will add events to it later
db <- paste0(baseDir,'Databases/',DriftID,' - Clean.sqlite3')


#Load in acoustic study from an RDS file and Gg times from log

data<-readRDS(paste0(baseDir,'AcousticStudies/',DriftID,'_Gg.rds'))

# Double check warning messages
print(getWarnings(data)$message)

#Import known Risso's dolphin Event Times
GgTimes<-read.csv(here('EventTimes','ADRIFT_GgDetections.csv'))
GgTimes$start<-as.POSIXct(GgTimes$UTC,format='%Y-%m-%d %H:%M:%S',tz='UTC')
GgTimes$end<-as.POSIXct(GgTimes$end,format='%Y-%m-%d %H:%M:%S',tz='UTC')
DriftTimes<-filter(GgTimes,DriftName==DriftID)
DriftTimes$id<-paste0(DriftID,'_',seq(1:nrow(DriftTimes)))

if(substr(DriftID,1,4)=='CCES'){
  DriftTimes$end<-DriftTimes$start + seconds(119) }

#Clean click detections

#Keep click detections from one channel (upper hydrophone = HTI-92-WB)
data<-filter(data, Channel==1)

#Omit click detectors 0 and 1
data<-filter(data,detectorName!='Click_Detector_0')
data<-filter(data,detectorName!='Click_Detector_1')
data<-filter(data, peak>15 & duration<2000)
ClickData<-getClickData(data)

ClickSummary <- ClickData %>%
  group_by(eventId)%>%
  summarise(ClickNum = n())





#only do this for 1st acoustic study
  #All_Useable_Events<-data.frame(ClickSummary)


  
All_Useable_Events<-rbind(ClickSummary,  All_Useable_Events)

#Save All_Useable_Events file to Github repository
  saveRDS(All_Useable_Events, paste0(here('data'), '/All_Useable_Events_Summary.rds'))

  
install.packages("writexl")
library(writexl)        
AllEvents <- readRDS('C:/Users/sarah/OneDrive/Documents/GitHub/Rissos_Geographic_Variation/Data/All_Useable_Events_Summary.rds')
write_xlsx(AllEvents, 'C:/Users/sarah/OneDrive/Documents/GitHub/Rissos_Geographic_Variation/data/All_Useable_Events_Summary.xlsx')



#To create histogram
Data <- readxl::read_xlsx('C:/Users/sarah/OneDrive/Documents/GitHub/Rissos_Geographic_Variation/data/All_Useable_Events_Summary.xlsx')

Data$species <- as.factor("Gg")

ggplot(Data,aes(species, ClickNum))+geom_violin()

ggplot(Data,aes(ClickNum))+geom_histogram(bins = 100)+ggtitle("Distribution of Clicks Per Useable Event")



















