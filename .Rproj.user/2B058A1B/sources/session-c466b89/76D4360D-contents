
# make sure you have Rtools installed
if(!require('devtools')) install.packages('devtools')
# install from GitHub
devtools::install_github('TaikiSan21/PAMpal')
# #PAMpal tutorial available here: https://taikisan21.github.io/PAMpal/

source(here('R','Matched Click Template Detector','matchTemplateFunctions.R'))

#load packages
library(here)
library(PAMpal)
library(dplyr)

###USER-DEFINED FIELDS####
baseDir <- 'D:/Analysis/RissosGeographicVariation/'  
DriftName<-'ADRIFT_032'
binFolder <- paste0(baseDir,'Binaries/',DriftName)    #Folder with binaries
# this database should be a COPY of the original because we will add events to it later
db <- paste0(baseDir,'Databases/',DriftName,' - Copy.sqlite3')
###########################

### Process with PAMpal ###
pps <- PAMpalSettings(db, binFolder, sr_hz='auto', filterfrom_khz=10, filterto_khz=80, winLen_sec=.0025)

#Use event times from CSV file
data <- processPgDetections(pps, mode='time',grouping=here('EventTimes','ADRIFT_032_GgTimes.csv'),
                            format='%m/%d/%y %H:%M',species='Gg',id=DriftName)

#Keep click detections from one channel
data<-filter(data, Channel==2)
#Omit click detectors 0 and 1
data<-filter(data,detectorName!='Click_Detector_0')
data<-filter(data,detectorName!='Click_Detector_1')

#Remove likely false positives based on following features:
#peak frequency<19 kHz
ClickData<-getClickData(data)
ClickData<-filter(ClickData, peak>19)

#Add events to database
addPgEvent(db,binary = data@files[["binaries"]],UIDs=ClickData$UID,eventType='GgClick')

#save
saveRDS(data, paste0(baseDir,'AcousticStudies/',DriftName,'_Gg.rds'))