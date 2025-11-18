#Need to re-move potential mixed species encounters with Risso's
#Determine which Risso's events overlap with other species,
#or occur within 20 minutes of another odontocete acoustic encounter

library(here)
library(tidyverse)
library(lubridate)


#Load all detection data
data<-read.csv(here('EventTimes','AllDetections_wGPS.csv'))
data$UTC<-as.POSIXct(data$UTC,format='%Y-%m-%d %H:%M:%S',tz='UTC')
data$end<-as.POSIXct(data$end,format='%Y-%m-%d %H:%M:%S',tz='UTC')

#Only keep odontocete events
#Remove all baleen whales and anthro detections, and Risso's
odontocetes<-c('Lo','MS','ZC','UO','BB','MC','BW','BW43','BWC')

OdontoceteData<-data %>%
  filter(species %in% odontocetes)

#Risso's events
GgData<-data %>%
  filter(species=="Gg")

#Loop through all drifts with Risso's
#For each Risso's event
#determine if there is any overlap with other odontocetes
#Generate a list of mixed species events to omit from analysis
Events_MixedSpp<-as.character()
Events_Near_GgStart<-as.character()
Events_Near_GgEnd<-as.character()
  
GgDrifts<-unique(GgData$DriftName)
for(d in 1:length(GgDrifts)){
  #Risso's events in this drift
  GgDriftEvents<-filter(GgData,DriftName==GgDrifts[d])

  #Odontocete events in this drift
  ODriftEvents<-filter(OdontoceteData,DriftName==GgDrifts[d])
  
  #loop through each Risso's event on this drift to look for overlap
  for(e in 1:nrow(GgDriftEvents)){
    overlap<-which(ODriftEvents$UTC>GgDriftEvents$UTC[e] & ODriftEvents$UTC<GgDriftEvents$end[e])
    #look for events within 20 minutes of start
    BeforeGgStart<-which(ODriftEvents$UTC>GgDriftEvents$UTC[e]-minutes(20)& ODriftEvents$UTC<GgDriftEvents$UTC[e])
    AfterGgEnd<-which(ODriftEvents$UTC<GgDriftEvents$end[e]+minutes(20) & ODriftEvents$end>GgDriftEvents$end[e])
  
    if (any(c( overlap, BeforeGgStart, AfterGgEnd) > 0)) {
      if (length(overlap) > 0){
        # print(paste('Event',e,'on Drift:',GgDrifts[d],'is a mixed species event'))
        Events_MixedSpp<-c(Events_MixedSpp,paste0(GgDrifts[d],'_00',e))}
      if (length(BeforeGgStart) > 0){
        # print(paste('Event',e,'on Drift:',GgDrifts[d],'has another odontocete event within 20 min of Gg Event Start'))
        Events_Near_GgStart<-c(Events_Near_GgStart,paste0(GgDrifts[d],'_00',e))}
      if (length(AfterGgEnd) > 0){
        # print(paste('Event',e,'on Drift:',GgDrifts[d],'has another odontocete event within 20 min of Gg Event End'))
        Events_Near_GgEnd<-c(Events_Near_GgEnd,paste0(GgDrifts[d],'_00',e))}
    }
    remove(overlap,BeforeGgStart,AfterGgEnd)
  }
}

EventsToRemove<-c(Events_MixedSpp,Events_Near_GgStart,Events_Near_GgEnd)
sort(EventsToRemove)
write.table(EventsToRemove,file=here('EventTimes','GgEventsToRemove.csv'),
          row.names = F,col.names='EventName')
