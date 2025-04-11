source('C:/Users/anne.simonis/Documents/GitHub/PAMpal/devel/matchTemplateFunctions.R')
library(PAMpal)


###USER-DEFINED FIELDS####
baseDir <- 'H:/Odontocetes/Beaked whales/'  
DriftName<- paste0('ADRIFT_',101:108)

# change if the set of templates changes
templateNames <- c("ZC","BW43","BW39V","MS","BB","BWC")
# keeping the match/reject values just in case they are useful, but
# can remove these in future if they aren't necessary
extraCols <- c(paste0(templateNames, '_match'))

for(d in 1:length(DriftName)){
  binFolder <- paste0(baseDir,'Binaries/',DriftName[d])    #Folder with binaries
  db <- paste0(baseDir,'Databases/',DriftName[d],' - Copy.sqlite3')
  saveFile <- paste0(baseDir,DriftName[d],'.rds')
  
  # allData <- loadTemplateFolder(binFolder, names=templateNames, extraCols=extraCols, file=saveFile)

  allData<-readRDS(saveFile)  #load saved binary info
  
  # these are in order of "templateNames" above. Can look at data and see if any of these need to
  # be raised/lowered
  threshVals <- c(.55, .45, .6, .5, .25, .6)
  allData <- addTemplateLabels(allData, db=NULL, templateNames, threshVals)
  # minDets is minimum detections to count as an event, maxSep is max time between detections
  # before an event is ended. maxLength is maximum length of an event
  allData <- markGoodEvents(allData, minDets=3, maxSep=120, maxLength=120)
  
  # adds events meeting nDets/nSeconds criteria to the database
  # make sure db is a COPY of the original for safety
  addTemplateEvents(db, binFolder, allData)
  
  ### OPTIONAL process again with PAMpal to do stuff ###
  pps <- PAMpalSettings(db, binFolder, sr_hz=288e3, filterfrom_khz=10, filterto_khz=NULL, winLen_sec=.0025)
  data <- processPgDetections(pps, mode='db', id=paste0('MatchTemp_',DriftName[d]))
  data <- setSpecies(data, method = 'pamguard')
 
  saveRDS(data, paste0(baseDir,'matchTemplateStudy_',DriftName[d],'.rds'))

  }
