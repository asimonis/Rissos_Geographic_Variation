source('C:/Users/anne.simonis/Documents/GitHub/PAMpal/devel/matchTemplateFunctions.R')
library(PAMpal)


###USER-DEFINED FIELDS####
baseDir <- 'H:/Odontocetes/Beaked whales/'  
DriftName<-'PG_2_02_09_CCES_020_288kHz'
binFolder <- paste0(baseDir,'Binaries/','PG_2_02_09_CCES_020_Dec288kHz')    #Folder with binaries
# this database should be a COPY of the original because we will add events to it later
db <- paste0(baseDir,'Databases/',DriftName,' - Copy.sqlite3')
OrigDB<-'H:/Odontocetes/Beaked whales/CCES_2018_Databases/PamGuard64 2_00_16e Drift-15_Final.sqlite3'
###########################


# change if the set of templates changes
templateNames <- c("ZC","BW43","BW39V","MS","BB","BWC")
# keeping the match/reject values just in case they are useful, but
# can remove these in future if they aren't necessary
extraCols <- c(paste0(templateNames, '_match'))

# the binary processing takes a really long time, this automatically saves to an RDS file
# so that you don't have to reprocess in future
saveFile <- paste0(baseDir,DriftName,'.rds')

allData <- loadTemplateFolder(binFolder, names=templateNames, extraCols=extraCols, file=saveFile)

allData<-readRDS(saveFile)  #load saved binary info

# these are in order of "templateNames" above. Can look at data and see if any of these need to
# be raised/lowered
threshVals <- c(.55, .45, .6, .5, .25, .6)
# threshVals <- rep(.01,6)
allData <- addTemplateLabels(allData, db=NULL, templateNames, threshVals)
# minDets is minimum detections to count as an event, maxSep is max time between detections
# before an event is ended. maxLength is maximum length of an event
allData <- markGoodEvents(allData, minDets=3, maxSep=120, maxLength=120)

# summary of how many of the detections in manually annotated events were tagged by template
manualSummary <- summariseManualEvents(allData, db=OrigDB)
# summary of how many detections tagged by template were present in manually annotated events
templateSummary <- summariseTemplateEvents(allData, db=OrigDB)

# Labeling "allData" by the template event results
# Don't want to use all FP template events (these are labeled "none"),
# subset randomly some fixed amount instead
whichNotBW <- which(templateSummary$overlapLabel == 'none')
# good practice to set seed before any randomness for reproducibility
set.seed(16)
# change this to however many NotBW events you want
nSubset <- 200
whichSubset <- sample(whichNotBW, size=nSubset, replace=FALSE)
# we'll only label the selected subset as "NotBW", leave the rest as "none"
templateSummary$overlapLabel[whichSubset] <- 'NotBW'
# attach these template labels to the data
allData <- left_join(allData, templateSummary[c('templateEvent', 'overlapLabel')], by='templateEvent')
# remove all "none" events so we only add the labeled and selected subset to database as events
addTemplateEvents(db, binFolder, filter(allData, overlapLabel != 'none'), labelCol='overlapLabel')


# ###
# allData$templateGood <- allData$templateEvent %in% templateSummary$templateEvent[templateSummary$pctOverlap > 0]
# 
# spd <- summaryPlotData(allData, templateNames, threshVals)        
# ggplot(spd, aes(x=nDets)) +
#   geom_point(aes(y=MeanScore, col=eventLabel), size=2) +
#   facet_wrap(templateGood~Template, ncol=6) +
#   ggtitle('Detections vs Mean Sore')
# spd %>% 
#   select(-Template, -MeanScore) %>% 
#   distinct() %>% 
#   ggplot(aes(x=nDets)) +
#   geom_histogram(aes(fill=eventLabel)) +
#   facet_wrap(~templateGood, scales='free') +
#   xlim(0, 20) +
#   ggtitle('Number of Detections')

tic()
### OPTIONAL process again with PAMpal to do stuff ###
pps <- PAMpalSettings(db, binFolder, sr_hz=288e3, filterfrom_khz=10, filterto_khz=NULL, winLen_sec=.0025)
data <- processPgDetections(pps, mode='db', id='MatchTemp_CCES_012')
data <- setSpecies(data, method = 'pamguard')
# new "FP" and "TP" events in addition to originals
table(species(data))
saveRDS(data, 'H:/Odontocetes/Beaked whales/MTC/matchTemplateStudy_CCES_012_subset_NonBW.rds')
toc()
