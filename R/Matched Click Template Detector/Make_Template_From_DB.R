#Make click template from PG Database using PAMpal

##########################
library(PAMmisc)
library(PAMpal)   #note, this needs v. 0.14 or higher of PAMpal
library(tuneR)
library(signal)

library(devtools)
library(remotes)
devtools::install_github('TaikiSan21/PamBinaries')
remotes::install_github("TaikiSan21/PAMmisc")
devtools::install_github('TaikiSan21/PAMPal')

############################################################################
# Input paths to PAMGuard files, sound files, etc. 
############################################################################

ClickDuration<-0.4
TemplDuration=0.8

# DB and Binaries for processing
db <- 'H:/Odontocetes/NBHF/TrainingData/PG2_02_09_CCES_022_Ksp.sqlite3'
bin <-'H:/Odontocetes/NBHF/TrainingData/PG2_02_09_CCES_022_Ksp'
# 
# db <- 'H:/Odontocetes/Beaked whales/MTC/ADRIFT_050 - Copy.sqlite3'
# bin <-'H:/Odontocetes/Beaked whales/MTC/ADRIFT_050 - Copy'


# Set path where wav clips will be stored 
WavClipFolder <- 'C:/Users/anne.simonis/Documents/GitHub/Odontocetes/data/Species_template_csvs'

# Run PAMpal (any settings here can be changed to your desired values)
Fs=384000
pps <- PAMpalSettings(db=db, binaries=bin, sr_hz=Fs, filterfrom_khz=100, filterto_khz=NULL, winLen_sec=.0025)

# Process PAMGuard data into the Acoustic Study Object 
data <- processPgDetections(pps, mode='db', id="Ksp")

# data@events<-data@events[19]

#Use PAMpal to extract click wave
# uids <- unique(data@events[["PG2_02_09_CCES_022_Ksp.OE25"]]@detectors[["Click_Detector_5"]][["UID"]])
uids <-unique(data@events[[1]]@detectors[["Click_Detector_5"]][["UID"]])
binData <- getBinaryData(data, UID=uids)
Click<-binData[[1]]$wave[,1]

#normalize
seewave::savewav(Click,Fs,channel=1,filename='C:/Users/anne.simonis/Documents/GitHub/Odontocetes/data/Species_template_csvs/CCES_022_Ksp_Click.wav')
wavclip<-readWave('C:/Users/anne.simonis/Documents/GitHub/Odontocetes/data/Species_template_csvs/CCES_022_Ksp_Click.wav')

# savewav(Click,Fs,channel=1,filename='C:/Users/anne.simonis/Documents/GitHub/Odontocetes/data/Species_template_csvs/ADRIFT_050_Zc.wav')
# wavclip<-readWave('C:/Users/anne.simonis/Documents/GitHub/Odontocetes/data/Species_template_csvs/ADRIFT_050_Zc.wav')

ClickNorm<-normalize(wavclip,center=TRUE,rescale=TRUE)

#Reduce to smaller clip
SampleDuration<-length(ClickNorm)/Fs
peak<-which.max(abs(ClickNorm@left))
ClickSamples<- Fs*ClickDuration/1000
SnipClip<-ClickNorm@left[(peak-ClickSamples/2):(peak+ClickSamples/2)]

SnipClip<- 0.25*SnipClip
# add zero padding to both ends of sample for a sample of length Fs*TemplDuration/1000
n= length(SnipClip)
add= ((Fs*TemplDuration/1000) - n)/2

# plot snip around maximum value 
plot(SnipClip,type='l')

# create first row of output csv from zero-padded data
row1= c(rep(0,add),SnipClip,rep(0,add))
if (length(row1) != (Fs*TemplDuration/1000)) row1= c(row1,0)   #add another zero if zero-padding should have been an odd number
plot(row1,type="l")

# second row of output csv includes only sample rate and duration
# row2= array(NA,(Fs*TemplDuration/1000))
row2= array(NA,length(row1))
row2[1]= Fs
# merge two rows and export as csv
output= rbind(row1,row2)
write.table(output,file='C:/Users/anne.simonis/Documents/GitHub/Odontocetes/data/Species_template_csvs/Kogia_384kHz.csv',col.names=F,row.names=F,sep=',')
