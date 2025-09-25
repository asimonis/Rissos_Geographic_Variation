


# 1. Decide which events will be included in analysis
# 2. Read in each acoustic study, clean out the unwanted detections, and then calculate the number of clicks in each event (Parts 1-3 from SpectrogramsAndSpectra code)
# 3. Aggregate all click event data into one data.frame
# 4. Calculate summary statistics (total number of clicks, mean, median, IQR)


ClickSummary <- ClickData %>%
  group_by(eventId)%>%
  summarise(ClickNum = n())

#only do this for 1st acoustic study
AllEvents<-data.frame(ClickSummary)

AllEvents<-rbind(ClickSummary,  AllEvents)

#Save AllEvents file to Github repository
saveRDS(data, paste0(here('data','AllEventSummary.rds'))
