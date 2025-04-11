
%% import a binary file and create an output template that can be used for testing
clear
file = 'H:\Odontocetes\Beaked whales\Banter_Training_Data\Binaries\ADRIFT_050\20230313\Click_Detector_Click_Detector_Clicks_20230313_231800.pgdf';
uid = 3882000058; % the UID of the click to slect a template. 
sr=288000; 
%load the clicks
[clicks, fileinfo]= loadPamguardBinaryFile(file);
%select the click with correct UID to be used as a template
template = clicks([clicks.UID]==uid); 
%create a filename for the template
filename = ['clicktemplate_ZC_ADRIFT_050_UID' num2str(uid) '.csv'];
%write template to csv file.
[clickstruct] = clickstruct2csv(template, sr, filename);


function [clickstruct] = clickstruct2csv(clickstruct, sr, filename)
%CLICKSTRUCT2CSV converts a click structure to a template that can be read
%by PAMGuard's matched click classifier.
A = clickstruct.wave(:,1)'; 
format long
dlmwrite(filename, A);
dlmwrite(filename, sr, '-append');
end