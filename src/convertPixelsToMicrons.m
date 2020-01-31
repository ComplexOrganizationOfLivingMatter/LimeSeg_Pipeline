function [meanFeatures,stdFeatures, gland3dFeatures, lumen3dFeatures,hollowGland3dFeatures] = convertPixelsToMicrons(files,numFile,meanFeatures,stdFeatures, gland3dFeatures, lumen3dFeatures,hollowGland3dFeatures)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
volumeSubstring={'Volume','volume'};
areaSubstring={'Area', 'area'};
lengthSubstring={'length','Length'};
fileName = strsplit(files(numFile).folder, {'/','\'});
fileName = convertCharsToStrings(strjoin({fileName{1,end-2},fileName{1,end-1}}, ' '));

load(fullfile(files(numFile).folder, 'zScaleOfGland'))

if exist(fullfile(files(numFile).folder, 'pixelScaleOfGland.mat'), 'file') == 0
        pixelScale = inputdlg(strcat('Insert pixel width of Gland',{'  '},'"',fileName,'"'));
        pixelScale = str2double(pixelScale{1});

        save(fullfile(files(numFile).folder, 'pixelScaleOfGland.mat'), 'pixelScale');
else
        load(fullfile(files(numFile).folder, 'pixelScaleOfGland.mat')); 
end
    

%% VolumeFeatures
CellsFeaturesVolumeIndexs = contains(meanFeatures.Properties.VariableNames,volumeSubstring);
TissueFeaturesVolumeIndexs = contains(gland3dFeatures.Properties.VariableNames,"Volume");

meanFeatures(:,CellsFeaturesVolumeIndexs) = splitvars(table(table2array(meanFeatures(:,CellsFeaturesVolumeIndexs)) * pixelScale^3),1);
stdFeatures(:,CellsFeaturesVolumeIndexs) = splitvars(table(table2array(stdFeatures(:,CellsFeaturesVolumeIndexs)) * pixelScale^3),1);
gland3dFeatures(:,TissueFeaturesVolumeIndexs) = splitvars(table(table2array(gland3dFeatures(:,TissueFeaturesVolumeIndexs)) * pixelScale^3),1);
lumen3dFeatures(:,TissueFeaturesVolumeIndexs) = splitvars(table(table2array(lumen3dFeatures(:,TissueFeaturesVolumeIndexs)) * pixelScale^3),1);
hollowGland3dFeatures(:,TissueFeaturesVolumeIndexs) = splitvars(table(table2array(hollowGland3dFeatures(:,TissueFeaturesVolumeIndexs)) * pixelScale^3),1);

%% Area Features
CellsFeaturesAreaIndexs = contains(meanFeatures.Properties.VariableNames,areaSubstring);
TissueFeaturesAreaIndexs = contains(gland3dFeatures.Properties.VariableNames,areaSubstring);

meanFeatures(:,CellsFeaturesAreaIndexs) = splitvars(table(table2array(meanFeatures(:,CellsFeaturesAreaIndexs)) * pixelScale^2),1);
stdFeatures(:,CellsFeaturesAreaIndexs) = splitvars(table(table2array(stdFeatures(:,CellsFeaturesAreaIndexs)) * pixelScale^2),1);
gland3dFeatures(:,TissueFeaturesAreaIndexs) = splitvars(table(table2array(gland3dFeatures(:,TissueFeaturesAreaIndexs)) * pixelScale^2),1);
lumen3dFeatures(:,TissueFeaturesAreaIndexs) = splitvars(table(table2array(lumen3dFeatures(:,TissueFeaturesAreaIndexs)) * pixelScale^2),1);
hollowGland3dFeatures(:,TissueFeaturesAreaIndexs) = splitvars(table(table2array(hollowGland3dFeatures(:,TissueFeaturesAreaIndexs)) * pixelScale^2),1);

%% Length Features
CellsFeaturesLengthIndexs = contains(meanFeatures.Properties.VariableNames,lengthSubstring);
TissueFeaturesLengthIndexs = contains(gland3dFeatures.Properties.VariableNames,lengthSubstring);

meanFeatures(:,CellsFeaturesLengthIndexs) = table(table2array(meanFeatures(:,CellsFeaturesLengthIndexs)) * pixelScale);
stdFeatures(:,CellsFeaturesLengthIndexs) = table(table2array(stdFeatures(:,CellsFeaturesLengthIndexs)) * pixelScale);
gland3dFeatures(:,TissueFeaturesLengthIndexs) = table(table2array(gland3dFeatures(:,TissueFeaturesLengthIndexs)) * pixelScale);
lumen3dFeatures(:,TissueFeaturesLengthIndexs) = table(table2array(lumen3dFeatures(:,TissueFeaturesLengthIndexs)) * pixelScale);
hollowGland3dFeatures(:,TissueFeaturesLengthIndexs) = table(table2array(hollowGland3dFeatures(:,TissueFeaturesLengthIndexs)) * pixelScale);

end

