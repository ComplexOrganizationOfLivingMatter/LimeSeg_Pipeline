addpath(genpath('src'))
addpath(genpath('lib'))
addpath(genpath('gui'))

close all
clear all

files = dir('**/data/Salivary gland/**/global_3DFeatures.mat');

allFilesName = [];
for numFiles=1:length(files)
    
load(fullfile(files(numFiles).folder, 'global_3DFeatures.mat'), 'totalMeanFeatures'); 

allHypothesis =[];
allPvalue = [];

for nFeatures=1:size(totalMeanFeatures,2)
    [hypothesis, pValue] = shapiroWilkTest(totalMeanFeatures{:,nFeatures});
    allHypothesis= [allHypothesis, hypothesis];   
    allPvalue = [allPvalue, pValue];    
end

fileName = strsplit(files(numFiles).folder, {'/','\'});
fileName = convertCharsToStrings(fileName{end});
allFilesName = [allFilesName ; fileName];

allShapiroWilkTests=allFilesName; 
    end