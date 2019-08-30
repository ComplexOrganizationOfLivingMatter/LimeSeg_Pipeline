%% Messy file to create the figures of the article, obtaining several measurements which were not captured directly.
% Use it at your own risk.
clear all
close all

addpath(genpath('src'))
addpath(genpath('lib'))
addpath(genpath(fullfile('..','Epithelia3D', 'InSilicoModels', 'TubularModel', 'src')));

files_WT = dir('**/data/Salivary gland_ExtractedVertices_Correct/**/Results/3d_layers_info.mat');
[averageGlandInfo_WT] = obtainBasicStatsPipelineFromGlands(files_WT, 7, 0, 'WT_30_08_2019');

files_Ecadhi = dir('**/data/Salivary gland/E-cadh Inhibited (flatten)/**/Results/3d_layers_info.mat');

minNumberOfSurfaceRatios = 5; %7WT; 5Ecadhi

%Colour of samples
% WT = 0
% Ecadhi Flatten = 2
colourSample = 2;
outputFolderName = 'EcadhiFlatten/';
[averageGlandInfo_Ecadhi] = obtainBasicStatsPipelineFromGlands(files_Ecadhi, minNumberOfSurfaceRatios, colourSample, outputFolderName);
