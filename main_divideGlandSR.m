%%main_divideGlandSR
clear all 
close all

addpath(genpath('src'))

path2save = 'E:\Pedro\LimeSeg_Pipeline\data\Salivary gland\Wildtype\2017-12-04\1a\Results\';
load(fullfile(path2save, 'layersTissue_v3.mat'),'labelledImage', 'apicalLayer', 'basalLayer')

finalSR = 5.28815299623935;
desiredSR = 1.5:0.5:finalSR;
interpolateImagesBySR(labelledImage, apicalLayer, basalLayer,finalSR,desiredSR,path2save)