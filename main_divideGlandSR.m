%%main_divideGlandSR
clear all 
close all

addpath(genpath('src'))

%%select phenotype
pathKindPhenotype = uigetdir();
pathGlands = dir(fullfile(pathKindPhenotype,'**','\layersTissue.mat'));

%%select excel with extracted features
pathExcelFeatures = uigetfile([pathKindPhenotype '/*xls'],'Get excel of extracted features');
T_features = readtable(fullfile(pathKindPhenotype,pathExcelFeatures));

parfor nGland = 1:size(pathGlands,1)
    
    idGland = cellfun(@(x) contains(pathGlands(nGland).folder,strrep(x,'/','\')),vertcat(T_features.ID_Glands(:)));    
    allImages = load(fullfile(pathGlands.folder(nGland), 'layersTissue_v3.mat'),'labelledImage', 'apicalLayer', 'basalLayer');
    labelledImage = allImages.labelledImage; basalLayer = allImages.basalLayer;apicalLayer = allImages.apicalLayer;
    finalSR = T_features.SurfaceRatio3D(idGland);
    desiredSR = 1.5:0.5:finalSR;
    interpolateImagesBySR(labelledImage, apicalLayer, basalLayer,finalSR,desiredSR,pathGlands.folder(nGland))

end


