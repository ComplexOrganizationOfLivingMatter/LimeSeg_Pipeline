%%main_divideGlandSR
clear all 
close all

addpath(genpath('src'))

%%select phenotype
pathKindPhenotype = uigetdir();
pathGlands = dir(fullfile(pathKindPhenotype,'**','layersTissue.mat'));

%%select excel with extracted features
pathExcelFeatures = uigetfile(fullfile(pathKindPhenotype, '*.xls'),'Get excel of extracted features');
T_features = readtable(fullfile(pathKindPhenotype,pathExcelFeatures));

finalDesiredSR=min(T_features.SurfaceRatio3D_radii);

parpool(10)
parfor nGland = 1:size(pathGlands,1)
    
    idGland = cellfun(@(x) contains(pathGlands(nGland).folder,strrep(x,'/','\')),vertcat(T_features.ID_Glands(:)));    
    allImages = load(fullfile(pathGlands(nGland).folder, 'layersTissue.mat'),'labelledImage', 'apicalLayer', 'basalLayer');
    labelledImage = allImages.labelledImage; basalLayer = allImages.basalLayer;apicalLayer = allImages.apicalLayer;
    if contains(pathGlands(nGland).folder,'Echinoid')
        atypCells=load(fullfile(pathGlands(nGland).folder, 'atypicalCells.mat'),'atypicalCells');
        basalLayer(ismember(basalLayer,atypCells.atypicalCells))=0;
        apicalLayer(ismember(apicalLayer,atypCells.atypicalCells))=0;
    end
    finalSR = T_features.SurfaceRatio3D_radii(idGland);
    desiredSR = 1.5:0.5:finalDesiredSR;
    interpolateImagesBySR(labelledImage, apicalLayer, basalLayer,finalSR,desiredSR,pathGlands(nGland).folder)

end


