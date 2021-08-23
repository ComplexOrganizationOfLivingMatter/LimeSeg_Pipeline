clear all
close all
addpath(genpath('..\NaturalVariation\Code\'))

%1. Load final segmented glands
pathKindPhenotype = uigetdir();
pathGlands = dir(fullfile(pathKindPhenotype,'**','\layersTissue_v2.mat'));

parfor nGland = 1:size(pathGlands,1)

    allImages = load(fullfile(pathGlands(nGland).folder, '\layersTissue_v2.mat'),'lateralLayer','lumenImage','labelledImage');
    labelledImage = allImages.labelledImage;lumenImage = allImages.lumenImage;lateralLayer = allImages.lateralLayer;

    flattenImage = convertLabelledImage2Flatten(labelledImage, lateralLayer);

%     volumeSegmenter(labelledImage,flattenImage);
    
    path2saveLayers = fullfile(pathGlands(nGland).folder, '\layersTissue_flatten.mat');
    [apicalLayer,basalLayer,lateralLayer,lumenSkeleton] = getApicalBasalLateralFromGlands(flattenImage,lumenImage,path2saveLayers);

end