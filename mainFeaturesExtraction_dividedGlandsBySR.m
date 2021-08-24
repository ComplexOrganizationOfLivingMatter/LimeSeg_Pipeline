%% Main features extraction of divided glands for each SR
clear all
close all
addpath(genpath('..\NaturalVariation\Code\'))

%1. Load final segmented glands
pathKindPhenotype = uigetdir();
pathGlands = dir(fullfile(pathKindPhenotype,'**','\layersTissue.mat'));

pathSRs = dir(fullfile(pathGlands(1).folder,'dividedGlandBySr','*mat'));

numOfSRs = size(pathSRs,1);
SRs = 1.5:0.5:(numOfSRs*0.5+1);
for nSR =1:numOfSRs
    %At least the 0.5% of lateral membrane contacting with other cell to be
    %considered as neighbor.
    % contactThreshold = 0.5;
    contactThreshold = 0.1;

    disp(['-----Surface ratio : '  num2str(SRs(nSR)) ' -----'])
    
    allGeneralInfo = cell(size(pathGlands,1),1);
    allTissues = cell(size(pathGlands,1),1);
    allLumens = cell(size(pathGlands,1),1);
    allHollowTissue3dFeatures = cell(size(pathGlands,1),1);
    allNetworkFeatures = cell(size(pathGlands,1),1);
    totalMeanCellsFeatures = cell(size(pathGlands,1),1);
    totalStdCellsFeatures = cell(size(pathGlands,1),1);

    path2saveSummary = [pathKindPhenotype '_' num2str(contactThreshold) '%_SR_' num2str(SRs(nSR)) ' _'];

    % parpool(10)
    realisticSR=zeros(size(pathGlands,1),1);
    parfor nGland = 1:size(pathGlands,1)

            splittedFolder = strsplit(pathGlands(nGland).folder,'\');
            disp([splittedFolder{end-2} '_' splittedFolder{end-1}])
            folderFeatures = [fullfile(pathGlands(nGland).folder,'dividedGlandBySr','Features_vx4_'), num2str(contactThreshold) '_sr' num2str(SRs(nSR))];

            if ~exist(folderFeatures,'dir')
                mkdir(folderFeatures);
            end

            
            if ~exist(fullfile(folderFeatures, 'global_3dFeatures.mat'),'file')
                allImages = load(fullfile(pathGlands(nGland).folder, '\layersTissue.mat'),'apicalLayer','lateralLayer','lumenImage','lumenSkeleton');
                
                dividedImage = load(fullfile(pathGlands(nGland).folder,'dividedGlandBySr',['sr_' num2str(SRs(nSR)) '.mat']),'interpImageCells');
                
                labelledImage = dividedImage.interpImageCells;lumenImage = allImages.lumenImage;lateralLayer = allImages.lateralLayer; apicalLayer = allImages.apicalLayer;lumenSkeleton=allImages.lumenSkeleton;
                lateralLayer(labelledImage==0)=0;    
                cystFilled = imfill(labelledImage>0 | lumenImage>0,'holes');
                perimCystFilled = bwperim(cystFilled);
                basalLayer = uint8(zeros(size(labelledImage)));
                basalLayer(perimCystFilled) = labelledImage(perimCystFilled);

            else
                allImages = load(fullfile(pathGlands(nGland).folder, '\layersTissue.mat'),'apicalLayer','lumenSkeleton');
                labelledImage = []; apicalLayer=allImages.apicalLayer; basalLayer = []; lateralLayer =[]; lumenImage=[]; lumenSkeleton=allImages.lumenSkeleton;
            end
            

            pixelScale=struct2array(load(fullfile(pathGlands(nGland).folder,'pixelScaleOfGland.mat'),'pixelScale'));    
            validNoValidCells = load(fullfile(pathGlands(nGland).folder,'valid_cells.mat'),'validCells','noValidCells'); 
            validCells = validNoValidCells.validCells;
            noValidCells = validNoValidCells.noValidCells;
            fileName = [splittedFolder{end-2} '/' splittedFolder{end-1}];

            [allGeneralInfo{nGland},allTissues{nGland},allLumens{nGland},allHollowTissue3dFeatures{nGland},allNetworkFeatures{nGland},totalMeanCellsFeatures{nGland},totalStdCellsFeatures{nGland}]=calculate3DMorphologicalFeatures(labelledImage,apicalLayer,basalLayer,lateralLayer,lumenImage,folderFeatures,fileName,pixelScale,contactThreshold,validCells,noValidCells);

            realisticSR(nGland) = calculateRealisticSR(totalMeanCellsFeatures{nGland},apicalLayer,lumenSkeleton,validCells,pixelScale,folderFeatures);
    end

    summarizeAllTissuesProperties(allGeneralInfo,allTissues,allLumens,allHollowTissue3dFeatures,allNetworkFeatures,totalMeanCellsFeatures,totalStdCellsFeatures,path2saveSummary);
    
end
