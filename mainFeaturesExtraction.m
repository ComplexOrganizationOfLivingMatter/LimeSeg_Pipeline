%%main features extraction glands
clear all
close all
addpath(genpath('..\NaturalVariation\Code\'))

%1. Load final segmented glands
pathKindPhenotype = uigetdir();
pathGlands = dir(fullfile(pathKindPhenotype,'**','\3d_layers_info.mat'));

%At least the 0.5% of lateral membrane contacting with other cell to be
%considered as neighbor.
% contactThreshold = 0.5;
contactThreshold = 0.1;

allGeneralInfo = cell(size(pathGlands,1),1);
allTissues = cell(size(pathGlands,1),1);
allLumens = cell(size(pathGlands,1),1);
allHollowTissue3dFeatures = cell(size(pathGlands,1),1);
allNetworkFeatures = cell(size(pathGlands,1),1);
totalMeanCellsFeatures = cell(size(pathGlands,1),1);
totalStdCellsFeatures = cell(size(pathGlands,1),1);

path2saveSummary = [pathKindPhenotype '_' num2str(contactThreshold) '%_'];

parpool(10)
realisticSR=zeros(size(pathGlands,1),1);
parfor nGland = 1:size(pathGlands,1)
        
        splittedFolder = strsplit(pathGlands(nGland).folder,'\');
        display([splittedFolder{end-2} '_' splittedFolder{end-1}])
        folderFeatures = [fullfile(pathGlands(nGland).folder,'Features_vx4_'), num2str(contactThreshold)];
        
        if ~exist(folderFeatures,'dir')
            mkdir(folderFeatures);
        end

        if ~exist(fullfile(pathGlands(nGland).folder, '\layersTissue_v2.mat'),'file')
            if exist(fullfile(pathGlands(nGland).folder,'realSize3dLayers.mat'),'file')
                realSizeImages = load(fullfile(pathGlands(nGland).folder,'realSize3dLayers.mat'),'labelledImage_realSize','lumenImage_realSize');
                labelledImage = realSizeImages.labelledImage_realSize;
                lumenImage = realSizeImages.lumenImage_realSize;
            else
                images=load(fullfile(pathGlands(nGland).folder,pathGlands(nGland).name),'labelledImage','lumenImage');
                zScale=struct2array(load(fullfile(pathGlands(nGland).folder,'zScaleOfGland.mat'),'zScale'));
                labelledImage = images.labelledImage;
                lumenImage = images.lumenImage;
                labelledImage = imresize3(labelledImage,[size(labelledImage,1),size(labelledImage,2),round(size(labelledImage,3)*zScale)],'nearest');
                lumenImage = imresize3(lumenImage,[size(lumenImage,1),size(lumenImage,2),round(size(lumenImage,3)*zScale)],'nearest');
                if size(labelledImage,3)/size(labelledImage,1)>0.5
                    labelledImage = imresize3(labelledImage,[size(labelledImage,1)*zScale,size(labelledImage,2)*zScale,round(size(labelledImage,3))],'nearest');
                    lumenImage = imresize3(lumenImage,[size(lumenImage,1)*zScale,size(lumenImage,2)*zScale,round(size(lumenImage,3))],'nearest');
                end

            end
            
            %%get apical and basal layers, and Lumen
            path2saveLayers = fullfile(pathGlands(nGland).folder, '\layersTissue_v2.mat');
            [apicalLayer,basalLayer,lateralLayer,lumenSkeleton] = getApicalBasalLateralFromGlands(labelledImage,lumenImage,path2saveLayers);
            
        else
            if ~exist(fullfile(folderFeatures, 'global_3dFeatures.mat'),'file')
                allImages = load(fullfile(pathGlands(nGland).folder, '\layersTissue_v2.mat'),'apicalLayer','basalLayer','lateralLayer','lumenImage','labelledImage','lumenSkeleton');
                labelledImage = allImages.labelledImage;lumenImage = allImages.lumenImage;lateralLayer = allImages.lateralLayer; basalLayer = allImages.basalLayer;apicalLayer = allImages.apicalLayer;lumenSkeleton=allImages.lumenSkeleton;
                
            else
                allImages = load(fullfile(pathGlands(nGland).folder, '\layersTissue_v2.mat'),'apicalLayer','lumenSkeleton');
                labelledImage = []; apicalLayer=allImages.apicalLayer; basalLayer = []; lateralLayer =[]; lumenImage=[]; lumenSkeleton=allImages.lumenSkeleton;
            end
        end
        
        pixelScale=struct2array(load(fullfile(pathGlands(nGland).folder,'pixelScaleOfGland.mat'),'pixelScale'));    
        validNoValidCells = load(fullfile(pathGlands(nGland).folder,'valid_cells.mat'),'validCells','noValidCells'); 
        validCells = validNoValidCells.validCells;
        noValidCells = validNoValidCells.noValidCells;
        fileName = [splittedFolder{end-2} '/' splittedFolder{end-1}];
        
        [allGeneralInfo{nGland},allTissues{nGland},allLumens{nGland},allHollowTissue3dFeatures{nGland},allNetworkFeatures{nGland},totalMeanCellsFeatures{nGland},totalStdCellsFeatures{nGland}]=calculate3DMorphologicalFeatures(labelledImage,apicalLayer,basalLayer,lateralLayer,lumenImage,folderFeatures,fileName,pixelScale,contactThreshold,validCells,noValidCells);
        
        realisticSR(nGland) = calculateRealisticSR(totalMeanCellsFeatures{nGland},apicalLayer,lumenSkeleton,validCells,pixelScale);
end

summarizeAllTissuesProperties(allGeneralInfo,allTissues,allLumens,allHollowTissue3dFeatures,allNetworkFeatures,totalMeanCellsFeatures,totalStdCellsFeatures,path2saveSummary);
