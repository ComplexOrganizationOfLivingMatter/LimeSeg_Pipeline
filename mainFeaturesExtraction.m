%%main features extraction segmented Cysts
clear all
close all
addpath(genpath('..\..\NaturalVariation\Code\'))

%1. Load final segmented glands
pathKindPhenotype = uigetdir();
pathGlands = dir(fullfile(pathKindPhenotype,'**','\3d_layers_info.mat'));

%At least the 0.5% of lateral membrane contacting with other cell to be
%considered as neighbor.
contactThreshold = 0.5;


allGeneralInfo = cell(size(pathGlands,1),1);
allTissues = cell(size(pathGlands,1),1);
allLumens = cell(size(pathGlands,1),1);
allHollowTissue3dFeatures = cell(size(pathGlands,1),1);
allNetworkFeatures = cell(size(pathGlands,1),1);
totalMeanCellsFeatures = cell(size(pathGlands,1),1);
totalStdCellsFeatures = cell(size(pathGlands,1),1);

path2saveSummary = [pathKindPhenotype '_' num2str(contactThreshold) '%_'];

for nGland = 1:size(pathGlands,1)
        
        splittedFolder = strsplit(pathGlands(nGland).folder,'\');
        display([splittedFolder{end-2} '_' splittedFolder{end-1}])
        folderFeatures = [fullfile(pathGlands(nGland).folder,'Features'), num2str(contactThreshold)];
        if ~exist(folderFeatures,'dir')
            mkdir(folderFeatures);
        end

        if ~exist(fullfile(pathGlands(nGland).folder, '\layersTissue.mat'),'file')
            if exist(fullfile(pathGlands(nGland).folder,'realSize3dLayers.mat'),'file')
                load(fullfile(pathGlands(nGland).folder,'realSize3dLayers.mat'),'labelledImage_realSize')
                labelledImage = labelledImage_realSize;
                clearvars labelledImage_realSize
            else
                load(fullfile(pathGlands(nGland).folder,pathGlands(nGland).name),'labelledImage')
                load(fullfile(pathGlands(nGland).folder,'zScaleOfGland.mat'),'zScale')

                labelledImage = imresize3(labelledImage,[size(labelledImage,1),size(labelledImage,2),round(size(labelledImage,3)*zScale)],'nearest');
                if size(labelledImage,3)>size(labelledImage,1)
                    labelledImage = imresize3(labelledImage,[size(labelledImage,1)*zScale,size(labelledImage,2)*zScale,round(size(labelledImage,3))],'nearest');
                end

            end

            load(fullfile(pathGlands(nGland).folder,'pixelScaleOfGland.mat'),'pixelScale')    
            fileName = [splittedFolder{end-2} '/' splittedFolder{end-1}];
            
            %%get apical and basal layers, and Lumen
            [apicalLayer,basalLayer,lateralLayer,lumenImage] = getApicalBasalLateralAndLumenFromCyst(labelledImage);
            save(fullfile(pathGlands(nGland).folder, '\layersTissue.mat'),'apicalLayer','basalLayer','lateralLayer','lumenImage','labelledImage','-v7.3')
        else
            if ~exist(fullfile(folderFeatures, 'global_3dFeatures.mat'),'file')
                load(fullfile(pathGlands(nGland).folder, '\layersTissue.mat'),'apicalLayer','basalLayer','lateralLayer','lumenImage','labelledImage')
            else
                labelledImage = []; apicalLayer=[]; basalLayer = []; lateralLayer =[]; lumenImage=[];
            end
        end

        [allGeneralInfo{nGland},allTissues{nGland},allLumens{nGland},allHollowTissue3dFeatures{nGland},allNetworkFeatures{nGland},totalMeanCellsFeatures{nGland},totalStdCellsFeatures{nGland}]=calculate3DMorphologicalFeatures(labelledImage,apicalLayer,basalLayer,lateralLayer,lumenImage,folderFeatures,fileName,pixelScale,contactThreshold);
end

summarizeAllTissuesProperties(allGeneralInfo,allTissues,allLumens,allHollowTissue3dFeatures,allNetworkFeatures,totalMeanCellsFeatures,totalStdCellsFeatures,path2saveSummary);
