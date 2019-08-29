%% Messy file to create the figures of the article, obtaining several measurements which were not captured directly.
% Use it at your own risk.
clear all
close all

addpath(genpath('src'))
addpath(genpath('lib'))
addpath(genpath(fullfile('..','Epithelia3D', 'InSilicoModels', 'TubularModel', 'src')));

files_WT = dir('**/data/Salivary gland_ExtractedVertices_Correct/**/Results/3d_layers_info.mat');
files_Ecadhi = dir('**/data/Salivary gland/E-cadh Inhibited (flatten)/**/Results/3d_layers_info.mat');

files = files_Ecadhi;
nonDiscardedFiles = cellfun(@(x) contains(lower(x), 'discarded') == 0, {files.folder});
files = files(nonDiscardedFiles);

%resultsFileName = '3d_layers_info.mat';
%resultsFileName = 'glandDividedInSurfaceRatios.mat';
%resultsFileName = 'glandDividedInSurfaceRatios_PredefinedSR.mat';
resultsFileName = 'glandDividedInSurfaceRatios_AllUnrollFeatures.mat';

minNumberOfSurfaceRatios = 5; %7WT; 5Ecadhi
%namesSR = arrayfun(@(x) ['sr' strrep(num2str(x),'.','_')],1:numberOfSurfaceRatios,'UniformOutput', false);
%namesSR = {'sr1' 'sr2'};

steps = 2.5/(7-1);
surfaceRatiosExtrapolatedFrom3D = 1:steps:((steps*(7-1))+1);

totalDataBasal = [];
totalDataAccum = [];

realSRBasal = [];
realSRAccum = [];
meanVolumeMicronsPerGland = zeros(length(files),2);



for numFile = 1:length(files)
    files(numFile).folder
    
    load(fullfile(files(numFile).folder, files(numFile).name));
    load(fullfile(files(numFile).folder, 'valid_cells.mat'));
    if exist(fullfile(files(numFile).folder, resultsFileName), 'file')
        load(fullfile(files(numFile).folder, resultsFileName))
    else
        continue
    end
    load([files(numFile).folder '\unrolledGlands\gland_SR_basal\final3DImg.mat'],'img3dComplete')
    
    %% Calculate final variables
    %Volume
    volume3d = regionprops3(img3dComplete,'Volume');
    volume3d = cat(1,volume3d.Volume);
    meanVolumeMicronsPerGland(numFile,1) = mean(volume3d(validCells));
    meanVolumeMicronsPerGland(numFile,2) = std(volume3d(validCells));
    
    %Surface ratio: final 3D, final 2D, selectedSR 3D, selectedSR 2D
    meanSR(numFile, 1) = infoPerSurfaceRatio{end, 2};
    meanSR(numFile, 2) = infoPerSurfaceRatio{end, 7};
    meanSR(numFile, 3) = infoPerSurfaceRatio{minNumberOfSurfaceRatios, 2};
    meanSR(numFile, 4) = infoPerSurfaceRatio{minNumberOfSurfaceRatios, 7};
    
    %Calculate initial parameters
    numberOfSurfaceRatios = size(infoPerSurfaceRatio, 1);
    [areaCellsPerSurfaceRealization, volumePerSurfaceRealization, neighsSurface, neighsAccumSurfaces, percentageScutoids, apicoBasalTransitions, numLostNeighsAccum, numWonNeighsAccum] = getBasicInformationToPlot(infoPerSurfaceRatio, neighboursOfAllSurfaces, numberOfSurfaceRatios);
    
    [numNeighOfNeighPerSurfacesRealization, numNeighPerSurfaceRealization] = getNumNeighsOfNeighs(neighsSurface, numberOfSurfaceRatios);
    [numNeighOfNeighAccumPerSurfacesRealization, numNeighAccumPerSurfacesRealization] = getNumNeighsOfNeighs(neighsAccumSurfaces, numberOfSurfaceRatios);
    
    %namesSR = surfaceRatioOfGland(1:minNumberOfSurfaceRatios);
    namesSR = infoPerSurfaceRatio.SR3D;
    namesSR = arrayfun(@(x) ['sr' strrep(num2str(x),'.','_')], namesSR, 'UniformOutput', false);
    namesSR = namesSR(1:minNumberOfSurfaceRatios);
    
    %Save Initial information
    numNeighPerSurface{numFile, 1} = array2table(numNeighPerSurfaceRealization(validCells, 1:minNumberOfSurfaceRatios),'VariableNames',namesSR);
    numNeighAccumPerSurfaces{numFile, 1} = array2table(numNeighAccumPerSurfacesRealization(validCells, 1:minNumberOfSurfaceRatios),'VariableNames',namesSR);
    numNeighOfNeighPerSurface{numFile, 1} = array2table(numNeighOfNeighPerSurfacesRealization(validCells, 1:minNumberOfSurfaceRatios),'VariableNames',namesSR);
    numNeighOfNeighAccumPerSurface{numFile, 1} = array2table(numNeighOfNeighAccumPerSurfacesRealization(validCells, 1:minNumberOfSurfaceRatios),'VariableNames',namesSR);
    areaCellsPerSurface{numFile, 1} = array2table(areaCellsPerSurfaceRealization(validCells,1:minNumberOfSurfaceRatios) ./ mean(areaCellsPerSurfaceRealization(validCells,1:minNumberOfSurfaceRatios)),'VariableNames',namesSR);
    volumePerSurface{numFile, 1} = array2table(volumePerSurfaceRealization(validCells,1:minNumberOfSurfaceRatios) ./ mean(volumePerSurfaceRealization(validCells,1:minNumberOfSurfaceRatios)),'VariableNames',namesSR);
    
    xlswrite(fullfile(files(numFile).folder, 'stats'), horzcat(validCells', numNeighPerSurfaceRealization(validCells, :)), 'numNeighPerSR');
    xlswrite(fullfile(files(numFile).folder, 'stats'), horzcat(validCells', numNeighAccumPerSurfacesRealization(validCells, :)), 'numNeighAccumPerSR');
    xlswrite(fullfile(files(numFile).folder, 'stats'), horzcat(validCells', numNeighOfNeighPerSurfacesRealization(validCells, :)), 'numNeighOfNeighPerSR');
    xlswrite(fullfile(files(numFile).folder, 'stats'), horzcat(validCells', numNeighOfNeighAccumPerSurfacesRealization(validCells, :)), 'numNeighOfNeighAccumPerSR');
    xlswrite(fullfile(files(numFile).folder, 'stats'), horzcat(validCells', areaCellsPerSurfaceRealization(validCells, :) ./ mean(areaCellsPerSurfaceRealization(validCells, :))), 'volumePerSurface');
    xlswrite(fullfile(files(numFile).folder, 'stats'), horzcat(validCells', volumePerSurfaceRealization(validCells, :) ./ mean(volumePerSurfaceRealization(validCells, :))), 'volumePerSurface');
    
    %% Additional average info
    meanNumNeighPerSurfaceRealization = mean(numNeighAccumPerSurfacesRealization(validCells, :), 1);
    stdNumNeighPerSurfaceRealization = std(numNeighAccumPerSurfacesRealization(validCells, :), 1);
    totalAreaPerSR = sum(areaCellsPerSurfaceRealization(validCells, :)) / sum(areaCellsPerSurfaceRealization(validCells, end));
    mean_PercScutoids = mean(percentageScutoids(validCells, :), 1);
    std_PercScutoids = std(percentageScutoids(validCells, :), 1);
    mean_apicoBasalTransitions = [0, mean(apicoBasalTransitions(validCells, :), 1)];
    std_apicoBasalTransitions = [0, std(apicoBasalTransitions(validCells, :), 1)];
    
    %Scutoids per number of sides
    [meanWinningPerSidePerFile{numFile, 1}, cellsPerSide{numFile}] = calculateMeanWinning3DNeighbours(numNeighAccumPerSurfacesRealization(:, 1:minNumberOfSurfaceRatios), validCells, minNumberOfSurfaceRatios);
    
    averageGlandInfo{numFile} = table(infoPerSurfaceRatio.SR3D(end), infoPerSurfaceRatio.SR2D(end), length(validCells), meanNumNeighPerSurfaceRealization(1), stdNumNeighPerSurfaceRealization(1), ...
        mean(numNeighPerSurfaceRealization(validCells, end)), std(numNeighPerSurfaceRealization(validCells, end)), meanNumNeighPerSurfaceRealization(end), stdNumNeighPerSurfaceRealization(end), ...
        mean_PercScutoids(end), std_PercScutoids(end), mean_apicoBasalTransitions(end), std_apicoBasalTransitions(end), ...
        'VariableNames', {'SurfaceRatio_3D', 'SurfaceRatio_2D', 'nCells', 'neighbours_apical_avg', 'neighbours_apical_std', 'neighbours_basal_avg', 'neighbours_basal_std', 'neighbours_total_avg', 'neighbours_total_std', ...
        'scutoids_avg', 'scutoids_std', 'apicoBasalTransitions_avg', 'apicoBasalTransitions_std'});
end

getStatsAndRepresentationsEulerLewis3D(numNeighOfNeighPerSurface,numNeighOfNeighAccumPerSurface,numNeighPerSurface,numNeighAccumPerSurfaces,areaCellsPerSurface,volumePerSurface,'Results\SalivaryGlands\', surfaceRatiosExtrapolatedFrom3D, 0);