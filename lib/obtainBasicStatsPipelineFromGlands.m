function [averageGlandInfo] = obtainBasicStatsPipelineFromGlands(files, minNumberOfSurfaceRatios, colourSample, outputFolderName)
%OBTAINBASICSTATSPIPELINEFROMGLANDS Summary of this function goes here
%   Detailed explanation goes here

nonDiscardedFiles = cellfun(@(x) contains(lower(x), 'discarded') == 0, {files.folder});
files = files(nonDiscardedFiles);

resultsFileName = 'glandDividedInSurfaceRatios_AllUnrollFeatures.mat';

WTMinNumberOfSurfaceRatios = 7;
steps = 2.5/(WTMinNumberOfSurfaceRatios-1);
surfaceRatiosExtrapolatedFrom3D = 1:steps:((steps*(WTMinNumberOfSurfaceRatios-1))+1);
surfaceRatiosExtrapolatedFrom3D = surfaceRatiosExtrapolatedFrom3D(1:minNumberOfSurfaceRatios);

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
    
    if exist(fullfile(files(numFile).folder, 'stats.mat'), 'file') == 0
        %% Calculate final variables
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

        [polygonDistribution_total] = calculate_polygon_distribution(numNeighAccumPerSurfacesRealization(:, end), validCells);

        polygonDistribution_total = cell2table(polygonDistribution_total(2, 1:15), 'VariableNames', strcat('total_', polygonDistribution_total(1, 1:15)));
        
        save(fullfile(files(numFile).folder, 'stats'), 'numNeighPerSurfaceRealization', 'numNeighAccumPerSurfacesRealization', 'numNeighOfNeighPerSurfacesRealization', ...
            'percentageScutoids', 'numNeighOfNeighAccumPerSurfacesRealization', 'areaCellsPerSurfaceRealization', ...
            'volumePerSurfaceRealization', 'polygonDistribution_total', 'apicoBasalTransitions');
    else
        load(fullfile(files(numFile).folder, 'stats'))
    end
    
    namesSR = surfaceRatiosExtrapolatedFrom3D;
    namesSR = arrayfun(@(x) ['sr' strrep(num2str(x),'.','_')], namesSR, 'UniformOutput', false);
    
    %Save Initial information
    numNeighPerSurface{numFile, 1} = array2table(numNeighPerSurfaceRealization(validCells, 1:minNumberOfSurfaceRatios),'VariableNames',namesSR);
    numNeighAccumPerSurfaces{numFile, 1} = array2table(numNeighAccumPerSurfacesRealization(validCells, 1:minNumberOfSurfaceRatios),'VariableNames',namesSR);
    numNeighOfNeighPerSurface{numFile, 1} = array2table(numNeighOfNeighPerSurfacesRealization(validCells, 1:minNumberOfSurfaceRatios),'VariableNames',namesSR);
    numNeighOfNeighAccumPerSurface{numFile, 1} = array2table(numNeighOfNeighAccumPerSurfacesRealization(validCells, 1:minNumberOfSurfaceRatios),'VariableNames',namesSR);
    areaCellsPerSurface{numFile, 1} = array2table(areaCellsPerSurfaceRealization(validCells,1:minNumberOfSurfaceRatios) ./ mean(areaCellsPerSurfaceRealization(validCells,1:minNumberOfSurfaceRatios)),'VariableNames',namesSR);
    volumePerSurface{numFile, 1} = array2table(volumePerSurfaceRealization(validCells,1:minNumberOfSurfaceRatios) ./ mean(volumePerSurfaceRealization(validCells,1:minNumberOfSurfaceRatios)),'VariableNames',namesSR);
    
    xlswrite(fullfile(files(numFile).folder, 'stats.xls'), horzcat(validCells', numNeighPerSurfaceRealization(validCells, :)), 'numNeighPerSR');
    xlswrite(fullfile(files(numFile).folder, 'stats'), horzcat(validCells', numNeighAccumPerSurfacesRealization(validCells, :)), 'numNeighAccumPerSR');
    xlswrite(fullfile(files(numFile).folder, 'stats'), horzcat(validCells', numNeighOfNeighPerSurfacesRealization(validCells, :)), 'numNeighOfNeighPerSR');
    xlswrite(fullfile(files(numFile).folder, 'stats'), horzcat(validCells', numNeighOfNeighAccumPerSurfacesRealization(validCells, :)), 'numNeighOfNeighAccumPerSR');
    xlswrite(fullfile(files(numFile).folder, 'stats'), horzcat(validCells', areaCellsPerSurfaceRealization(validCells, :)), 'areaCellsPerSurfaceRealization');
    xlswrite(fullfile(files(numFile).folder, 'stats'), horzcat(validCells', volumePerSurfaceRealization(validCells, :)), 'volumePerSurface');
    
    %% Additional average info
    meanNumNeighPerSurfaceRealization = mean(numNeighAccumPerSurfacesRealization(validCells, :), 1);
    stdNumNeighPerSurfaceRealization = std(numNeighAccumPerSurfacesRealization(validCells, :), 1);
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
    
    averageGlandInfo{numFile} = [averageGlandInfo{numFile}, polygonDistribution_total];
    
    clearvars 'meanNeighsScutoidsPerSF_ValidCells' 'neighbours'
end
averageGlandInfo = vertcat(averageGlandInfo{:});
save(fullfile(files(numFile).folder, 'averageGlandInfo'), 'averageGlandInfo');
getStatsAndRepresentationsEulerLewis3D(numNeighOfNeighPerSurface,numNeighOfNeighAccumPerSurface,numNeighPerSurface,numNeighAccumPerSurfaces,areaCellsPerSurface,volumePerSurface, fullfile('Results/SalivaryGlands', outputFolderName), surfaceRatiosExtrapolatedFrom3D, colourSample);
end

