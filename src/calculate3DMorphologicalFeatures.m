function [] = calculate3DMorphologicalFeatures(folderName)
%CALCULATE3DMORPHOLOGICALFEATURES Summary of this function goes here
%   Detailed explanation goes here

if contains(folderName, '\data\Salivary gland\')
    files = dir(fullfile(folderName, '/3d_layers_info.mat'));
else
    files = dir(fullfile('**/data/Salivary gland/', folderName, '/**/Results/3d_layers_info.mat'));
end
nonDiscardedFiles = cellfun(@(x) contains(lower(x), 'discarded') == 0, {files.folder});
files = files(nonDiscardedFiles);

selpath = dir(fullfile('**/data/Salivary gland/', folderName));

totalMeanFeatures = cell([length(files) 16]);
totalStdFeatures = cell([length(files) 16]);
allGlands = cell([length(files) 41]);
allLumens = cell([length(files) 17]);
allGeneralInfo = cell(length(files), 4);
totalSTD3DNeighsFeatures = cell(length(files), 6);
totalMean3DNeighsFeatures = cell(length(files), 6);

parfor numFile=1:length(files)
    if exist(fullfile(files(numFile).folder, 'unrolledGlands/gland_SR_basal/verticesInfo.mat'), 'file') == 0
        continue
    end
    [cells3dFeatures, gland3dFeatures, lumen3dFeatures, polygon_distribution_apical, polygon_distribution_basal, cellularFeatures, numCells, surfaceRatio2D, surfaceRatio3D, validCells, polygon_distribution_total] = obtainAllFeatures(files,numFile);
    %% Calculate mean and std of 3D features
    meanFeatures = varfun(@(x) mean(x),cells3dFeatures(:, 2:end));
    stdFeatures = varfun(@(x) std(x),cells3dFeatures(:, 2:end));
    
    cellularFeaturesNoCells = varfun(@(x) cell2mat(x), cellularFeatures(:, 4:5));
    cellularFeatures.Total_neighbours = cellularFeaturesNoCells.Fun_Total_neighbours;
    cellularFeatures.Apicobasal_neighbours = cellularFeaturesNoCells.Fun_Apicobasal_neighbours;
    mean3DNeighsFeatures = varfun(@(x) mean(x), cellularFeatures(validCells, 2:7));
    std3DNeighsFeatures = varfun(@(x) std(x), cellularFeatures(validCells, 2:7));
    
    totalMeanFeatures(numFile, :) = table2cell(meanFeatures);
    totalStdFeatures(numFile, :) = table2cell(stdFeatures);
    allGlands(numFile, :) = table2cell([gland3dFeatures, cell2table(polygon_distribution_apical(2, :), 'VariableNames', strcat('apical_', polygon_distribution_apical(1, :))), cell2table(polygon_distribution_basal(2, :), 'VariableNames', strcat('basal_', polygon_distribution_basal(1, :))), cell2table(polygon_distribution_total(2, :), 'VariableNames', strcat('total_', polygon_distribution_total(1, :)))]);
    allLumens(numFile, :) = table2cell(lumen3dFeatures);
    totalMean3DNeighsFeatures(numFile, :) = table2cell(mean3DNeighsFeatures);
    totalSTD3DNeighsFeatures(numFile, :) = table2cell(std3DNeighsFeatures);
    
    fileName = strsplit(files(numFile).folder, {'/','\'});
    fileName = convertCharsToStrings(strjoin({fileName{1,end-2},fileName{1,end-1}}, ' '));
    allGeneralInfo(numFile, :) = [{fileName}, {surfaceRatio3D}, {surfaceRatio2D}, {numCells}];
end

if contains(folderName, 'Salivary gland') == 0
    allGeneralInfo = cell2table(allGeneralInfo, 'VariableNames', {'ID_Glands', 'SurfaceRatio3D', 'SurfaceRatio2D', 'NCells'});
    allGlands = cell2table(allGlands, 'VariableNames', {'ID_Cell','Volume','EquivDiameter','PrincipalAxisLength','ConvexVolume','Solidity','SurfaceArea','aspectRatio','sphericity','normalizedVolume','basalNumNeighs','basal_area_cells2D','apicalNumNeighs','apical_area_cells2D','percentageScutoids','totalNeighs','apicoBasalTransitions','apical_triangles','apical_squares','apical_pentagons','apical_hexagons','apical_heptagons','apical_octogons','apical_nonagons','apical_decagons','basal_triangles','basal_squares','basal_pentagons','basal_hexagons','basal_heptagons','basal_octogons','basal_nonagons','basal_decagons', 'total_triangles','total_squares','total_pentagons','total_hexagons','total_heptagons','total_octogons','total_nonagons','total_decagons'});
    allLumens = cell2table(allLumens, 'VariableNames', {'ID_Cell','Volume','EquivDiameter','PrincipalAxisLength','ConvexVolume','Solidity','SurfaceArea','aspectRatio','sphericity','normalizedVolume','basalNumNeighs','basal_area_cells2D','apicalNumNeighs','apical_area_cells2D','percentageScutoids','totalNeighs','apicoBasalTransitions'});
    totalMeanFeatures = cell2table(totalMeanFeatures, 'VariableNames', {'Fun_Volume','Fun_EquivDiameter','Fun_PrincipalAxisLength','Fun_ConvexVolume','Fun_Solidity','Fun_SurfaceArea','Fun_aspectRatio','Fun_sphericity','Fun_normalizedVolume','Fun_apicalNumNeighs','Fun_apical_area_cells2D','Fun_basalNumNeighs','Fun_basal_area_cells2D','Fun_percentageScutoids','Fun_totalNeighs','Fun_apicoBasalTransitions'});
    totalStdFeatures = cell2table(totalStdFeatures, 'VariableNames', {'Fun_Volume','Fun_EquivDiameter','Fun_PrincipalAxisLength','Fun_ConvexVolume','Fun_Solidity','Fun_SurfaceArea','Fun_aspectRatio','Fun_sphericity','Fun_normalizedVolume','Fun_apicalNumNeighs','Fun_apical_area_cells2D','Fun_basalNumNeighs','Fun_basal_area_cells2D','Fun_percentageScutoids','Fun_totalNeighs','Fun_apicoBasalTransitions'});
    totalMean3DNeighsFeatures = cell2table(totalMean3DNeighsFeatures, 'VariableNames', {'Fun_Apical_sides','Fun_Basal_sides','Fun_Total_neighbours','Fun_Apicobasal_neighbours','Fun_Scutoids','Fun_apicoBasalTransitions'});
    totalSTD3DNeighsFeatures = cell2table(totalSTD3DNeighsFeatures, 'VariableNames', {'Fun_Apical_sides','Fun_Basal_sides','Fun_Total_neighbours','Fun_Apicobasal_neighbours','Fun_Scutoids','Fun_apicoBasalTransitions'});
    
    allGlands.Properties.VariableNames = cellfun(@(x) strcat('Gland_', x), allGlands.Properties.VariableNames, 'UniformOutput', false);
    allLumens.Properties.VariableNames = cellfun(@(x) strcat('Lumen_', x), allLumens.Properties.VariableNames, 'UniformOutput', false);
    totalMeanFeatures.Properties.VariableNames = cellfun(@(x) strcat('AverageCell_', x(5:end)), totalMeanFeatures.Properties.VariableNames, 'UniformOutput', false);
    totalStdFeatures.Properties.VariableNames = cellfun(@(x) strcat('STDCell_', x(5:end)), totalStdFeatures.Properties.VariableNames, 'UniformOutput', false);

    totalMean3DNeighsFeatures.Properties.VariableNames = cellfun(@(x) strcat('AverageCell_3D', x(5:end)), totalMean3DNeighsFeatures.Properties.VariableNames, 'UniformOutput', false);
    totalSTD3DNeighsFeatures.Properties.VariableNames = cellfun(@(x) strcat('STDCell_3D', x(5:end)), totalSTD3DNeighsFeatures.Properties.VariableNames, 'UniformOutput', false);
    
    save(fullfile(selpath(1).folder, 'global_3dFeatures.mat'), 'allGeneralInfo', 'totalMeanFeatures','totalStdFeatures', 'allLumens', 'allGlands')
    writetable([allGeneralInfo,totalMeanFeatures,totalStdFeatures, totalMean3DNeighsFeatures, totalSTD3DNeighsFeatures, allGlands, allLumens], fullfile(selpath(1).folder,'global_3dFeatures.xls'),'Range','B2');
end
end

