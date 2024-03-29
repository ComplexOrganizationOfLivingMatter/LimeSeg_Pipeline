function [] = calculate3DMorphologicalFeatures(folderName, TissueName)
%CALCULATE3DMORPHOLOGICALFEATURES Summary of this function goes here
%   Detailed explanation goes here

if contains(folderName, strcat('\data\',TissueName,'/'))
    files = dir(fullfile(folderName, '3d_layers_info.mat'));
else
    files = dir(fullfile('**/data/',TissueName,'/', folderName, '/**/Results/3d_layers_info.mat'));
end
nonDiscardedFiles = cellfun(@(x) contains(lower(x), 'discarded') == 0, {files.folder});
files = files(nonDiscardedFiles);

selpath = dir(fullfile('**/data/',TissueName,'/', folderName));

totalMeanFeatures = cell([length(files) 18]);
totalStdFeatures = cell([length(files) 18]);
allGlands = cell([length(files) 43]);
allLumens = cell([length(files) 19]);
allGeneralInfo = cell(length(files), 4);
allNetworkFeatures = cell(length(files),5);
totalSTD3DNeighsFeatures = cell(length(files), 6);
totalMean3DNeighsFeatures = cell(length(files), 6);
allHollowGland3dFeatures = cell([length(files) 19]);

for numFile=1:length(files)
   [cells3dFeatures, gland3dFeatures, lumen3dFeatures,hollowGland3dFeatures, polygon_distribution_apical, polygon_distribution_basal, cellularFeatures, numCells, surfaceRatio2D, surfaceRatio3D, validCells, polygon_distribution_total,apicoBasalNeighs] = obtain3dFeatures(files,numFile);
    %% Calculate mean and std of 3D features
    cells3dFeatures((cells3dFeatures.ID_Cell == "Lumen" | cells3dFeatures.ID_Cell == "Gland and Lumen"),:)=[];
    meanFeatures = varfun(@(x) mean(x),cells3dFeatures(:, 2:end));
    stdFeatures = varfun(@(x) std(x),cells3dFeatures(:, 2:end));
    
    [meanFeatures,stdFeatures, gland3dFeatures, lumen3dFeatures,hollowGland3dFeatures] = convertPixelsToMicrons(files,numFile,meanFeatures,stdFeatures, gland3dFeatures, lumen3dFeatures,hollowGland3dFeatures);
    
    validCells=setdiff(validCells,find(cellularFeatures.Apical_sides==0 & cellularFeatures.Basal_sides==0));
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
    [degreeNodesCorrelation,coefCluster,betweennessCentrality] = obtainNetworksFeatures(apicoBasalNeighs,validCells, fullfile(files(numFile).folder, 'network3dFeatures.mat'));
    allNetworkFeatures(numFile,:) = [{mean(coefCluster)}, {mean(betweennessCentrality)},{degreeNodesCorrelation} {std(coefCluster)},{std(betweennessCentrality)}];
    allHollowGland3dFeatures(numFile, :) = table2cell(hollowGland3dFeatures);
end

if contains(folderName, 'Salivary gland') == 0
    allGeneralInfo = cell2table(allGeneralInfo, 'VariableNames', {'ID_Glands', 'SurfaceRatio3D', 'SurfaceRatio2D', 'NCells'});
    allGlands = cell2table(allGlands, 'VariableNames', {'ID_Cell','Volume','EquivDiameter','PrincipalAxisLength','ConvexVolume','Solidity','SurfaceArea','aspectRatio','sphericity','normalizedVolume','basalNumNeighs','basal_area_cells2D','basal_area_cells3D', 'apicalNumNeighs','apical_area_cells2D','apical_area_cells3D','percentageScutoids','totalNeighs','apicoBasalTransitions','apical_triangles','apical_squares','apical_pentagons','apical_hexagons','apical_heptagons','apical_octogons','apical_nonagons','apical_decagons','basal_triangles','basal_squares','basal_pentagons','basal_hexagons','basal_heptagons','basal_octogons','basal_nonagons','basal_decagons', 'total_triangles','total_squares','total_pentagons','total_hexagons','total_heptagons','total_octogons','total_nonagons','total_decagons'});
    allLumens = cell2table(allLumens, 'VariableNames', {'ID_Cell','Volume','EquivDiameter','PrincipalAxisLength','ConvexVolume','Solidity','SurfaceArea','aspectRatio','sphericity','normalizedVolume','basalNumNeighs','basal_area_cells2D','basal_area_cells3D', 'apicalNumNeighs','apical_area_cells2D','apical_area_cells3D','percentageScutoids','totalNeighs','apicoBasalTransitions'});
    allHollowGland3dFeatures = cell2table(allHollowGland3dFeatures, 'VariableNames',  {'ID_Cell','Volume','EquivDiameter','PrincipalAxisLength','ConvexVolume','Solidity','SurfaceArea','aspectRatio','sphericity','normalizedVolume','basalNumNeighs','basal_area_cells2D','basal_area_cells3D', 'apicalNumNeighs','apical_area_cells2D','apical_area_cells3D','percentageScutoids','totalNeighs','apicoBasalTransitions'});
    allNetworkFeatures= cell2table(allNetworkFeatures, 'VariableNames', {'Coeficient_Cluster', 'Betweenness_Centrality', 'Assortativity', 'coeficient_Cluster_STD','Betweenness_Centrality_STD'});
    totalMeanFeatures = cell2table(totalMeanFeatures, 'VariableNames', {'Fun_Volume','Fun_EquivDiameter','Fun_PrincipalAxisLength','Fun_ConvexVolume','Fun_Solidity','Fun_SurfaceArea','Fun_aspectRatio','Fun_sphericity','Fun_normalizedVolume','Fun_apicalNumNeighs','Fun_apical_area_cells2D','Fun_apical_area_cells3D','Fun_basalNumNeighs','Fun_basal_area_cells2D','Fun_basal_area_cells3D','Fun_percentageScutoids','Fun_totalNeighs','Fun_apicoBasalTransitions'});
    totalStdFeatures = cell2table(totalStdFeatures, 'VariableNames', {'Fun_Volume','Fun_EquivDiameter','Fun_PrincipalAxisLength','Fun_ConvexVolume','Fun_Solidity','Fun_SurfaceArea','Fun_aspectRatio','Fun_sphericity','Fun_normalizedVolume','Fun_apicalNumNeighs','Fun_apical_area_cells2D','Fun_apical_area_cells3D','Fun_basalNumNeighs','Fun_basal_area_cells2D','Fun_basal_area_cells3D','Fun_percentageScutoids','Fun_totalNeighs','Fun_apicoBasalTransitions'});
    totalMean3DNeighsFeatures = cell2table(totalMean3DNeighsFeatures, 'VariableNames', {'Fun_Apical_sides','Fun_Basal_sides','Fun_Total_neighbours','Fun_Apicobasal_neighbours','Fun_Scutoids','Fun_apicoBasalTransitions'});
    totalSTD3DNeighsFeatures = cell2table(totalSTD3DNeighsFeatures, 'VariableNames', {'Fun_Apical_sides','Fun_Basal_sides','Fun_Total_neighbours','Fun_Apicobasal_neighbours','Fun_Scutoids','Fun_apicoBasalTransitions'});
    
    allGlands.Properties.VariableNames = cellfun(@(x) strcat('Gland_', x), allGlands.Properties.VariableNames, 'UniformOutput', false);
    allHollowGland3dFeatures.Properties.VariableNames = cellfun(@(x) strcat('HollowGland_', x), allHollowGland3dFeatures.Properties.VariableNames, 'UniformOutput', false);
    allLumens.Properties.VariableNames = cellfun(@(x) strcat('Lumen_', x), allLumens.Properties.VariableNames, 'UniformOutput', false);
    totalMeanFeatures.Properties.VariableNames = cellfun(@(x) strcat('AverageCell_', x(5:end)), totalMeanFeatures.Properties.VariableNames, 'UniformOutput', false);
    totalStdFeatures.Properties.VariableNames = cellfun(@(x) strcat('STDCell_', x(5:end)), totalStdFeatures.Properties.VariableNames, 'UniformOutput', false);
    
    totalMean3DNeighsFeatures.Properties.VariableNames = cellfun(@(x) strcat('AverageCell_3D', x(5:end)), totalMean3DNeighsFeatures.Properties.VariableNames, 'UniformOutput', false);
    totalSTD3DNeighsFeatures.Properties.VariableNames = cellfun(@(x) strcat('STDCell_3D', x(5:end)), totalSTD3DNeighsFeatures.Properties.VariableNames, 'UniformOutput', false);

    save(fullfile(selpath(1).folder, 'global_3dFeatures.mat'), 'allGeneralInfo', 'totalMeanFeatures','totalStdFeatures', 'allLumens', 'allGlands', 'totalMean3DNeighsFeatures', 'totalSTD3DNeighsFeatures', 'allNetworkFeatures', 'allHollowGland3dFeatures')
    allFeatures = [allGeneralInfo,totalMeanFeatures,totalStdFeatures, allGlands, allLumens, totalMean3DNeighsFeatures, totalSTD3DNeighsFeatures];
    % The order of the head is the following: nCell, (Basal, apical area and SR 2D), (Basal, apical area and SR 3D), (Basal, apical and
    % apicobasal N-2D), (basal apical and apicobasal N-3D),Scutoids2D and 3D,apicobasalTransition 2D and 3D, 2D poligon
    % distribution, (Volumen,ConvexVolume and Solidity Cells),(Surface Area cell, 
    % Sphericity cells), AxisLength cells, AspectRatio cells, (Volumen,ConvexVolume and Solidity Gland),SurfaceArea and sphericity Gland, AxisLength Gland, AspectRatio Gland
    % (Volumen,ConvexVolume and Solidity Lumen), SurfaceArea and sphericity Lumen, 
    % AxisLength Lumen, AspectRatio Lumen,percentageLumenSpace,hollowGland Features, networkFeatures(coeficient clustering,betweenness centrality and assortativity)
    PercentageLumenSpace = table(table2array(allFeatures(:,85)) ./ table2array(allFeatures(:,42))); 
    PercentageLumenSpace.Properties.VariableNames = {'PercentageLumenSpace'};
    FeaturesPerCell=array2table([table2array(allFeatures(:,42))./ table2array(allFeatures(:,4)), table2array(allHollowGland3dFeatures(:,2))./ table2array(allFeatures(:,4)), table2array(allFeatures(:,85))./ table2array(allFeatures(:,4))]); 
    FeaturesPerCell.Properties.VariableNames = {'GlandVolume_perCell','HollowGlandVolume_perCell','LumenVolume_perCell'};
    
    finalTable = [allFeatures(:,1), allFeatures(:,4),allFeatures(:,15), allFeatures(:,18), allFeatures(:,3),allFeatures(:,16), allFeatures(:,19),allFeatures(:,2),allFeatures(:,14),allFeatures(:,17),allFeatures(:,21),totalMean3DNeighsFeatures(:,1),totalMean3DNeighsFeatures(:,2), totalMean3DNeighsFeatures(:,4), allFeatures(:,20),totalMean3DNeighsFeatures(:,5),allFeatures(:,22),totalMean3DNeighsFeatures(:,6), allFeatures(:,60:83),allFeatures(:,5),allFeatures(:,8:10),allFeatures(:,12),allFeatures(:,7),allFeatures(:,11),allFeatures(:,42),allFeatures(:,45:47),allFeatures(:,49),allFeatures(:,44), allFeatures(:,48),allFeatures(:,85),allFeatures(:,88:90),allFeatures(:,92),allFeatures(:,87),allFeatures(:,91),PercentageLumenSpace,allHollowGland3dFeatures(:,2),allHollowGland3dFeatures(:,5:7),allHollowGland3dFeatures(:,9),allHollowGland3dFeatures(:,4),allHollowGland3dFeatures(:,8),FeaturesPerCell,allNetworkFeatures(:,1:3)];
    finalSTDTable = ([totalStdFeatures(:,12),totalStdFeatures(:,15),totalSTD3DNeighsFeatures(:,1:2),totalSTD3DNeighsFeatures(:,4:5),allNetworkFeatures(:,4:5),totalStdFeatures(:,1),totalStdFeatures(:,4:5),totalStdFeatures(:,6),totalStdFeatures(:,8),totalStdFeatures(:,3),totalStdFeatures(:,7)]);
    
    writetable(finalTable, fullfile(selpath(1).folder,'global_3dFeatures.xls'),'Range','B2');
    writetable(finalSTDTable, fullfile(selpath(1).folder,'global_3dFeatures.xls'),'Sheet', 2,'Range','B2');
end
end

