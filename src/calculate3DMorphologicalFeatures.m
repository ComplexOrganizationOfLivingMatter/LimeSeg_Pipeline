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

totalMeanFeatures = [];
totalStdFeatures = [];
allGlands = [];
allLumens = [];
allGeneralInfo = {};
totalSTD3DNeighsFeatures = [];
totalMean3DNeighsFeatures = [];

for numFiles=1:length(files)
    load(fullfile(files(numFiles).folder, 'valid_cells.mat'));
    if exist(fullfile(files(numFiles).folder, 'morphological3dFeatures.mat'), 'file') ~= 0
        files(numFiles).folder
        if exist(fullfile(files(numFiles).folder, 'unrolledGlands/gland_SR_basal/verticesInfo.mat'), 'file') == 0
            continue
        end
        load(fullfile(files(numFiles).folder, '3d_layers_info.mat'))%, 'labelledImage_realSize', 'lumenImage_realSize');
        
        cellularFeatures = calculate_CellularFeatures(apical3dInfo,basal3dInfo,apicalLayer,basalLayer,labelledImage,noValidCells,validCells,[]);
        
        surfaceRatio3D = sum(basalLayer(:)) / sum(apicalLayer(:));
        
        %% Basal features
        load(fullfile(files(numFiles).folder, 'unrolledGlands/gland_SR_basal/verticesInfo.mat'), 'newVerticesNeighs2D', 'cylindre2DImage');
        basalNeighs = getNeighboursFromVertices(newVerticesNeighs2D);
        basalNumNeighs = cellfun(@(x) length(x), basalNeighs)';
        [polygon_distribution_basal] = calculate_polygon_distribution(basalNumNeighs, validCells);
        basalNumNeighs = basalNumNeighs(validCells);
        basal_area_cells2D=cell2mat(struct2cell(regionprops(cylindre2DImage,'Area'))).';
        basal_area_cells2D = basal_area_cells2D(validCells);
        
        basalInfo = table(basalNumNeighs, basal_area_cells2D);
        
        %% Apical features
        load(fullfile(files(numFiles).folder, 'unrolledGlands/gland_SR_1/verticesInfo.mat'), 'newVerticesNeighs2D', 'cylindre2DImage');
        apicalNeighs = getNeighboursFromVertices(newVerticesNeighs2D);
        apicalNumNeighs = cellfun(@(x) length(x), apicalNeighs)';
        [polygon_distribution_apical] = calculate_polygon_distribution(apicalNumNeighs, validCells);
        apicalNumNeighs = apicalNumNeighs(validCells);
        apical_area_cells2D=cell2mat(struct2cell(regionprops(cylindre2DImage,'Area'))).';
        apical_area_cells2D = apical_area_cells2D(validCells);
        
        apicalInfo = table(apicalNumNeighs, apical_area_cells2D);
        
        %% Total features
        percentageScutoids = cellfun(@(x, y) ~isempty(setxor(x,y)), apicalNeighs(validCells), basalNeighs(validCells))';
        totalNeighs = cellfun(@(x,y) length(unique([x;y])), apicalNeighs(validCells), basalNeighs(validCells))';
        apicoBasalTransitions = cellfun(@(x, y) length(unique(vertcat(setdiff(x, y), setdiff(y, x)))), apicalNeighs(validCells), basalNeighs(validCells))';
        
        %% Extract each cell and calculate 3D features
        [cells3dFeatures] = extract3dDescriptors(labelledImage_realSize, validCells');
        cells3dFeatures = horzcat(cells3dFeatures, apicalInfo, basalInfo, table(percentageScutoids, totalNeighs, apicoBasalTransitions));

        %% Lumen features
        [lumen3dFeatures] = extract3dDescriptors(lumenImage_realSize, 1);
        lumen3dFeatures.ID_Cell = 'Lumen';
        lumen3dFeatures.basalNumNeighs = -1;
        lumen3dFeatures.basal_area_cells = -1;
        lumen3dFeatures.apicalNumNeighs = -1;
        lumen3dFeatures.apical_area_cells = -1;
        lumen3dFeatures.percentageScutoids = -1;
        lumen3dFeatures.totalNeighs = -1;
        lumen3dFeatures.apicoBasalTransitions = -1;

        %% Global Gland
        % We need calculate thickness of the glands or number of cell in
        % transversal axis
        [gland3dFeatures] = extract3dDescriptors(labelledImage_realSize>0, 1);
        gland3dFeatures.ID_Cell = 'Gland';
        gland3dFeatures.basalNumNeighs = -1;
        gland3dFeatures.basal_area_cells = -1;
        gland3dFeatures.apicalNumNeighs = -1;
        gland3dFeatures.apical_area_cells = -1;
        gland3dFeatures.percentageScutoids = -1;
        gland3dFeatures.totalNeighs = -1;
        gland3dFeatures.apicoBasalTransitions = -1;
        
        numCells = length(validCells);
        surfaceRatio2D = sum(basal_area_cells2D) / sum(apical_area_cells2D);
        
        allFeatures = vertcat(cells3dFeatures, gland3dFeatures, lumen3dFeatures);
        %% Save variables and export to excel
        writetable(allFeatures,fullfile(files(numFiles).folder,'3dFeatures_LimeSeg3DSegmentation.xls'), 'Range','B2');
        save(fullfile(files(numFiles).folder, 'morphological3dFeatures.mat'), 'cells3dFeatures', 'gland3dFeatures', 'lumen3dFeatures', 'polygon_distribution_apical', 'polygon_distribution_basal', 'cellularFeatures', 'numCells', 'surfaceRatio2D', 'surfaceRatio3D');
    else
        load(fullfile(files(numFiles).folder, 'morphological3dFeatures.mat'));
    end
    
    %% Calculate mean and std of 3D features
    meanFeatures = varfun(@(x) mean(x),cells3dFeatures(:, 2:end));
    stdFeatures = varfun(@(x) std(x),cells3dFeatures(:, 2:end));
    
    cellularFeaturesNoCells = varfun(@(x) cell2mat(x), cellularFeatures(:, 4:5));
    cellularFeatures.Total_neighbours = cellularFeaturesNoCells.Fun_Total_neighbours;
    cellularFeatures.Apicobasal_neighbours = cellularFeaturesNoCells.Fun_Apicobasal_neighbours;
    mean3DNeighsFeatures = varfun(@(x) mean(x), cellularFeatures(validCells, 2:6));
    std3DNeighsFeatures = varfun(@(x) std(x), cellularFeatures(validCells, 2:6));
    
    totalMeanFeatures = vertcat(totalMeanFeatures, meanFeatures);
    totalStdFeatures = vertcat(totalStdFeatures, stdFeatures);
    allGlands = vertcat(allGlands, [gland3dFeatures, cell2table(polygon_distribution_apical(2, :), 'VariableNames', strcat('apical_', polygon_distribution_apical(1, :))), cell2table(polygon_distribution_basal(2, :), 'VariableNames', strcat('basal_', polygon_distribution_basal(1, :)))]);
    allLumens = vertcat(allLumens, lumen3dFeatures);
    totalMean3DNeighsFeatures = vertcat(totalMean3DNeighsFeatures, mean3DNeighsFeatures);
    totalSTD3DNeighsFeatures = vertcat(totalSTD3DNeighsFeatures, std3DNeighsFeatures);
    
    fileName = strsplit(files(numFiles).folder, {'/','\'});
    fileName = convertCharsToStrings(strjoin({fileName{1,end-2},fileName{1,end-1}}, ' '));
    allGeneralInfo = [allGeneralInfo ; {fileName}, {surfaceRatio3D}, {surfaceRatio2D}, {numCells}];
end

if contains(folderName, 'Salivary gland') == 0
    allGeneralInfo = cell2table(allGeneralInfo, 'VariableNames', {'ID_Glands', 'SurfaceRatio3D', 'SurfaceRatio2D', 'NCells'});

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

