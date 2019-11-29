function [cells3dFeatures, gland3dFeatures, lumen3dFeatures, polygon_distribution_apical, polygon_distribution_basal, cellularFeatures, numCells, surfaceRatio2D, surfaceRatio3D, validCells] = obtainAllFeatures(files,numFile)
%OBTAINALLFEATURES Summary of this function goes here
%   Detailed explanation goes here
    load(fullfile(files(numFile).folder, 'valid_cells.mat'));
    if exist(fullfile(files(numFile).folder, 'morphological3dFeatures.mat'), 'file') ~= 0
        files(numFile).folder
        load(fullfile(files(numFile).folder, '3d_layers_info.mat'))%, 'labelledImage_realSize', 'lumenImage_realSize');
        
        cellularFeatures = calculate_CellularFeatures(apical3dInfo,basal3dInfo,apicalLayer,basalLayer,labelledImage,noValidCells,validCells,[]);
        
        load(fullfile(files(numFile).folder, 'unrolledGlands/gland_SR_basal/final3DImg.mat'), 'img3d');
        basalLayer = img3d;        
        load(fullfile(files(numFile).folder, 'unrolledGlands/gland_SR_1/final3DImg.mat'), 'img3d');
        surfaceRatio3D = sum(ismember(basalLayer(:), validCells)) / sum(ismember(img3d(:), validCells));
        
        %% Basal features
        load(fullfile(files(numFile).folder, 'unrolledGlands/gland_SR_basal/verticesInfo.mat'), 'newVerticesNeighs2D', 'cylindre2DImage');
        basalNeighs = getNeighboursFromVertices(newVerticesNeighs2D);
        basalNumNeighs = cellfun(@(x) length(x), basalNeighs)';
        [polygon_distribution_basal] = calculate_polygon_distribution(basalNumNeighs, validCells);
        basalNumNeighs = basalNumNeighs(validCells);
        basal_area_cells2D=cell2mat(struct2cell(regionprops(cylindre2DImage,'Area'))).';
        basal_area_cells2D = basal_area_cells2D(validCells);
        
        basalInfo = table(basalNumNeighs, basal_area_cells2D);
        
        %% Apical features
        load(fullfile(files(numFile).folder, 'unrolledGlands/gland_SR_1/verticesInfo.mat'), 'newVerticesNeighs2D', 'cylindre2DImage');
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
        lumen3dFeatures.basal_area_cells2D = -1;
        lumen3dFeatures.apicalNumNeighs = -1;
        lumen3dFeatures.apical_area_cells2D = -1;
        lumen3dFeatures.percentageScutoids = -1;
        lumen3dFeatures.totalNeighs = -1;
        lumen3dFeatures.apicoBasalTransitions = -1;

        %% Global Gland
        % We need calculate thickness of the glands or number of cell in
        % transversal axis
        [gland3dFeatures] = extract3dDescriptors(labelledImage_realSize>0, 1);
        gland3dFeatures.ID_Cell = 'Gland';
        gland3dFeatures.basalNumNeighs = -1;
        gland3dFeatures.basal_area_cells2D = -1;
        gland3dFeatures.apicalNumNeighs = -1;
        gland3dFeatures.apical_area_cells2D = -1;
        gland3dFeatures.percentageScutoids = -1;
        gland3dFeatures.totalNeighs = -1;
        gland3dFeatures.apicoBasalTransitions = -1;
        
        numCells = length(validCells);
        surfaceRatio2D = sum(basal_area_cells2D) / sum(apical_area_cells2D);
        
        allFeatures = vertcat(cells3dFeatures, gland3dFeatures, lumen3dFeatures);
        %% Save variables and export to excel
        writetable(allFeatures,fullfile(files(numFile).folder,'3dFeatures_LimeSeg3DSegmentation.xls'), 'Range','B2');
        save(fullfile(files(numFile).folder, 'morphological3dFeatures.mat'), 'cells3dFeatures', 'gland3dFeatures', 'lumen3dFeatures', 'polygon_distribution_apical', 'polygon_distribution_basal', 'cellularFeatures', 'numCells', 'surfaceRatio2D', 'surfaceRatio3D');
    else
        load(fullfile(files(numFile).folder, 'morphological3dFeatures.mat'));
    end
end

