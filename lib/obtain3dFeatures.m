function [cells3dFeatures, gland3dFeatures, lumen3dFeatures, polygon_distribution_apical, polygon_distribution_basal, cellularFeatures, numCells, surfaceRatio2D, surfaceRatio3D, validCells, polygon_distribution_total] = obtainAllFeatures(files,numFile)
%OBTAINALLFEATURES Summary of this function goes here
%   Detailed explanation goes here
    load(fullfile(files(numFile).folder, 'valid_cells.mat'));
    if exist(fullfile(files(numFile).folder, 'morphological3dFeatures.mat'), 'file') == 0
        files(numFile).folder
        load(fullfile(files(numFile).folder, '3d_layers_info.mat'))%, 'labelledImage_realSize', 'lumenImage_realSize');
        
        if contains(files(numFile).folder, 'flatten')
            basalLayer = getBasalFrom3DImage(labelledImage_realSize, lumenImage_realSize, 0, labelledImage_realSize == 0 & lumenImage_realSize == 0);
            [apicalLayer] = getApicalFrom3DImage(lumenImage_realSize, labelledImage_realSize);
            [apical3dInfo] = calculateNeighbours3D(apicalLayer, 2, apicalLayer == 0);
            apical3dInfo = apical3dInfo.neighbourhood';

            [basal3dInfo] = calculateNeighbours3D(basalLayer, 2, basalLayer == 0);
            basal3dInfo = basal3dInfo.neighbourhood';
            
            cellularFeatures = calculate_CellularFeatures(apical3dInfo,basal3dInfo,apicalLayer,basalLayer,labelledImage_realSize,noValidCells,validCells,[]);
        else
            labelledImage_realSize = imresize3(labelledImage, 4.06, 'nearest');
            lumenImage_realSize = imresize3(double(lumenImage), 4.06, 'nearest');  
            cellularFeatures = calculate_CellularFeatures(apical3dInfo,basal3dInfo,apicalLayer,basalLayer,labelledImage,noValidCells,validCells,[]);
        end
        
        [apicalLayer] = getApicalFrom3DImage(lumenImage_realSize, labelledImage_realSize);    

        apicalLayer_onlyValidCells = ismember(apicalLayer, validCells) .* apicalLayer;
        apical_area_cells3D = cell2mat(struct2cell(regionprops(apicalLayer_onlyValidCells,'Area'))).';
        apical_area_cells3D = apical_area_cells3D(validCells);
        
        basalLayer_onlyValidCells = ismember(basalLayer, validCells) .* basalLayer;
        basal_area_cells3D = cell2mat(struct2cell(regionprops(basalLayer_onlyValidCells,'Area'))).';
        basal_area_cells3D = basal_area_cells3D(validCells);
        
        surfaceRatio3D = sum(ismember(basalLayer(:), validCells)) / sum(ismember(apicalLayer(:), validCells));
        
        %% Basal features
        basalNumNeighs = cellfun(@(x) length(x), basal3dInfo)'; 
        [polygon_distribution_basal] = calculate_polygon_distribution(basalNumNeighs, validCells);
        basalNumNeighs = -1 * ones(size(basal_area_cells3D));
        basal_area_cells2D = -1 * ones(size(basal_area_cells3D));
        
        basalInfo = table(basalNumNeighs, basal_area_cells2D, basal_area_cells3D);
        
        %% Apical features
        apicalNumNeighs = cellfun(@(x) length(x), apical3dInfo)';
        [polygon_distribution_apical] = calculate_polygon_distribution(apicalNumNeighs, validCells);
        apicalNumNeighs =  -1 * ones(size(apical_area_cells3D));
        apical_area_cells2D= -1 * ones(size(apical_area_cells3D));
        
        apicalInfo = table(apicalNumNeighs, apical_area_cells2D, apical_area_cells3D);
        
        %% Total features
        percentageScutoids = -1 * ones(size(apical_area_cells3D));
        totalNeighs = -1 * ones(size(apical_area_cells3D));
        apicoBasalTransitions = -1 * ones(size(apical_area_cells3D));
        polygon_distribution_total = calculate_polygon_distribution(cellfun(@(x,y) length(unique([x;y])), cells3dFeatures.Apical_sides, cells3dFeatures.Basal_sides), validCells);

        %% Extract each cell and calculate 3D features
        [cells3dFeatures] = extract3dDescriptors(labelledImage_realSize, validCells');
        cells3dFeatures = horzcat(cells3dFeatures, apicalInfo, basalInfo, table(percentageScutoids, totalNeighs, apicoBasalTransitions));
        
        %% Lumen features
        [lumen3dFeatures] = extract3dDescriptors(lumenImage_realSize, 1);
        lumen3dFeatures.ID_Cell = 'Lumen';
        lumen3dFeatures.basalNumNeighs = -1;
        lumen3dFeatures.basal_area_cells2D = -1;
        lumen3dFeatures.basal_area_cells3D = -1;
        lumen3dFeatures.apicalNumNeighs = -1;
        lumen3dFeatures.apical_area_cells2D = -1;
        lumen3dFeatures.apical_area_cells3D = -1;
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
        gland3dFeatures.basal_area_cells3D = -1;
        gland3dFeatures.apicalNumNeighs = -1;
        gland3dFeatures.apical_area_cells2D = -1;
        gland3dFeatures.apical_area_cells3D = -1;
        gland3dFeatures.percentageScutoids = -1;
        gland3dFeatures.totalNeighs = -1;
        gland3dFeatures.apicoBasalTransitions = -1;
        
        numCells = length(validCells);
        surfaceRatio2D = sum(basal_area_cells2D) / sum(apical_area_cells2D);
        
        allFeatures = vertcat(cells3dFeatures, gland3dFeatures, lumen3dFeatures);
        %% Save variables and export to excel
        writetable(allFeatures,fullfile(files(numFile).folder,'3dFeatures_LimeSeg3DSegmentation.xls'), 'Range','B2');
        save(fullfile(files(numFile).folder, 'morphological3dFeatures.mat'), 'cells3dFeatures', 'gland3dFeatures', 'lumen3dFeatures', 'polygon_distribution_apical', 'polygon_distribution_basal', 'cellularFeatures', 'numCells', 'surfaceRatio2D', 'surfaceRatio3D', 'polygon_distribution_total');
    else
        load(fullfile(files(numFile).folder, 'morphological3dFeatures.mat'));
    end
end
