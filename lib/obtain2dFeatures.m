function [allFeatures,polygon_distribution_basal,polygon_distribution_apical,polygon_distribution_total,surfaceRatio2D] = obtain2dFeatures(files,numFile, allFeatures, validCells)  
%OBTAIN2DFEATURES Summary of this function goes here
%   Detailed explanation goes here

        %% Basal features
        load(fullfile(files(numFile).folder, 'unrolledGlands/gland_SR_basal/verticesInfo.mat'), 'newVerticesNeighs2D', 'cylindre2DImage');
        basal_area_cells2D=cell2mat(struct2cell(regionprops(cylindre2DImage,'Area'))).';
        allFeatures.basal_area_cells2D(1:end-2,:) = basal_area_cells2D(validCells);
        basalNeighs = getNeighboursFromVertices(newVerticesNeighs2D);
        
        %% Apical features
        load(fullfile(files(numFile).folder, 'unrolledGlands/gland_SR_1/verticesInfo.mat'), 'newVerticesNeighs2D', 'cylindre2DImage');
        apicalNeighs = getNeighboursFromVertices(newVerticesNeighs2D);
        
        if size(apicalNeighs, 2) > size(basalNeighs, 2)
            basalNeighs(size(apicalNeighs, 2)) = {[]};
        elseif  size(apicalNeighs, 2) < size(basalNeighs, 2)
            apicalNeighs(size(basalNeighs, 2)) = {[]};
        end
        
        basalNumNeighs = cellfun(@(x) length(x), basalNeighs)';
        [polygon_distribution_basal] = calculate_polygon_distribution(basalNumNeighs, validCells);
        allFeatures.basalNumNeighs(1:end-2,:) = basalNumNeighs(validCells);
        
        apicalNumNeighs = cellfun(@(x) length(x), apicalNeighs)';
        [polygon_distribution_apical] = calculate_polygon_distribution(apicalNumNeighs, validCells);
        allFeatures.apicalNumNeighs(1:end-2,:) = apicalNumNeighs(validCells);
        apical_area_cells2D=cell2mat(struct2cell(regionprops(cylindre2DImage,'Area'))).';
        allFeatures.apical_area_cells2D(1:end-2,:) = apical_area_cells2D(validCells);
        
        %% Total features
        allFeatures.percentageScutoids(1:end-2,:) = cellfun(@(x, y) ~isempty(setxor(x,y)), apicalNeighs(validCells), basalNeighs(validCells))';
        allFeatures.totalNeighs(1:end-2,:) = cellfun(@(x,y) length(unique([x;y])), apicalNeighs(validCells), basalNeighs(validCells))';
        allFeatures.apicoBasalTransitions(1:end-2,:) = cellfun(@(x, y) length(unique(vertcat(setdiff(x, y), setdiff(y, x)))), apicalNeighs(validCells), basalNeighs(validCells))';
        polygon_distribution_total = calculate_polygon_distribution(cellfun(@(x,y) length(unique([x;y])), apicalNeighs, basalNeighs), validCells);
        
        %% Global Gland
        surfaceRatio2D = sum(allFeatures.basal_area_cells2D) / sum(allFeatures.apical_area_cells2D);
        
end

