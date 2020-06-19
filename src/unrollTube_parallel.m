function unrollTube_parallel(selpath)
%UNROLLTUBE_PARALLEL Perform unrolling all the layers
%   
    minNumberOfSurfaceRatios = 7;
    steps3D = 0.5;
    
    idName_splitted = strsplit(selpath, filesep);
    idName = strjoin(idName_splitted(end-3:end-1), '_');
    
    load(fullfile(selpath, 'valid_cells.mat'), 'validCells', 'noValidCells');
    
    if exist(strcat(selpath, '\', idName ,'_samirasFormat.xls'), 'file') == 0
        
        %% Unroll all the intermediate layers, although apical and basal may have already been unrolled
        if exist(fullfile(selpath, 'dividedGland' ,'glandDividedInSurfaceRatios.mat'), 'file') > 0
            %% Apical and basal and all the artificial surface ratios
            load(fullfile(selpath, 'dividedGland', 'glandDividedInSurfaceRatios.mat'));
            %%mkdir(fullfile(selpath, 'unrolledGlands', 'gland_SR_1'));
            infoPerSurfaceRatio(:, 1) = {[]};
            %% Get gland rotation using the gland lumen
            load(fullfile(selpath, 'unrolledGlands', 'gland_SR_1','apicalRotationsOriginal.mat'),'rotationsOriginal');
                      
            [samiraTablePerSR{1}, apicalAreaValidCells, ~] = unrollTube(infoPerSurfaceRatio{1, 3},[], fullfile(selpath, 'unrolledGlands', 'gland_SR_1'), fullfile(selpath, 'valid_cells.mat'), fullfile(selpath, '3d_layers_info.mat'),rotationsOriginal);
            areaValidCells{1} = apicalAreaValidCells;
            %addAttachedFiles(gcp, fullfile(selpath, 'valid_cells.mat'))
            
            surfaceRatioOfGland = vertcat(infoPerSurfaceRatio{:,2})';
            
            infoPerSurfaceRatio(1, :) = {[]};
            nSR = length(surfaceRatioOfGland);
            
            for numPartition = 2:nSR
                infoPerSurfaceRatio(numPartition, 3) = {ndSparse(infoPerSurfaceRatio{numPartition, 3})};
            end
            
            load(fullfile(selpath, 'realSize3dLayers.mat'), 'labelledImage_realSizeFlatten')
            for numPartition = 2:nSR
                if numPartition ~= nSR
                    [samiraTablePerSR{numPartition}, areaValidCells{numPartition}] = unrollTube(full(infoPerSurfaceRatio{numPartition, 3}),labelledImage_realSizeFlatten, fullfile(selpath, 'unrolledGlands', ['gland_SR_' num2str(numPartition)]), fullfile(selpath, 'valid_cells.mat'), fullfile(selpath, '3d_layers_info'), rotationsOriginal);
                else
                    [samiraTablePerSR{numPartition}, areaValidCells{numPartition}] = unrollTube(full(infoPerSurfaceRatio{numPartition, 3}),labelledImage_realSizeFlatten, fullfile(selpath, 'unrolledGlands', 'gland_SR_basal'), fullfile(selpath, 'valid_cells.mat'), fullfile(selpath, '3d_layers_info'), rotationsOriginal);
                end
            end
            
            %% Calculate theoretical surface ratio
            steps = 2.5/(minNumberOfSurfaceRatios-1);
            surfaceRatioOfGland = 1:steps:((steps*(minNumberOfSurfaceRatios-1))+1);
            
            if minNumberOfSurfaceRatios > length(samiraTablePerSR)
                minNumberOfSurfaceRatios = length(samiraTablePerSR);
                surfaceRatioOfGland = surfaceRatioOfGland(1:minNumberOfSurfaceRatios);
                realFinalStep = steps * ((infoPerSurfaceRatio{end, 2} - infoPerSurfaceRatio{end-1, 2}) / steps3D);
                surfaceRatioOfGland(length(samiraTablePerSR)) = surfaceRatioOfGland(end-1) + realFinalStep;
            end
            
            for numPartition = 2:minNumberOfSurfaceRatios
                sT_Actual = samiraTablePerSR{numPartition};
                sT_Actual(:, 1) = {surfaceRatioOfGland(numPartition)};
                samiraTablePerSR{numPartition} = sT_Actual;
            end
            
            samiraTable = vertcat(samiraTablePerSR{1:minNumberOfSurfaceRatios});
        end
        
        %% Creating samira table
        samiraTableT = cell2table(samiraTable, 'VariableNames',{'Radius', 'CellIDs', 'TipCells', 'BorderCell','verticesValues_x_y'});
        
        newCrossesTable = lookFor4cellsJunctionsAndExportTheExcel(samiraTableT);
        
        writetable(samiraTableT, strcat(selpath, '\', idName ,'_samirasFormat.xls'));
        writetable(newCrossesTable, strcat(selpath, '\', idName ,'_VertCrosses.xls'));
    end
    if exist('infoPerSurfaceRatio', 'var') 
        clearvars infoPerSurfaceRatio
    end
    close all
    %% Saving final information with the right order
    filesOf2DUnroll = dir(fullfile(selpath, '**', 'verticesInfo.mat'));
    if exist(fullfile(selpath, 'glandDividedInSurfaceRatios_AllUnrollFeatures.mat'), 'file') == 0
        load(fullfile(selpath, 'dividedGland' ,'glandDividedInSurfaceRatios.mat'), 'infoPerSurfaceRatio');
        
        load(fullfile(filesOf2DUnroll(1).folder, 'verticesInfo.mat'));
        apicalLayer = cylindre2DImage;
        infoPerSurfaceRatio{1, 6} = cylindre2DImage;
        [neighboursApical] = getNeighboursFromVertices(newVerticesNeighs2D);
        
        load(fullfile(filesOf2DUnroll(3).folder, 'verticesInfo.mat'));
        basalLayer = cylindre2DImage;
        nSR = size(infoPerSurfaceRatio, 1);
        infoPerSurfaceRatio{nSR, 6} = cylindre2DImage;
        [neighboursBasal] = getNeighboursFromVertices(newVerticesNeighs2D);
        
        neighboursOfAllSurfaces = cell(nSR, 1);
        for numSR = 1:nSR
            load(fullfile(filesOf2DUnroll(numSR).folder, 'final3DImg.mat'), 'img3d');
            load(fullfile(filesOf2DUnroll(numSR).folder, 'verticesInfo.mat'), 'cylindre2DImage', 'newVerticesNeighs2D')
            [total_neighbours3D] = getNeighboursFromVertices(newVerticesNeighs2D);
            midLayer = cylindre2DImage;
            
            splittedFolder = strsplit(filesOf2DUnroll(numSR).folder, '_');
            if isequal(splittedFolder{end}, 'basal') == 0
                idToSave = str2double(splittedFolder{end});
            else
                idToSave = nSR;
            end
            
            disp(['id to save: ' num2str(idToSave)])
            
            infoPerSurfaceRatio{idToSave, 6} = cylindre2DImage;
            [neighboursMid] = getNeighboursFromVertices(newVerticesNeighs2D);
            
            neighbours_data = table(neighboursApical, neighboursMid);
            neighbours_data.Properties.VariableNames = {'Apical','Basal'};
            neighboursOfAllSurfaces{idToSave} = neighboursMid;
            
            [infoPerSurfaceRatio{idToSave, 8}, infoPerSurfaceRatio{idToSave, 7}] = calculate_CellularFeatures(neighboursApical, neighboursMid, apicalLayer, midLayer, img3d, noValidCells, validCells, [], total_neighbours3D);
            
        end
        
        %% Calculate theoretical surface ratio
        surfaceRatioOfGland_real = vertcat(infoPerSurfaceRatio{:, 7})';
%         totalPartitions = nSR;
%         initialPartitions = (1:(totalPartitions-1))/totalPartitions;
%         surfaceRatioOfGland = surfaceRatioOfGland_real;
%         surfaceRatioOfGland(2:10) = initialPartitions * (surfaceRatioOfGland_real(end) - 1) + 1;
        infoPerSurfaceRatio(:, 7) = num2cell(surfaceRatioOfGland_real)';
        
        infoPerSurfaceRatio = cell2table(infoPerSurfaceRatio, 'VariableNames', {'Image3DWithVolumen', 'SR3D', 'Layer3D', 'ApicalBasalCellFeatures3D', 'BasalApicalCellFeatures3D', 'UnrolledLayer2D', 'SR2D', 'ApicalBasalCellFeatures2D'});
        %             save(fullfile(selpath, 'glandDividedInSurfaceRatios_AllUnrollFeatures.mat'), 'cellularFeatures_BasalToApical', 'cellularFeatures_ApicalToBasal', 'meanSurfaceRatio');
        save(fullfile(selpath, 'glandDividedInSurfaceRatios_AllUnrollFeatures.mat'), 'infoPerSurfaceRatio', 'neighboursOfAllSurfaces', '-v7.3');
    end
end

