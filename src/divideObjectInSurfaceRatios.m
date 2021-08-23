function [infoPerSurfaceRatio] = divideObjectInSurfaceRatios(selpath, testing)
%DIVIDEOBJECTINSURFACERATIOS Divide cylindrical object into different
%surface ratios
%   Obtain the cells’ neighborhood relation of salivary glands for
%   different values of the radial expansion.
    infoPerSurfaceRatio = [];
    if (exist(fullfile(selpath, 'dividedGland', 'glandDividedInSurfaceRatios.mat'), 'file') == 0 && exist(fullfile(selpath, 'unrolledGlands', 'gland_SR_1', 'verticesInfo.mat'), 'file')>0) || exist('testing', 'var')
        %% Loading variables
        variablesOfFile = who('-file', fullfile(selpath, 'realSize3dLayers.mat'));
        if any(cellfun(@(x) isequal('labelledImage_realSizeFlatten', x), variablesOfFile))
            load(fullfile(selpath, 'realSize3dLayers.mat'), 'labelledImage_realSizeFlatten', 'lumenImage_realSize');
            labelledImage_realSize = labelledImage_realSizeFlatten;
            clearvars labelledImage_realSizeFlatten
        else
            load(fullfile(selpath, 'realSize3dLayers.mat'), 'labelledImage_realSize', 'lumenImage_realSize');
        end
        load(fullfile(selpath, 'valid_cells.mat'));
        load(fullfile(selpath, '3d_layers_info.mat'),'colours');
       
        %% Obtain layers on its real 3D size
        basalLayer = getBasalFrom3DImage(labelledImage_realSize, lumenImage_realSize, 0, labelledImage_realSize == 0 & lumenImage_realSize == 0);
        [apicalLayer] = getApicalFrom3DImage(lumenImage_realSize, labelledImage_realSize);
        lumenImage = lumenImage_realSize>0;
        startingSurface = basalLayer;
        endSurface = apicalLayer;
        
        %% Divide object in surface ratios
        apical3dInfo = calculateNeighbours3D(endSurface, 2);
        neighboursApical = checkPairPointCloudDistanceCurateNeighbours(endSurface, apical3dInfo.neighbourhood);
        numNeighsApical = cellfun(@length, neighboursApical);
        
        basal3dInfo = calculateNeighbours3D(startingSurface, 2);
        neighboursBasal_init = checkPairPointCloudDistanceCurateNeighbours(startingSurface, basal3dInfo.neighbourhood);
        numNeighsBasal = cellfun(@length, neighboursBasal_init);
        
        neighbours_data = table(neighboursApical', neighboursBasal_init');
        neighbours_data.Properties.VariableNames = {'Apical','Basal'};

        %%  Calculate basal surface ratio
        apical_area_cells=cell2mat(struct2cell(regionprops(apicalLayer,'Area'))).';
        basal_area_cells=cell2mat(struct2cell(regionprops(basalLayer,'Area'))).';
        if length(apical_area_cells) > length(basal_area_cells)
            basal_area_cells(length(apical_area_cells)) = 0;
        elseif length(apical_area_cells) < length(basal_area_cells)
            apical_area_cells(length(basal_area_cells)) = 0;
        end
        %meanSurfaceRatio = mean(surfaceRatioValidCells);
        apicoBasal_SurfaceRatio = sum(basal_area_cells(validCells)) / sum(apical_area_cells(validCells));
        
        %% Split in 10 pieces
        %             totalPartitions = 10;
        %             initialPartitions = (1:(totalPartitions-1))/totalPartitions;
        
        %% Split with given surface ratios
        initialPartitions = ((1.5:0.5:apicoBasal_SurfaceRatio) - 1) / (apicoBasal_SurfaceRatio - 1);
        totalPartitions = length(initialPartitions)+1;
        
        realSurfaceRatio = 1:0.5:apicoBasal_SurfaceRatio;
        
        imageOfSurfaceRatios = cell(totalPartitions, 1);
        imageOfSurfaceRatios(:) = {zeros(size(labelledImage_realSize))};
        
               
        for numCell = unique([validCells, noValidCells])
            % First, We calculated the cell height by estimating the
            % distance between the average voxel positions of the
            % apical surface with respect to the average voxel
            % positions of its basal surface
            [xStarting, yStarting, zStarting] = ind2sub(size(startingSurface), find(startingSurface == numCell));
            [xEnd, yEnd, zEnd] = ind2sub(size(endSurface), find(endSurface == numCell));
            
            [allXs, allYs, allZs] = ind2sub(size(labelledImage_realSize), find(labelledImage_realSize == numCell));
            
            allPixels = [allXs, allYs, allZs];
            startingPixels = [xStarting, yStarting, zStarting];
            endPixels = [xEnd, yEnd, zEnd];
            %[distanceEndingStarting] = pdist2(endPixels, startingPixels);
            %[distanceStartingAllPixels] = pdist2(allPixels, startingPixels);
            [distanceEndingAllPixels] = pdist2(allPixels, endPixels);
            
            %surfaceRatioDistance = mean(min(distanceEndingStarting, [], 2));
            surfaceRatioDistance =  pdist2(mean(endPixels), mean(startingPixels));
            %                 surfaceRatioDistance = mean(distanceEndingStarting(:));
            
            partitions = surfaceRatioDistance * initialPartitions;
            
            clearvars endPixels startingPixels allXs allYs allZs xStarting yStarting zStarting
            
            imageOfSurfaceRatios{1, 1} = endSurface;
            imageOfSurfaceRatios{1, 2} = 1;
            for numPartition = 1:(totalPartitions-1)
                %             upperBoundStarting = ((ceil(partitions(totalPartitions - numPartition)*roundingFactor)/roundingFactor)+(1/roundingFactor)) >= distanceStartingAllPixels;
                %             lowerBoundStarting = ((floor(partitions(totalPartitions - numPartition)*roundingFactor)/roundingFactor)-(1/roundingFactor)) <= distanceStartingAllPixels;
                %
                %             upperBoundEnd = ((ceil(partitions(numPartition)*roundingFactor)/roundingFactor)+(1/roundingFactor)) >= distanceEndingAllPixels;
                %             lowerBoundEnd = ((floor(partitions(numPartition)*roundingFactor)/roundingFactor)-(1/roundingFactor)) <= distanceEndingAllPixels;
                %
                %             pixelsOfCurrentPartitionSurfaceRatioFromStarting = any(upperBoundStarting & lowerBoundStarting, 2);
                %             pixelsOfCurrentPartitionSurfaceRatioFromEnd = any(upperBoundEnd & lowerBoundEnd, 2);
                %             pixelsOfSR = allPixels(pixelsOfCurrentPartitionSurfaceRatioFromStarting & pixelsOfCurrentPartitionSurfaceRatioFromEnd, :);
                
                %pixelsCloserToEndSurface = partitions(numPartition) >= min(distanceEndingAllPixels, [], 2);
                pixelsCloserToEndSurface = partitions(numPartition) > distanceEndingAllPixels;
                pixelsOfSR = allPixels(any(pixelsCloserToEndSurface, 2), :);
                
                if isempty(pixelsOfSR) == 0
                    x = pixelsOfSR(:, 1);
                    y = pixelsOfSR(:, 2);
                    z = pixelsOfSR(:, 3);
                    
                    indicesCell = sub2ind(size(labelledImage_realSize), x, y, z);
                    
                    actualPartition = imageOfSurfaceRatios{numPartition+1};
                    
                    actualPartition(indicesCell) = numCell;
                    
                    imageOfSurfaceRatios{numPartition+1} = actualPartition;
                    
                    %                 figure; paint3D(actualSurface, [], colours(2:end, :));
                    %                 hold on; paint3D(startingSurface == 1, [], colours);
                    %
                    %                 for numIndex = 1:length(x)
                    %                    plot3(x(numIndex), y(numIndex), z(numIndex), '*r');
                    %                 end
                end
            end
        end
        
        clearvars allPixels
        
        imageOfSurfaceRatios(:, 2) = num2cell(realSurfaceRatio);
        imageOfSurfaceRatios{totalPartitions+1, 1} = uint16(labelledImage_realSize);
        imageOfSurfaceRatios{totalPartitions+1, 2} = uint16(apicoBasal_SurfaceRatio);
        
        %% Generate again the images from the previous information
        for numPartition = 1:(totalPartitions+1)
            if numPartition > 1
                initialImage = imageOfSurfaceRatios{numPartition, 1};
                %initialBasalImage = getBasalFrom3DImage(initialImage, lumenImage>0, 4);
                
                if numPartition > 2 || contains(selpath, 'e-cadh')==0
                    basalImage_closed_initial = imclose(initialImage>0, strel('sphere', 2));
                else
                    basalImage_closed_initial = imclose(imdilate(initialImage>0, strel('sphere', 1)), strel('sphere', 5));
                end
                
                basalImage_closed = basalImage_closed_initial;
                [~, y, ~] = ind2sub(size(initialImage),find(initialImage>0));
                
                downSide = initialImage(:, min(y), :);
                upSide = initialImage(:, max(y), :);
                if sum(downSide(:)>0) > sum(upSide(:)>0)
                    basalImage_closed(:, 1:(min(y)+3), :) = 1;
                else
                    basalImage_closed(:, (max(y)-3):end, :) = 1;
                end
                basalImage_closed(lumenImage) = 1;
                basalImage_filled = imfill(double(basalImage_closed), 18, 'holes');
                if sum(downSide(:)>0) > sum(upSide(:)>0)
                    basalImage_filled(:, 1:(min(y)+3), :) = basalImage_closed_initial(:, 1:(min(y)+3), :);
                else
                    basalImage_filled(:, (max(y)-3):end, :) = basalImage_closed_initial(:, (max(y)-3):end, :);
                end
                
                se = strel('sphere', 1);
                finalObjectEroded = imerode(basalImage_filled, se);
                basalLayer = basalImage_filled - finalObjectEroded;
                
                basalLayer(:, :, end) = initialImage(:, :, end);
                basalLayer(:, :, 1) = initialImage(:, :, 1);
                
                
                if sum(downSide(:)>0) > sum(upSide(:)>0)
                    basalLayer(:, 1:(min(y)+3), :) = 0;
                else
                    basalLayer(:, (max(y)-3):end, :) = 0;
                end
                
                regionsFound = regionprops3(basalLayer>0, {'Volume', 'VoxelIdxList'});
                if size(regionsFound, 1) > 1
                    [~, biggestRegion] = max(regionsFound.Volume);
                    smallerRegions = setdiff(1:size(regionsFound, 1), biggestRegion);
                    badIds = vertcat(regionsFound.VoxelIdxList{smallerRegions});
                    basalLayer(badIds) = 0;
                end
                %basalLayer(imfill(lumenImage, 'holes')) = 0;
                
                [imageOfSurfaceRatios{numPartition, 3}] = uint16(fill0sWithCells(uint16(initialImage).*uint16(basalLayer), labelledImage_realSize, basalLayer == 0));
                %unrollTube(imageOfSurfaceRatios{numPartition, 3}, fullfile(selpath, ['gland_SR_' num2str(imageOfSurfaceRatios{numPartition, 2})]), noValidCells, colours, 1);
            else
                imageOfSurfaceRatios{numPartition, 3} = endSurface;
            end
            
            %% Calculate and export information of each concentric layer
            basal3dInfo = calculateNeighbours3D(imageOfSurfaceRatios{numPartition, 3}, 2);
            neighboursBasal = checkPairPointCloudDistanceCurateNeighbours(imageOfSurfaceRatios{numPartition, 3}, basal3dInfo.neighbourhood);
            [imageOfSurfaceRatios{numPartition, 4}, meanSurfaceRatio(numPartition)] = calculate_CellularFeatures(neighboursApical', neighboursBasal', endSurface, imageOfSurfaceRatios{numPartition, 3}, imageOfSurfaceRatios{numPartition, 1}, noValidCells, validCells, [], []);
            
            [imageOfSurfaceRatios{numPartition, 5}, ~] = calculate_CellularFeatures(neighboursBasal_init', neighboursBasal', startingSurface, imageOfSurfaceRatios{numPartition, 3}, imageOfSurfaceRatios{numPartition, 1}, noValidCells, validCells, [], []);
            %figure; paint3D( imageOfSurfaceRatios{numPartition, 1}, [], colours);
            if ~exist('testing', 'var')
                h = figure('Visible', 'off');
                paint3D( ismember(imageOfSurfaceRatios{numPartition, 3}, validCells) .* double(imageOfSurfaceRatios{numPartition, 3}), [], colours, 2);
                mkdir(fullfile(selpath, 'dividedGland'));
                savefig(h, fullfile(selpath, 'dividedGland' , ['gland_SR' num2str(meanSurfaceRatio(numPartition)), '.fig']))
                print(h, fullfile(selpath, 'dividedGland' , ['gland_SR' num2str(meanSurfaceRatio(numPartition)), '.jpeg']),'-djpeg','-r300')
            end
        end
        close all
        
        infoPerSurfaceRatio = imageOfSurfaceRatios;
        if ~exist('testing', 'var')
            save(fullfile(selpath, 'dividedGland', 'glandDividedInSurfaceRatios.mat'), 'infoPerSurfaceRatio', '-v7.3');
        end
    end
end