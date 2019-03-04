function divideObjectInSurfaceRatios(selpath)
%DIVIDEOBJECTINSURFACERATIOS Summary of this function goes here
%   Detailed explanation goes here
    

    %% Loading variables
    load(fullfile(selpath, '3d_layers_info.mat'));
    load(fullfile(selpath, 'valid_cells.mat'));
    
    load(fullfile(selpath, 'apical', 'verticesInfo.mat'), 'newVerticesNeighs2D');
    apicalNeighs = newVerticesNeighs2D;
    apicalRealNeighs = arrayfun(@(x) sum(any(ismember(apicalNeighs, x), 2)), 1:max(validCells));

    load(fullfile(selpath, 'basal', 'verticesInfo.mat'), 'newVerticesNeighs2D');
    basalNeighs = newVerticesNeighs2D;
    basalRealNeighs = arrayfun(@(x) sum(any(ismember(basalNeighs, x), 2)), 1:max(validCells));
    
    obj_img = labelledImage;
    startingSurface = basalLayer; 
    endSurface = apicalLayer;
    
    %% Divide object in surface ratios
    apical3dInfo = calculateNeighbours3D(endSurface, 2);
    neighboursApical = checkPairPointCloudDistanceCurateNeighbours(endSurface, apical3dInfo.neighbourhood);
    numNeighsApical = cellfun(@length, neighboursApical);
    sum(numNeighsApical(validCells) ~= apicalRealNeighs(validCells)')
    
    basal3dInfo = calculateNeighbours3D(startingSurface, 2);
    neighboursBasal_init = checkPairPointCloudDistanceCurateNeighbours(startingSurface, basal3dInfo.neighbourhood);
    numNeighsBasal = cellfun(@length, neighboursBasal_init);
    sum(numNeighsBasal(validCells) ~= basalRealNeighs(validCells)')
    
    neighbours_data = table(neighboursApical', neighboursBasal_init');
    neighbours_data.Properties.VariableNames = {'Apical','Basal'};
    [~, apicoBasal_SurfaceRatio] = calculate_CellularFeatures(neighbours_data, neighboursApical', neighboursBasal_init', endSurface, startingSurface, obj_img, noValidCells, validCells, [], []);
    
    %% Split in 10 pieces
    totalPartitions = 10;
    initialPartitions = (1:(totalPartitions-1))/totalPartitions;
    
    %% Split with given surface ratios
%     definedSurfaceRatio = 1./(0.9:-0.1:0.3);
%     totalPartitions = length(definedSurfaceRatio)+1;
%     initialPartitions = (definedSurfaceRatio - 1) / (apicoBasal_SurfaceRatio - 1);
    
    realSurfaceRatio = initialPartitions * (apicoBasal_SurfaceRatio - 1) + 1;
    
    imageOfSurfaceRatios = cell(totalPartitions, 1);
    imageOfSurfaceRatios(:) = {zeros(size(obj_img))};
    for numCell = unique([validCells, noValidCells])
        [xStarting, yStarting, zStarting] = ind2sub(size(startingSurface), find(startingSurface == numCell));
        [xEnd, yEnd, zEnd] = ind2sub(size(endSurface), find(endSurface == numCell));

        [allXs, allYs, allZs] = ind2sub(size(obj_img), find(obj_img == numCell));

        allPixels = [allXs, allYs, allZs];
        startingPixels = [xStarting, yStarting, zStarting];
        endPixels = [xEnd, yEnd, zEnd];
        [distanceEndingStarting] = pdist2(endPixels, startingPixels);
        %[distanceStartingAllPixels] = pdist2(allPixels, startingPixels);
        [distanceEndingAllPixels] = pdist2(allPixels, endPixels);

        %surfaceRatioDistance = mean(min(distanceEndingStarting, [], 2));
        surfaceRatioDistance = mean(distanceEndingStarting(:));
        
        partitions = surfaceRatioDistance * initialPartitions;

        
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

                indicesCell = sub2ind(size(obj_img), x, y, z);

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
    
    imageOfSurfaceRatios(:, 2) = num2cell([1 realSurfaceRatio]);
    imageOfSurfaceRatios{totalPartitions+1, 1} = obj_img;
    imageOfSurfaceRatios{totalPartitions+1, 2} = apicoBasal_SurfaceRatio;
    
    for numPartition = 1:(totalPartitions+1)
        if numPartition > 1
            tipAdded = 26;
            initialImage_Tips = double(addTipsImg3D(tipAdded, imageOfSurfaceRatios{numPartition, 1}));

            initialBasalImage = getBasalFrom3DImage(initialImage_Tips, addTipsImg3D(tipAdded, lumenImage)>0, 4);
            
            initialBasalImage_closed = imclose(double(initialBasalImage>0), strel('sphere', tipAdded-1));
            initialBasalImage_filled = imfill(initialBasalImage_closed);
            
%             initialBasalImage_filled = permute(initialBasalImage_filled, [1 3 2]);
%             figure; imshow(initialBasalImage_filled(:, :, 100))
%             initialBasalImage_Tips = permute(initialBasalImage_Tips, [1 3 2]);
            
            innerRegion = initialBasalImage_filled - imerode(initialBasalImage_filled, strel('sphere', 1));
%             test = permute(innerRegion, [1 3 2]);
%             figure; imshow(test(:, :, 100))

            
            finalBasalImage = fill0sWithCells(double(initialBasalImage_Tips), innerRegion == 0);
            %figure; paint3D(finalBasalImage, [], colours);
            
            [imageOfSurfaceRatios{numPartition, 3}] = finalBasalImage(tipAdded+1:(size(finalBasalImage, 1) - tipAdded), tipAdded+1:(size(finalBasalImage, 2) - tipAdded), tipAdded+1:(size(finalBasalImage, 3) - tipAdded));
        else
            imageOfSurfaceRatios{numPartition, 3} = endSurface;
        end
        basal3dInfo = calculateNeighbours3D(imageOfSurfaceRatios{numPartition, 3}, 2);
        neighboursBasal = checkPairPointCloudDistanceCurateNeighbours(imageOfSurfaceRatios{numPartition, 3}, basal3dInfo.neighbourhood);
        neighbours_data = table(neighboursApical', neighboursBasal');
        neighbours_data.Properties.VariableNames = {'Apical','Basal'};
        [imageOfSurfaceRatios{numPartition, 4}, meanSurfaceRatio(numPartition)] = calculate_CellularFeatures(neighbours_data, neighboursApical', neighboursBasal', endSurface, imageOfSurfaceRatios{numPartition, 3}, imageOfSurfaceRatios{numPartition, 1}, noValidCells, validCells, [], []);
        
        neighbours_data = table(neighboursBasal_init', neighboursBasal');
        neighbours_data.Properties.VariableNames = {'Apical','Basal'};
        [imageOfSurfaceRatios{numPartition, 5}, ~] = calculate_CellularFeatures(neighbours_data, neighboursBasal_init', neighboursBasal', startingSurface, imageOfSurfaceRatios{numPartition, 3}, imageOfSurfaceRatios{numPartition, 1}, noValidCells, validCells, [], []);
        neighbours{numPartition} = neighboursBasal;
        %figure; paint3D( imageOfSurfaceRatios{numPartition, 1}, [], colours);
        h = figure('Visible', 'off'); paint3D( ismember(imageOfSurfaceRatios{numPartition, 3}, validCells) .* imageOfSurfaceRatios{numPartition, 3}, [], colours);
        mkdir(fullfile(selpath, 'dividedGland'));
        print(h, fullfile(selpath, 'dividedGland' , ['gland_SR' num2str(meanSurfaceRatio(numPartition)), '.tif']),'-dtiff','-r300')
    end
    close all
    
    infoPerSurfaceRatio = imageOfSurfaceRatios;
    save(fullfile(selpath, 'dividedGland', 'glandDividedInSurfaceRatios.mat'), 'infoPerSurfaceRatio', 'neighbours');
end