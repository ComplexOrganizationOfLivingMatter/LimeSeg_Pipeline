function [samiraTable, areaOfValidCells, rotationsOriginal] = unrollTube(img3d_original, outputDir, noValidCells, colours, apicalArea, apicalRotationsOriginal)
%UNROLLTUBE Summary of this function goes here
%   Detailed explanation goes here
    colours = vertcat([1 1 1], colours);
    
    %% Unroll
    pixelSizeThreshold = 10;
    
    previousSizeLabels = -1;
    if exist(fullfile(outputDir, 'final3DImg.mat'), 'file')
        load(fullfile(outputDir, 'final3DImg.mat'));
    else
%         %Option 1: too much neighbours
%         insideGland = imdilate(img3d_original>0, strel('sphere', 1));
%         edgePixels = find(insideGland & img3d_original==0);
%         img2Dilate = zeros(size(img3d_original));
%         sphereDilated = strel('sphere', 2);
%         neighbours3d = cell(max(img3d_original(:)), 1);
%         for numPixel = edgePixels'
%             img2Dilate(numPixel) = 1;
%             neighbours = unique(img3d_original(imdilate(img2Dilate, sphereDilated)>0));
%             if length(neighbours) > 2
%                 neighbours(neighbours==0) = [];
%                 for numNeighbour = neighbours'
%                     neighbours3d{numNeighbour} = unique([neighbours', neighbours3d{numNeighbour}]);
%                 end
%             end
%             img2Dilate(numPixel) = 0;
%         end
%         neighbours = neighbours3d;
% 
%         %Option 2: use the closest distance to a pixel
%     %     insideGland = imclose(img3d_original>0, strel('sphere', 3));
%     %     img3d_OnlyLateralWalls = img3d_original<=0 & insideGland;
%     %     img3d_OnlyLateralWalls = img3d_original .* imdilate(img3d_OnlyLateralWalls, strel('sphere', 1));
%         [x, y, z] = ind2sub(size(img3d_original), find(img3d_original>0));
%         cellsPointCloud = pointCloud([x, y, z], 'Intensity', img3d_original(img3d_original>0));
%         neighbours = cell(max(img3d_original(:)), 1);
%         for numCell = 1:max(img3d_original(:))
%             actualCellPerim = img3d_original==numCell;
%             actualPointCloud = select(cellsPointCloud, find(cellsPointCloud.Intensity ~= numCell));
%             for numPointPerim = find(actualCellPerim)'
%                 [x, y, z] = ind2sub(size(img3d_original), numPointPerim);
%                 [indices, ~] = findNearestNeighbors(actualPointCloud, [x, y, z], 5, 'MaxLeafChecks', 200, 'Sort', true);
%                 cellNeighbours = actualPointCloud.Intensity(indices);
%                 cellNeighbours(cellNeighbours == numCell) = [];
%                 if isempty(cellNeighbours) == 0
%                     neighbours{numCell} = unique([cellNeighbours', neighbours{numCell}]);
%                 end
%             end
%         end
% 
%         %Option 3: Original version
%         [neighbours] = calculateNeighbours3D(img3d_original); %Correct neighbours
%         neighbours = neighbours.neighbourhood;

        %Option 4: making projections
        if exist('apicalRotationsOriginal', 'var') == 0
            [img3d_original, rotationsOriginal] = rotateImg3(img3d_original);
            save(fullfile(outputDir, 'apicalRotationsOriginal.mat'), 'rotationsOriginal');
        else
            [img3d_original] = rotateImg3(img3d_original, apicalRotationsOriginal);
        end

        [allX,allY,allZ]=ind2sub(size(img3d_original),find(img3d_original>0));
        img3d_originalCropped = img3d_original(min(allX):max(allX),min(allY):max(allY),min(allZ):max(allZ));

        % Check which side is the longest
        img3d_closed = fill0sWithCells(img3d_original, imclose(img3d_original>0, strel('sphere', 3)) == 0);
        neighbours = getNeighboursFromFourProjectedPlanesFrom3Dgland(img3d_closed, colours);
        neighbours = checkPairPointCloudDistanceCurateNeighbours(img3d_closed, neighbours);


        sizeImg3d = size(img3d_originalCropped);
        [~, indices] = sort(sizeImg3d);
        img3d_original = permute(img3d_originalCropped, indices);
        img3d_original = fill0sWithCells(img3d_original, imclose(img3d_original>0, strel('sphere', 2)) == 0);
        tipValue = 21;
        img3d_original = addTipsImg3D(tipValue, img3d_original);
        [verticesInfo] = getVertices3D(img3d_original, neighbours);
        vertices3D = vertcat(verticesInfo.verticesPerCell{:});

    %     [colours] = exportAsImageSequence_WithVertices(img3d, outputDir, colours, -1, vertices3D);
    %     figure; paint3D(img3d_original, [], colours);
    %     for numVertex = 1:size(vertices3D, 1)
    %         hold on; plot3(vertices3D(numVertex, 1), vertices3D(numVertex, 2), vertices3D(numVertex, 3), 'rx')
    %     end

        vertices3D_Neighbours = verticesInfo.verticesConnectCells;
        vertices3D_Neighbours(cellfun(@isempty, verticesInfo.verticesPerCell), :) = [];
        cellNumNeighbours = cellfun(@length, neighbours);

        resizeImg = 1/4.06; %TODO: REFACTOR REGARDING THE MUTANT
        imgSize = round(size(img3d_original)/resizeImg);
        img3d = double(imresize3(img3d_original, imgSize, 'nearest'));
        validRegion_filled = imfill(double(imclose(img3d>0, strel('sphere', 20))), 26);
        validRegion = validRegion_filled & imerode(validRegion_filled, strel('sphere', 1))==0;
        img3d = fill0sWithCells(img3d .* double(validRegion), validRegion==0);
        vertices3D = round(vertices3D / resizeImg);
        mkdir(outputDir);
        save(fullfile(outputDir, 'final3DImg.mat'), 'img3d', 'vertices3D_Neighbours', 'vertices3D', 'cellNumNeighbours', 'neighbours', '-v7.3');
    end
    
    if exist(fullfile(outputDir, 'verticesInfo.mat'), 'file') == 0
    
        imgFinalCoordinates=cell(size(img3d,3),1);
        imgFinalCoordinates3x=cell(size(img3d,3),1);
        for coordZ = 1 : size(img3d,3)
            if sum(sum(img3d(:, :, coordZ) > 0)) < pixelSizeThreshold || sum(sum(img3d(:, :, coordZ))) < pixelSizeThreshold
                continue
            end
            %figure; imshow(img3d(:, :, coordZ)+2, colorcube)

            %% Remove pixels surrounding the boundary
            filledImage = imfill(double(img3d(:, :, coordZ)>0));
            %imshow(double(img3d(:, :, coordZ)>0))
            filledImage = bwareafilt(filledImage>0, 1, 4);
            finalPerimImage = bwperim(filledImage);
            solidityOfObjects = regionprops(filledImage, 'Solidity');
            
            solidityThreshold = 0.6;
            %Check if there is a hole
            if solidityOfObjects.Solidity < solidityThreshold
%                 convexPerimImage = regionprops(imclose(finalPerimImage, strel('disk', 5)), 'convexHull');
%                 convexPerimImage = convexPerimImage.ConvexHull;
%                 
%                 validRegion = zeros(size(finalPerimImage));
%                 [xq, yq] = find(validRegion==0);
%                 in = inpolygon(xq,yq, round(convexPerimImage(:, 2)), round(convexPerimImage(:, 1)));
%                 
%                 indicesInsideCell = sub2ind(size(finalPerimImage), xq, yq);
%                 
%                 validRegion(indicesInsideCell(in)) = 1;
%                 
%                 finalPerimImage = bwperim(validRegion);

                
                finalPerimImage = bwskel(filledImage);
                %fill0sWithCells(img3d(:, :, coordZ) ,validRegion);
                
                [X,Y] = meshgrid(1:size(filledImage,2), 1:size(filledImage,1));
                
                s = regionprops(filledImage, 'BoundingBox');
                
                bb = floor(s.BoundingBox); %// Could be floating point, so floor it
                cenx = bb(1) + (bb(3) / 2.0); %// Get the centre of the bounding box
                ceny = bb(2) + (bb(4) / 2.0);
                
                radi = max(bb(3), bb(4)) / 2; %// Find the best radius
                perimeterNew = ((X - cenx).^2 + (Y - ceny).^2) <= radi^2; %// Draw our circle and place in a temp. image
                perimeterNew = bwperim(perimeterNew); %// Add this circle on top of our output image
                %figure; imshow(perimeterNew)
            end

            %imshow(finalPerimImage)
            
            %% Obtaining the center of the cylinder
            [x, y] = find(finalPerimImage > 0);
            centroidCoordZ = mean([x, y], 1); % Centroid of each real Y of the cylinder
            centroidX = centroidCoordZ(1);
            centroidY = centroidCoordZ(2);

            [x, y] = find(finalPerimImage > 0);

            %% labelled mask
            if solidityOfObjects.Solidity < solidityThreshold
                centroidX = ceny;
                centroidY = cenx;
                
                [x_PerimeterNew, y_PerimeterNew] = find(perimeterNew > 0);
                angleLabelCoord_NewPerimeter = atan2(y_PerimeterNew - centroidY, x_PerimeterNew - centroidX);
                angleLabelCoord_NewPerimeter_Sorted = sort(angleLabelCoord_NewPerimeter);
                minDistance = abs(angleLabelCoord_NewPerimeter_Sorted(1) - angleLabelCoord_NewPerimeter_Sorted(2));
                img3dPerimFilled = fill0sWithCells(img3d(:, :, coordZ), filledImage==0);
                maskLabel=finalPerimImage.*img3dPerimFilled;
            else
                maskLabel=finalPerimImage.*img3d(:, :, coordZ);
            end
            %angles label coord regarding centroid
            angleLabelCoord = atan2(y - centroidY, x - centroidX);
            [angleLabelCoordSort, orderedIndices] = sort(angleLabelCoord);
            
            %% Completing the missing parts of the circle perim
            if solidityOfObjects.Solidity < solidityThreshold
                distanceToNextPoint = angleLabelCoordSort([2:end 1]) - angleLabelCoordSort;
                distanceToNextPoint(end) = distanceToNextPoint(end) + 6;
                if max(distanceToNextPoint) > minDistance*3
                    [~, positionsToFill] = max(distanceToNextPoint);
                    if positionsToFill+1 > length(distanceToNextPoint)
                        newAngles = angleLabelCoord_NewPerimeter_Sorted(angleLabelCoordSort(1) > angleLabelCoord_NewPerimeter_Sorted);
                        
                        angleLabelCoordSort = [newAngles; angleLabelCoordSort];
                        orderedIndices = [zeros(size(newAngles)); orderedIndices];
                    else
                        newAngles = angleLabelCoord_NewPerimeter_Sorted(angleLabelCoordSort(positionsToFill) < angleLabelCoord_NewPerimeter_Sorted & angleLabelCoordSort(positionsToFill+1) > angleLabelCoord_NewPerimeter_Sorted);
                        angleLabelCoordSort = [angleLabelCoordSort(1:positionsToFill); newAngles; angleLabelCoordSort(positionsToFill+1:end)];
                        orderedIndices = [orderedIndices(1:positionsToFill); zeros(size(newAngles)); orderedIndices(positionsToFill+1:end)];
                    end
                end
            end
            
            
            %% Assing label to pixels of perimeters
            %If a perimeter coordinate have no label pixels in a range of pi/45 radians, it label is 0
            orderedLabels = zeros(1,length(angleLabelCoordSort));
            for nCoord = 1:length(angleLabelCoordSort)
                if orderedIndices(nCoord) ~= 0
                    indicesClosest = sub2ind(size(maskLabel), x(orderedIndices(nCoord)), y(orderedIndices(nCoord)));
                    pixelLabel = maskLabel(indicesClosest);
                    orderedLabels(nCoord) = pixelLabel;
                else
                    orderedLabels(nCoord) = 0;
                end
            end
            
%             if abs(previousSizeLabels - length(angleLabelCoordSort)) > 150 && previousSizeLabels ~= -1
%                 orderedLabels = imresize(orderedLabels, [1 0.1*length(orderedLabels) + 0.9*previousSizeLabels], 'nearest');
%             end
%             previousSizeLabels = length(orderedLabels);
%             
            hold off;

            imgFinalCoordinates3x{coordZ} = repmat(orderedLabels, 1, 3);
            
            %If we find that a cell does not continue in the other side and
            %finish on one side, we put the same pixel on the other side.
            if orderedLabels(1) ~= orderedLabels(end)
                orderedLabels(end+1) = orderedLabels(1);
            end
            imgFinalCoordinates{coordZ} = orderedLabels;
        end

        %% Reconstruct deployed img
        ySize=max(cellfun(@length, imgFinalCoordinates3x));
        deployedImg3x = zeros(size(img3d,3),ySize);
        deployedImg = zeros(size(img3d,3),ySize);

        for coordZ = 1 : size(img3d,3)
            rowOfCoord3x = imgFinalCoordinates3x{coordZ};
            rowOfCoord = imgFinalCoordinates{coordZ};

            nEmptyPixels3x = 0;
            if length(rowOfCoord3x) < ySize
                nEmptyPixels3x = floor((ySize - length(rowOfCoord3x)) / 2);
                nEmptyPixels = floor((ySize - length(rowOfCoord)) / 2);
            end
            deployedImg3x(coordZ, 1 + nEmptyPixels3x : length(rowOfCoord3x) + nEmptyPixels3x) = rowOfCoord3x;
            deployedImg(coordZ, 1 + nEmptyPixels : length(rowOfCoord) + nEmptyPixels) = rowOfCoord;
        end

        %% Getting correct border cells, valid cells and no valid cells
         %cylindre2DImage = fillEmptySpacesByWatershed2D(deployedImg, imclose(deployedImg>0, strel('disk', 20)) == 0 , colours);
         cylindre2DImage = deployedImg;
         [wholeImage] = fillEmptySpacesByWatershed2D(deployedImg3x, imclose(deployedImg3x>0, strel('disk', 20)) == 0 , colours);

    %     figure;imshow(finalImage,colours)

        %% We only keep the cells in the middle
        relabelFinalImage = bwlabel(wholeImage,4);
        labelsFinal = unique(relabelFinalImage(deployedImg>0));
        midSectionImage = wholeImage;
        midSectionImage(~ismember(relabelFinalImage,labelsFinal))=0;

        %% We keep the valid cells from that middle image
        validCellsFinal  = setdiff(1:max(midSectionImage(:)), noValidCells);
        finalImageWithValidCells = ismember(midSectionImage, validCellsFinal).*midSectionImage;
    %     figure;imshow(finalImageWithValidCells,colours)

        h = figure ('units','normalized','outerposition',[0 0 1 1], 'visible', 'off');        
        imshow(midSectionImage+1, colours);
        hold on;
        midSectionNewLabels = bwlabel(midSectionImage, 4);
        centroids = regionprops(midSectionNewLabels, 'Centroid');
        centroids = round(vertcat(centroids.Centroid));

        ax = get(h, 'Children');
        set(ax,'Units','normalized')
        set(ax,'Position',[0 0 1 1])
        for numCentroid = 1:size(centroids, 1)
            labelSeed = midSectionImage(midSectionNewLabels == numCentroid);
            labelSeed = labelSeed(1);
            if mean(colours(labelSeed, :)) < 0.4
                text(ax, centroids(numCentroid, 1), centroids(numCentroid, 2), num2str(labelSeed), 'HorizontalAlignment', 'center', 'Color', 'white', 'FontSize', 6);
            else
                text(ax, centroids(numCentroid, 1), centroids(numCentroid, 2), num2str(labelSeed), 'HorizontalAlignment', 'center', 'FontSize', 6);
            end
        end
        h.InvertHardcopy = 'off';
        saveas(h, fullfile(outputDir, 'img_MidSection.tif'));
        imwrite(finalImageWithValidCells+1, colours, fullfile(outputDir, 'img_MidSection_ValidCells.tif'));
        imwrite(wholeImage+1, colours, fullfile(outputDir, 'img_WholeImage.tif'));
        imwrite(deployedImg+1, colours, fullfile(outputDir, 'img_original.tif'));

        %% Calculating surface ratio
        areaOfValidCells = sum(deployedImg(:)>0);

        if exist('apicalArea', 'var') == 0
            surfaceRatio = 1;
        else
            surfaceRatio = areaOfValidCells / apicalArea;
        end

        %% Connect vertices to obtain an image from the vertices
        newNeighbours2D = calculateNeighbours(deployedImg);
        newNeighbours2D_Checked = checkPairPointCloudDistanceCurateNeighbours(img3d, newNeighbours2D);

        newVertices2D = getVertices(deployedImg, newNeighbours2D_Checked);
        newVerticesNeighs2D = vertcat(newVertices2D.verticesConnectCells);
        newVerticesNeighs2D_empty = cellfun(@isempty, newVertices2D.verticesPerCell);
        newVerticesNeighs2D(newVerticesNeighs2D_empty, :) = [];
        newVertices2D = vertcat(newVertices2D.verticesPerCell{:});

        toGetBorderCells = cylindre2DImage+1;
        background = imopen(cylindre2DImage==0, strel('sphere', 1));
        toGetBorderCells(~background & cylindre2DImage == 0) = 0;
        backgroundNeighs = calculateNeighbours(toGetBorderCells);
        borderCells = backgroundNeighs{1} - 1;
        borderCells = intersect(validCellsFinal, borderCells);
        if exist('apicalArea', 'var') == 0
            nameOfSimulation = 'Apical';
        else
            nameOfSimulation = 'Basal';
        end
        
        occurrencesOfCells = tabulate(newVerticesNeighs2D(:));
        side_cells = occurrencesOfCells(:, 2);
        [polygon_distribution] = calculate_polygon_distribution(side_cells, validCellsFinal);
        polygon_distribution_T = cell2table(polygon_distribution(2:end, :));
        polygon_distribution_T.Properties.VariableNames = polygon_distribution(1, :);
        writetable(polygon_distribution_T, fullfile(outputDir, 'polygon_distribution.xls'))
        
        save(fullfile(outputDir,  'verticesInfo.mat'), 'cylindre2DImage', 'newVerticesNeighs2D', 'newVertices2D', 'centroids', 'validCellsFinal', 'borderCells', 'surfaceRatio', 'outputDir', 'nameOfSimulation', 'areaOfValidCells');
        
        save(fullfile(outputDir, 'allInfo.mat'), 'deployedImg', 'deployedImg3x', 'wholeImage', 'polygon_distribution');
    else
        load(fullfile(outputDir, 'allInfo.mat'));
        outputDirOri = outputDir;
        load(fullfile(outputDir, 'verticesInfo.mat'));
        outputDir = outputDirOri;
    end

    if exist(fullfile(outputDir, 'samiraTable.mat'), 'file') == 0
        outputDirSplitted = strsplit(outputDir, 'Results');
        load(fullfile(outputDirSplitted{1}, 'Results', 'valid_cells.mat'))
        samiraTable = connectVerticesOf2D(deployedImg, newVerticesNeighs2D, newVertices2D, centroids, validCells, borderCells, surfaceRatio, outputDir, nameOfSimulation, deployedImg3x, img3d);
        save(fullfile(outputDir, 'samiraTable.mat'), 'samiraTable');
    else
        load(fullfile(outputDir, 'samiraTable.mat'));
    end
end

