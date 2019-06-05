function [] = calculateSurfaceRatioMutant(basalLayer, apicalLayer, labelledImage, validCells)
%CALCULATESURFACERATIOMUNTANT Summary of this function goes here
%   We ought to get the surface ratio (RadiusBasal/RadiusApical). In this
%   case we want to get RadiusBasal as a cylindre amid two basal surfaces:
%   the one forming the outer gland (the maximum surface formed by the max
%   cell height of each cell) and an inner basal surface formed by the
%   minimum points of the cells intersecting.
    
    %% Get outer basal layer
    neighbours = calculateNeighbours3D(basalLayer, 2, basalLayer == 0);
    vertices3dBasal = getVertices3D(basalLayer, neighbours.neighbourhood);
    neighbours = calculateNeighbours3D(apicalLayer, 2, apicalLayer == 0);
    vertices3dApical = getVertices3D(apicalLayer, neighbours.neighbourhood);
    
    labelledImageToCalculateSR = zeros(size(labelledImage));
    for numCell = 1:max(labelledImage(:))
        verticesApical = vertcat(vertices3dApical.verticesPerCell{any(ismember(vertices3dApical.verticesConnectCells, numCell), 2)});
        verticesBasal = vertcat(vertices3dBasal.verticesPerCell{any(ismember(vertices3dBasal.verticesConnectCells, numCell), 2)});
        
        plotCell = alphaShape([verticesBasal; verticesApical], 1);
        plotCell.Alpha = criticalAlpha(plotCell, 'one-region')+10;
        
        indicesCell = find(labelledImage == numCell);
        [possibleIndicesX, possibleIndicesY, possibleIndicesZ] = ind2sub(size(labelledImage), indicesCell);
        pixelsIn = inShape(plotCell, possibleIndicesX, possibleIndicesY, possibleIndicesZ);
        labelledImageToCalculateSR(indicesCell(pixelsIn)) = numCell;
    end 
    
    labelledImageRealSR = double(imclose(labelledImageToCalculateSR > 0,strel('sphere',2)));
    labelledImageRealSR = labelledImageRealSR .* labelledImage;
    
    [apicalLayer] = getApicalFrom3DImage(imfill(labelledImageRealSR), labelledImageRealSR);
    [basalLayer] = getBasalFrom3DImage(labelledImage, lumenImage, tipValue, outsideGland);
    
    %% Create the cylinder between both surfaces
    load('D:\Pablo\LimeSeg_Pipeline\tmp\sr_info.mat')
    %TODO: REMOVE NO-VALID CELLS AREA
    plane2d = apicalLayer(:, :, 1) == -1;
    [allCoordinates2D_x, allCoordinates2D_y] = find(plane2d == 0);
    paint3D
    plane2d_outerBasal = plane2d;
    plane2d_innerBasal = plane2d;
    plane2d_apical = plane2d;
    for numCoordZ = unique(vertices3D(:, 3))'%1:size(apicalLayer, 3)
        %Outer basal layer
        inCoordinates = inShape(surfaceShapeOutter, allCoordinates2D_x, allCoordinates2D_y, repmat(numCoordZ, length(allCoordinates2D_x), 1));
        coordinatesIndices = sub2ind(size(apicalLayer(:, :, 1)), allCoordinates2D_x, allCoordinates2D_y);
        plane2d_outerBasal(coordinatesIndices(inCoordinates)) = 1;
        perim_outerBasal(:, :, numCoordZ) = bwperim(plane2d_outerBasal);
        [xOuterBasal, yOuterBasal] = find(perim_outerBasal(:, :, numCoordZ));
        
        %Inner basal layer
        inCoordinates = inShape(surfaceShapeInner, allCoordinates2D_x, allCoordinates2D_y, repmat(numCoordZ, length(allCoordinates2D_x), 1));
        coordinatesIndices = sub2ind(size(apicalLayer(:, :, 1)), allCoordinates2D_x, allCoordinates2D_y);
        plane2d_innerBasal(coordinatesIndices(inCoordinates)) = 1;
        perim_innerBasal(:, :, numCoordZ) = bwperim(plane2d_innerBasal);
        [xInnerBasal, yInnerBasal] = find(perim_innerBasal(:, :, numCoordZ));
        
        %Apical
        [xApical, yApical] = find(ismember(apicalLayer(:, :, numCoordZ), validCells));
        centroidApical = mean([xApical, yApical]);
        if isempty(xApical) == 0
            [distanceFromApicalToCentroidApical(numCoordZ)] = pdist2([xApical, yApical], centroidApical, 'euclidean', 'Largest', 1);
%             [distanceFromApicalToCentroidApical(numCoordZ)] = mean(pdist2([xApical, yApical], centroidApical));
        
            %Get distance between the centre of apical layer and the mid
            %basallayer to obtain the Radius basal
            [distanceFromInnerBasalToCentroidApical(numCoordZ)] = pdist2([xInnerBasal, yInnerBasal], centroidApical, 'euclidean', 'Largest', 1);
            [distanceFromOuterBasalToCentroidApical(numCoordZ)] = pdist2([xOuterBasal, yOuterBasal], centroidApical, 'euclidean', 'Largest', 1);
%             [distanceFromInnerBasalToCentroidApical] = mean(pdist2([xInnerBasal, yInnerBasal], centroidApical));
%             [distanceFromOuterBasalToCentroidApical] = mean(pdist2([xOuterBasal, yOuterBasal], centroidApical));
            distanceBasalToCentroid(numCoordZ) = mean([distanceFromInnerBasalToCentroidApical(numCoordZ), distanceFromOuterBasalToCentroidApical(numCoordZ)]);
            %basalPoints(numCoordZ, :) = [xInnerBasal(basalInnerPixelId(1)), yInnerBasal(basalInnerPixelId(1))];

            surfaceRatio(numCoordZ) = distanceBasalToCentroid(numCoordZ) / distanceFromApicalToCentroidApical(numCoordZ);
%             figure;imshow(double(plane2d_innerBasal)+double(plane2d_outerBasal)+apicalLayer(:, :, numCoordZ), parula(5))
%             hold on;
%             plot(centroidApical(2), centroidApical(1), 'r*')
        end
        %imshow(plane2d)
    end
    colours = colorcube(200);
    exportAsImageSequence(perim_outerBasal, 'outerbasalLayer', colours, 4)
    exportAsImageSequence(perim_innerBasal, 'innerbasalLayer', colours, 4)
    mean(surfaceRatio(round(size(apicalLayer, 3)/3):2*round(size(apicalLayer, 3)/3)))
    median(surfaceRatio(round(size(apicalLayer, 3)/3):2*round(size(apicalLayer, 3)/3)))
end

