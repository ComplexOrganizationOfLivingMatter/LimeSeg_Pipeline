function [] = calculateSurfaceRatioMuntant(basalLayer, apicalLayer, validCells)
%CALCULATESURFACERATIOMUNTANT Summary of this function goes here
%   We ought to get the surface ratio (RadiusBasal/RadiusApical). In this
%   case we want to get RadiusBasal as a cylindre amid two basal surfaces:
%   the one forming the outer gland (the maximum surface formed by the max
%   cell height of each cell) and an inner basal surface formed by the
%   minimum points of the cells intersecting.

%     %% Outer basal gland surface
%     convexImage3d = regionprops3(basalLayer>0, 'ConvexImage');
%     outterBoundary = convexImage3d.ConvexImage{1,1};

    %Another options is to get the maximum pixel distance and create the convex
    %image
    
%     basalPoints = zeros(max(basalLayer(:)), 3);
%     for numCell = 1:max(basalLayer(:))
%         centroidApical = regionprops3(apicalLayer == numCell, 'Centroid');
%         centroidApical = centroidApical.Centroid;
%         
%         indicesRealBasal = find(basalLayer == numCell);
%         [x, y, z] = ind2sub(size(basalLayer), indicesRealBasal);
%         [~, basalPixelId] = pdist2([x, y, z], centroidApical, 'euclidean', 'Largest', 1);
%         basalPoints(numCell, :) = [x(basalPixelId(1)), y(basalPixelId(1)), z(basalPixelId(1))];
%     end
%     
%     surfaceShapeOutter = alphaShape(basalPoints, 100);
%     %surfaceShape.Alpha = criticalAlpha(surfaceShape, 'one-region');
%     figure; plot(surfaceShapeOutter)
% 
%     
%     %% Inner basal surface
%     % Using vertices of the cells that are supposed to be the minima points
%     neighbours = getNeighboursFromFourProjectedPlanesFrom3Dgland(basalLayer, []);
%     neighbours = checkPairPointCloudDistanceCurateNeighbours(basalLayer, neighbours);
%     [verticesInfo] = getVertices3D(basalLayer, neighbours);
%     vertices3D = vertcat(verticesInfo.verticesPerCell{:});
% 
% %     indicesOfEdges = sub2ind(size(basalLayer), vertices3D(:, 1), vertices3D(:, 2), vertices3D(:, 3));
% %     innerSurfaceEdges = zeros(size(basalLayer));
% %     innerSurfaceEdges(indicesOfEdges) = 1;
% %     innerSurface3d = regionprops3(innerSurfaceEdges, {'ConvexImage', 'BoundingBox'});
% %     innerBasalSurface = innerSurface3d.ConvexImage{1,1};
% %     figure; paint3D(innerSurface3d.ConvexImage{1,1})
%     
%     surfaceShapeInner = alphaShape(vertices3D, 100);
%     %surfaceShape.Alpha = criticalAlpha(surfaceShape, 'one-region');
%     hold on; plot(surfaceShapeInner)
    
    %% Create the cylinder between both surfaces
    
    disp('creat');
    load('D:\Pablo\LimeSeg_Pipeline\tmp\sr_info.mat')
    %TODO: REMOVE NO-VALID CELLS AREA
    plane2d = apicalLayer(:, :, 1) == -1;
    [allCoordinates2D_x, allCoordinates2D_y] = find(plane2d == 0);
    
    plane2d_outerBasal = plane2d;
    plane2d_innerBasal = plane2d;
    plane2d_apical = plane2d;
    for numCoordZ = 1:size(apicalLayer, 3)
        %Outer basal layer
        inCoordinates = inShape(surfaceShapeOutter, allCoordinates2D_x, allCoordinates2D_y, repmat(numCoordZ, length(allCoordinates2D_x), 1));
        coordinatesIndices = sub2ind(size(apicalLayer(:, :, 1)), allCoordinates2D_x, allCoordinates2D_y);
        plane2d_outerBasal(coordinatesIndices(inCoordinates)) = 1;
        perim_outerBasal = bwperim(plane2d_outerBasal);
        [xOuterBasal, yOuterBasal] = find(perim_outerBasal);
        
        %Inner basal layer
        inCoordinates = inShape(surfaceShapeInner, allCoordinates2D_x, allCoordinates2D_y, repmat(numCoordZ, length(allCoordinates2D_x), 1));
        coordinatesIndices = sub2ind(size(apicalLayer(:, :, 1)), allCoordinates2D_x, allCoordinates2D_y);
        plane2d_innerBasal(coordinatesIndices(inCoordinates)) = 1;
        perim_innerBasal = bwperim(plane2d_innerBasal);
        [xInnerBasal, yInnerBasal] = find(perim_innerBasal);
        
        %Apical
        [xApical, yApical] = find(apicalLayer(:, :, numCoordZ));
        centroidApical = mean([xApical, yApical]);
        
        [~, basalInnerPixelId] = pdist2([xInnerBasal, yInnerBasal], centroidApical, 'euclidean', 'Largest', 1);
        basalPoints(numCoordZ, :) = [xInnerBasal(basalPixelId(1)), yInnerBasal(basalPixelId(1))];
        
        %Get distance between the centre of apical layer and the mid
        %basallayer to obtain the Radius basal
        
        %imshow(plane2d)
        plane2d(coordinatesIndices(inCoordinates)) = 0;
    end
end

