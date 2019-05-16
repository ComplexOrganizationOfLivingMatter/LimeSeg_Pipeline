function [] = calculateSurfaceRatioMuntant(basalLayer, apicalLayer)
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
    
    basalPoints = zeros(max(basalLayer(:)), 3);
    for numCell = 1:max(basalLayer(:))
        centroidApical = regionprops3(apicalLayer == numCell, 'Centroid');
        centroidApical = centroidApical.Centroid;
        
        indicesRealBasal = find(basalLayer == numCell);
        [x, y, z] = ind2sub(size(basalLayer), indicesRealBasal);
        [~, basalPixelId] = pdist2([x, y, z], centroidApical, 'euclidean', 'Largest', 1);
        basalPoints(numCell, :) = [x(basalPixelId(1)), y(basalPixelId(1)), z(basalPixelId(1))];
    end
    
    surfaceShapeOutter = alphaShape(basalPoints, 100);
    %surfaceShape.Alpha = criticalAlpha(surfaceShape, 'one-region');
    figure; plot(surfaceShapeOutter)

    
    %% Inner basal surface
    % Using vertices of the cells that are supposed to be the minima points
    neighbours = getNeighboursFromFourProjectedPlanesFrom3Dgland(basalLayer, []);
    neighbours = checkPairPointCloudDistanceCurateNeighbours(basalLayer, neighbours);
    [verticesInfo] = getVertices3D(basalLayer, neighbours);
    vertices3D = vertcat(verticesInfo.verticesPerCell{:});

%     indicesOfEdges = sub2ind(size(basalLayer), vertices3D(:, 1), vertices3D(:, 2), vertices3D(:, 3));
%     innerSurfaceEdges = zeros(size(basalLayer));
%     innerSurfaceEdges(indicesOfEdges) = 1;
%     innerSurface3d = regionprops3(innerSurfaceEdges, {'ConvexImage', 'BoundingBox'});
%     innerBasalSurface = innerSurface3d.ConvexImage{1,1};
%     figure; paint3D(innerSurface3d.ConvexImage{1,1})
    
    surfaceShapeInner = alphaShape(vertices3D, 100);
    %surfaceShape.Alpha = criticalAlpha(surfaceShape, 'one-region');
    hold on; plot(surfaceShapeInner)
    
    %% Create the cylinder between both surfaces
    
    disp('creat');
end

