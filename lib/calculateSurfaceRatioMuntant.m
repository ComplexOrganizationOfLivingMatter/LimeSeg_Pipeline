function [] = calculateSurfaceRatioMuntant(img3d)
%CALCULATESURFACERATIOMUNTANT Summary of this function goes here
%   We ought to get the surface ratio (RadiusBasal/RadiusApical). In this
%   case we want to get RadiusBasal as a cylindre amid two basal surfaces:
%   the one forming the outer gland (the maximum surface formed by the max
%   cell height of each cell) and an inner basal surface formed by the
%   minimum points of the cells intersecting.

%% Outer basal gland surface
convexImage3d = regionprops3(img3d>0, 'ConvexImage');
outterBoundary = convexImage3d.ConvexImage{1,1};

%Another options is to get the maximum pixel distance and create the convex
%image

%% Inner basal surface
% Using vertices of the cells that are supposed to be the minima points
neighbours = getNeighboursFromFourProjectedPlanesFrom3Dgland(img3d, []);
neighbours = checkPairPointCloudDistanceCurateNeighbours(img3d, neighbours);
[verticesInfo] = getVertices3D(img3d, neighbours);
vertices3D = vertcat(verticesInfo.verticesPerCell{:});

end

