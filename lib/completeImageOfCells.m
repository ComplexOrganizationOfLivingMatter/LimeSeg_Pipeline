function [labelMask] = completeImageOfCells(labelMask, invalidRegion)
%COMPLETEIMAGEOFCELLS Summary of this function goes here
%   Detailed explanation goes here

    edgePixels = find(labelMask == 0 & invalidRegion == 0);
    [xEdge, yEdge, zEdge] = ind2sub(size(labelMask), edgePixels);
    pixelsIndices = find(labelMask>0);
    [x, y, z] = ind2sub(size(labelMask), pixelsIndices);
    [~, indices] = pdist2([x, y, z], [xEdge, yEdge, zEdge], 'euclidean', 'Smallest', 1);
    labelMask(edgePixels) = labelMask(pixelsIndices(indices));
    %figure; paint3D(labelMask);
end

