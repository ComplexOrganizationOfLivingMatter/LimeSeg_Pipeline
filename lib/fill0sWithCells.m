function [labelMask] = fill0sWithCells(labelMask, invalidRegion)
%COMPLETEIMAGEOFCELLS Summary of this function goes here
%   Detailed explanation goes here

    missingRegions = labelMask == 0 & invalidRegion == 0;
    edgePixels = find(missingRegions);
    [xEdge, yEdge, zEdge] = ind2sub(size(labelMask), edgePixels);
    pixelsIndices = find((imdilate(missingRegions, strel('sphere', 4)).*labelMask)>0);
    [x, y, z] = ind2sub(size(labelMask), pixelsIndices);
    [~, indices] = pdist2([x, y, z], [xEdge, yEdge, zEdge], 'euclidean', 'Smallest', 1);
    labelMask(edgePixels) = labelMask(pixelsIndices(indices));
    %figure; paint3D(labelMask);
end

