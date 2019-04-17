function [labelMask] = fill0sWithCells(labelMask, img3dComplete, invalidRegion)
%COMPLETEIMAGEOFCELLS Summary of this function goes here
%   Detailed explanation goes here

    missingRegions = labelMask == 0 & invalidRegion == 0;
    edgePixels = find(missingRegions);
    labelMask(missingRegions) = img3dComplete(missingRegions);

    missingRegions = labelMask == 0 & invalidRegion == 0;
    edgePixels = find(missingRegions);
    if isempty(edgePixels) == 0
        [xEdge, yEdge, zEdge] = ind2sub(size(labelMask), edgePixels);
        pixelsIndices = find((imdilate(missingRegions, strel('sphere', 4)).*labelMask)>0);
        [x, y, z] = ind2sub(size(labelMask), pixelsIndices);

        numPartitions = 1000;
        indices = cell(numPartitions, 1);
        distances = cell(numPartitions, 1);
        partialPxs = floor(length(x)/numPartitions);
        for nPart = 1 : numPartitions
            subIndCoord = (1 + (nPart-1) * partialPxs) : (nPart * partialPxs);
            if nPart == numPartitions
                subIndCoord = (1 + (nPart-1) * partialPxs) : length(x);
            end
            [distances{nPart}, indices_nPart] = pdist2([x(subIndCoord), y(subIndCoord), z(subIndCoord)], [xEdge, yEdge, zEdge], 'euclidean', 'Smallest', 1);
            indices{nPart} = subIndCoord(indices_nPart);
        end
        distances_all = vertcat(distances{:});
        indices_all = vertcat(indices{:});

        [~, indices_Partitions] = min(distances_all);

        for numEdgePixel = 1:length(indices_Partitions)
            labelMask(edgePixels(numEdgePixel)) = labelMask(pixelsIndices(indices_all(indices_Partitions(numEdgePixel), numEdgePixel)));
        end
    end
    %figure; paint3D(labelMask);
end

