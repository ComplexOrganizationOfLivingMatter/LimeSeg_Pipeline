function [labelledImageLayer, verticesLayer] = obtainLayerByVertices(layer, labelledImage, validCells)
%OBTAINLAYERBYVERTICES Summary of this function goes here
%   Detailed explanation goes here
%% Get inner basal layer
    verticesLayer = [];
    for numZ = 1:size(layer, 3)
        neighs_real = calculateNeighbours(layer(:, :, numZ), 2);
        for numCell = 1:length(neighs_real)
            if isempty(neighs_real{numCell}) == 0
                neighbours = neighs_real{numCell};
                numCellDilated = imdilate(layer(:, :, numZ) == numCell, strel('disk', 2));
                for actualNeighbour = neighbours'
                    numCellDilated_Neighbour = imdilate(layer(:, :, numZ) == actualNeighbour, strel('disk', 2));
                    [x, y] = find(numCellDilated & numCellDilated_Neighbour);
                    verticesLayer = [verticesLayer; mean(x) mean(y) numZ];
                end
            end
        end
    end
    verticesLayerUnique = unique(round(verticesLayer), 'rows');
    
    surfaceShape = alphaShape(verticesLayerUnique, 1);
    surfaceShape.Alpha = criticalAlpha(surfaceShape, 'one-region')+5;
%     figure; plot(surfaceShapeInner)

    possibleIndices = find(labelledImage > 0);
    [possibleIndicesX, possibleIndicesY, possibleIndicesZ] = ind2sub(size(labelledImage), possibleIndices);
    pixelsIn = inShape(surfaceShape, possibleIndicesX, possibleIndicesY, possibleIndicesZ);

    labelledImageLayer = labelledImage;
    labelledImageLayer(possibleIndices(pixelsIn==0)) = 0;
end

