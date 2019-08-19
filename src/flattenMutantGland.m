function [labelledImageRealSR] = flattenMutantGland(apicalLayer, basalLayer, labelledImage, lumenImage)
%FLATTENMUTANTGLAND Obtain the gland
%   Detailed explanation goes here

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
    labelledImageRealSR(lumenImage>0) = 0;
end

