function [labelledImageRealSR] = flattenMutantGland(apicalLayer, basalLayer, labelledImage, lumenImage)
%FLATTENMUTANTGLAND Obtain the gland
%   Detailed explanation goes here

    neighbours = calculateNeighbours3D(basalLayer, 2, basalLayer == 0);
    vertices3dBasal = getVertices3D(basalLayer, neighbours.neighbourhood);
    neighbours = calculateNeighbours3D(apicalLayer, 2, apicalLayer == 0);
    vertices3dApical = getVertices3D(apicalLayer, neighbours.neighbourhood);
    
    labelledImageToCalculateSR = zeros(size(labelledImage));
    uniqCells = unique(labelledImage(:))';
    for numCell = uniqCells(2:end)
        verticesApical = vertcat(vertices3dApical.verticesPerCell{any(ismember(vertices3dApical.verticesConnectCells, numCell), 2)});
        verticesBasal = vertcat(vertices3dBasal.verticesPerCell{any(ismember(vertices3dBasal.verticesConnectCells, numCell), 2)});
        try 
            plotCell = alphaShape([verticesBasal; verticesApical], 1);
            plotCell.Alpha = criticalAlpha(plotCell, 'one-region')+20;

            indicesCell = find(labelledImage == numCell);
            [possibleIndicesX, possibleIndicesY, possibleIndicesZ] = ind2sub(size(labelledImage), indicesCell);
            pixelsIn = inShape(plotCell, possibleIndicesX, possibleIndicesY, possibleIndicesZ);
            labelledImageToCalculateSR(indicesCell(pixelsIn)) = numCell;
        catch
            %try/catch for non-valid cells with strange shape (avoid flatten)
            print(['Cell ' num2str(numCell) ' without flatten'])
            labelledImageToCalculateSR(labelledImage==numCell)=numCell;
        end
    end 
    labelledImageRealSR = double(imclose(labelledImageToCalculateSR > 0,strel('sphere',2)));
    labelledImageRealSR = labelledImageRealSR .* labelledImage;
    labelledImageRealSR(lumenImage>0) = 0;
end

