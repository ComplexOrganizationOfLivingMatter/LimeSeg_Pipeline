function [apicalLayer,basalLayer,labelledImageRealSR, lumenImageRealSR] = flattenMutantGland(apicalLayer, basalLayer, labelledImage)
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
    
    tipValue = 4;
    labelledImageRealSR = double(imclose(labelledImageToCalculateSR > 0,strel('sphere',2)));
    filledGland = imfill(labelledImageRealSR);
    labelledImageRealSR = labelledImageRealSR .* labelledImage;
    lumenImageRealSR = bwlabeln(labelledImageRealSR==0 & (filledGland>0));
    volumeLumenImage = regionprops3(lumenImageRealSR,'Volume');
    [~,indexLumenImage] = max(volumeLumenImage.Volume);
    [apicalLayer] = getApicalFrom3DImage(lumenImageRealSR==indexLumenImage, labelledImageRealSR);
    [basalLayer] = getBasalFrom3DImage(labelledImage, lumenImageRealSR==indexLumenImage, tipValue, filledGland==0);
    
    lumenImageRealSR = lumenImageRealSR==indexLumenImage;
end

