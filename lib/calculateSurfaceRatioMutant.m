function [realSurfaceRatio] = calculateSurfaceRatioMutant(basalLayer, apicalLayer, labelledImage, validCells)
%CALCULATESURFACERATIOMUNTANT Summary of this function goes here
%   We ought to get the surface ratio (RadiusBasal/RadiusApical). In this
%   case we want to get RadiusBasal as a cylindre amid two basal surfaces:
%   the one forming the outer gland (the maximum surface formed by the max
%   cell height of each cell) and an inner basal surface formed by the
%   minimum points of the cells intersecting.
    
    %% Get real basal layer
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
    
    %% Calculate real Surface Ratio
    apicalAreaCells=regionprops(apicalLayer,'Area');
    basalAreaCells=regionprops(basalLayer,'Area');
    apicalAreaCells = apicalAreaCells(validCells);
    basalAreaCells = basalAreaCells(validCells);
    realSurfaceRatio = mean([basalAreaCells.Area]) / mean([apicalAreaCells.Area]);
    
end

