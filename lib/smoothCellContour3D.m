function [labelledImage] = smoothCellContour3D(labelledImage, idCell, zToSmooth, lumenImage)
%SMOOTHCELLCONTOUR3D Summary of this function goes here
%   Detailed explanation goes here

    [x,y,z] = ind2sub(size(labelledImage),find(labelledImage == idCell));
    
    cellShape = alphaShape(x,y,z);
    newAlpha = criticalAlpha(cellShape, 'one-region');
    cellShape.Alpha = newAlpha;
    %zToSmooth = unique(z);
    
    if getappdata(0, 'canModifyOutsideGland') == 1
        [x,y,z] = ind2sub(size(labelledImage),find(labelledImage ~= idCell));
    else
        [x,y,z] = ind2sub(size(labelledImage),find(labelledImage ~= idCell & labelledImage ~= 0 & lumenImage == 0));
    end
    pixelsToCheck = ismember(z, zToSmooth);

    xToCheck = x(pixelsToCheck);
    yToCheck = y(pixelsToCheck);
    zToCheck = z(pixelsToCheck);
    
    inPixels = inShape(cellShape, xToCheck, yToCheck, zToCheck);
    
    for numInPixels = find(inPixels)'
        labelledImage(xToCheck(numInPixels), yToCheck(numInPixels), zToCheck(numInPixels)) = idCell;
    end
end

