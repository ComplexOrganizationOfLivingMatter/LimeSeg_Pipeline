function [labelledImageWithGoodIds] = matchingIDs(originalImageToRename, imageGoodIds)
%MATCHINGIDS Summary of this function goes here
%   Detailed explanation goes here

    allCells = unique(originalImageToRename(:))';
    allCells(allCells==0) = [];
    centroidsGoodsIds = regionprops(originalImageToRename, 'Centroid');
    centroidsGoodsIds = round(vertcat(centroidsGoodsIds.Centroid));
    
    labelledImageWithGoodIds = zeros(size(originalImageToRename));
    for numCell = allCells
        possibleCellIds = unique(imageGoodIds(originalImageToRename == numCell));
        possibleCellIds(possibleCellIds == 0) = [];
        
        modeValue = mode(imageGoodIds(originalImageToRename == numCell));
        centroidValue = imageGoodIds(centroidsGoodsIds(numCell, 2), centroidsGoodsIds(numCell, 1), centroidsGoodsIds(numCell, 3));
        if modeValue ~= centroidValue
            warning('Incoherence between centroid and mode'); 
        end
        
        labelledImageWithGoodIds(originalImageToRename == numCell) = modeValue;
    end
    
end

