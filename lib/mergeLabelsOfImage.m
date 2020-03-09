function [labelledImage] = mergeLabelsOfImage(labelledImage, cellsToMerge)
%MERGELABELSOFIMAGE Summary of this function goes here
%   Detailed explanation goes here
    for numCellsToMerge = 2:length(cellsToMerge)
        labelledImage(cellsToMerge(numCellsToMergemain)) = cellsToMerge(1);
    end
end

