function [newVerticesNeighs2D, newVertices] = obtainVerticesOfBorderCells(deployedImg, deployedImg3x, newVerticesNeighs2D, newVertices2D, borderCells)
%OBTAINVERTICESOFBORDERCELLS Summary of this function goes here
%   Detailed explanation goes here

    for numBorderCell = borderCells
        labelledActualCells3x = bwlabel(deployedImg3x==numBorderCell);
        %calculateNeighbours(
    end
end

