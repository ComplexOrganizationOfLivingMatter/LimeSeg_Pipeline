function [verticesNeighs2D, vertices2D] = obtainVerticesOfBorderCells(deployedImg, deployedImg3x, verticesNeighs2D, vertices2D, borderCells)
%OBTAINVERTICESOFBORDERCELLS Summary of this function goes here
%   Detailed explanation goes here

    newVertices{1} = [];
    newVertices{2} = [];
    newVerticesNeighs2D{1} = [];
    newVerticesNeighs2D{2} = [];
    for numBorderCell = borderCells'
        labelledActualCells3x = bwlabel(deployedImg3x==numBorderCell);
        newLabelBorderCells = unique(labelledActualCells3x .* ismember(deployedImg, numBorderCell));
        newLabelBorderCells(newLabelBorderCells==0) = [];
        
        cellsToCalculateNeighs = ismember(labelledActualCells3x, newLabelBorderCells);
        regionToCalculateNeighs = imdilate(cellsToCalculateNeighs, strel('disk', 3));
        
        leftAndRightSide = bwlabel(regionToCalculateNeighs);
        
        for numSide = 1:2
            img_sideOfBorderCell = (leftAndRightSide == numSide) .* deployedImg3x;
            neighbours = calculateNeighbours(img_sideOfBorderCell);
            [newVerticesActual] = getVertices(img_sideOfBorderCell, neighbours);
            newVertices{numSide} = [newVertices{numSide}; vertcat(newVerticesActual.verticesPerCell{:})];
            newVerticesNeighs2D{numSide} = [newVerticesNeighs2D{numSide}; newVerticesActual.verticesConnectCells];
        end
    end
    [newVerticesNeighs2D{1}, IA,IC] = unique(newVerticesNeighs2D{1}, 'rows');
    verticesToRemove = any(ismember(verticesNeighs2D, borderCells), 2);
    verticesNeighs2D(verticesToRemove, :) = [];
    vertices2D(verticesToRemove, :) = [];
    
    vertices2D = [vertices2D; newVertices];
    verticesNeighs2D = [verticesNeighs2D; newVerticesNeighs2D];
end

