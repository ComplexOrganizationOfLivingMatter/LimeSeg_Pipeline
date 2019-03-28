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
        regionToCalculateNeighs = double(imdilate(cellsToCalculateNeighs, strel('disk', 1)));
        
        [~, y] = find(regionToCalculateNeighs);
        indicesMidImage = find(regionToCalculateNeighs);
        
        regionToCalculateNeighs(indicesMidImage(y < (size(deployedImg3x, 2) / 2))) = 1;
        regionToCalculateNeighs(indicesMidImage(y >= (size(deployedImg3x, 2) / 2))) = 2;
        
        for numSide = 1:2
            img_sideOfBorderCell = (regionToCalculateNeighs == numSide) .* deployedImg3x;
            neighbours = calculateNeighbours(img_sideOfBorderCell);
            [newVerticesActual] = getVertices(img_sideOfBorderCell, neighbours);
            newVertices{numSide} = [newVertices{numSide}; vertcat(newVerticesActual.verticesPerCell{:})];
            newVerticesNeighs2D{numSide} = [newVerticesNeighs2D{numSide}; newVerticesActual.verticesConnectCells];
        end
    end
    
    [newVerticesNeighs2D{1}, indicesUnique_1] = unique(newVerticesNeighs2D{1}, 'rows');
    newVertices{1} = newVertices{1}(indicesUnique_1, :);
    [newVerticesNeighs2D{2}, indicesUnique_2] = unique(newVerticesNeighs2D{2}, 'rows');
    newVertices{2} = newVertices{2}(indicesUnique_2, :);
    
    verticesToRemove = any(ismember(verticesNeighs2D, borderCells), 2);
    verticesNeighs2D(verticesToRemove, :) = [];
    vertices2D(verticesToRemove, :) = [];
    
    vertices2D = [vertices2D; newVertices{1}; newVertices{2}];
    verticesNeighs2D = [verticesNeighs2D; newVerticesNeighs2D{1}; newVerticesNeighs2D{2}];
end

