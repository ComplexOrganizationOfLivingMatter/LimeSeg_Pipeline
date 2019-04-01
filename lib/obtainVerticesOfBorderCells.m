function [verticesNeighs2D, vertices2D] = obtainVerticesOfBorderCells(deployedImg, deployedImg3x, img3d, numBorderCell)
%OBTAINVERTICESOFBORDERCELLS Summary of this function goes here
%   Detailed explanation goes here

    newVertices{1} = [];
    newVertices{2} = [];
    newVerticesNeighs2D{1} = [];
    newVerticesNeighs2D{2} = [];
    labelledActualCells3x = bwlabel(deployedImg3x==numBorderCell);
    newLabelBorderCells = unique(labelledActualCells3x .* ismember(deployedImg, numBorderCell));
    newLabelBorderCells(newLabelBorderCells==0) = [];

    cellsToCalculateNeighs = ismember(labelledActualCells3x, newLabelBorderCells);
    regionToCalculateNeighs = double(imdilate(cellsToCalculateNeighs, strel('disk', 4)));

    [~, y] = find(regionToCalculateNeighs);
    indicesMidImage = find(regionToCalculateNeighs);

    regionToCalculateNeighs(indicesMidImage(y < (size(deployedImg3x, 2) / 2))) = 1;
    regionToCalculateNeighs(indicesMidImage(y >= (size(deployedImg3x, 2) / 2))) = 2;

    for numSide = 1:2
        img_sideOfBorderCell = (regionToCalculateNeighs == numSide) .* deployedImg3x;
        neighbours = calculateNeighbours(img_sideOfBorderCell);
        %neighbours = checkPairPointCloudDistanceCurateNeighbours(img3d, neighbours);
        [newVerticesActual] = getVertices(img_sideOfBorderCell, neighbours);
        emptyCells = any(ismember(newVerticesActual.verticesConnectCells, numBorderCell), 2)==0;
        %emptyCells = cellfun(@isempty, newVerticesActual.verticesPerCell);
        newVertices{numSide} = [newVertices{numSide}; vertcat(newVerticesActual.verticesPerCell{emptyCells==0})];
        newVerticesNeighs2D{numSide} = [newVerticesNeighs2D{numSide}; newVerticesActual.verticesConnectCells(emptyCells==0, :)];
%         newVertices{numSide} = [newVertices{numSide}; vertcat(newVerticesActual.verticesPerCell{:})];
%         newVerticesNeighs2D{numSide} = [newVerticesNeighs2D{numSide}; newVerticesActual.verticesConnectCells];
    end

    [newVerticesNeighs2D{1}, indicesUnique_1] = unique(newVerticesNeighs2D{1}, 'rows');
    newVertices{1} = newVertices{1}(indicesUnique_1, :);
    [newVerticesNeighs2D{2}, indicesUnique_2] = unique(newVerticesNeighs2D{2}, 'rows');
    newVertices{2} = newVertices{2}(indicesUnique_2, :);

    vertices2D = [newVertices{1}; newVertices{2}];
    verticesNeighs2D = [newVerticesNeighs2D{1}; newVerticesNeighs2D{2}];
end

