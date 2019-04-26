function [verticesNeighs2D, vertices2D, borderCell] = obtainVerticesOfBorderCells(deployedImg, deployedImg3x, img3d, numBorderCell)
%OBTAINVERTICESOFBORDERCELLS Summary of this function goes here
%   Detailed explanation goes here

    newVertices{1} = [];
    newVertices{2} = [];
    newVerticesNeighs2D{1} = [];
    newVerticesNeighs2D{2} = [];
    
    vertices2D = [];
    verticesNeighs2D = [];
    
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
        neighbourOfBorderCell = any(ismember(newVerticesActual.verticesConnectCells, numBorderCell), 2)==0;
        emptyCells = cellfun(@isempty, newVerticesActual.verticesPerCell);
        
        if any(emptyCells)
            numRadius = 3;
            for numEmptyCell = find(emptyCells)'
                newVertexX = [];
                while isempty(newVertexX)
                    dilatedImage_1 = imdilate(ismember(img_sideOfBorderCell, newVerticesActual.verticesConnectCells(numEmptyCell, 1)), strel('disk', numRadius));
                    dilatedImage_2 = imdilate(ismember(img_sideOfBorderCell, newVerticesActual.verticesConnectCells(numEmptyCell, 2)), strel('disk', numRadius));
                    dilatedImage_3 = imdilate(ismember(img_sideOfBorderCell, newVerticesActual.verticesConnectCells(numEmptyCell, 3)), strel('disk', numRadius));
                    
                    [newVertexX, newVertexY] = find(dilatedImage_1 & dilatedImage_2 & dilatedImage_3);
                    numRadius = numRadius + 1;
                end
                newVertex = round(mean([newVertexX, newVertexY], 1));
                newVerticesActual.verticesPerCell{numEmptyCell} = newVertex;
            end
        end

        newVertices{numSide} = [newVertices{numSide}; vertcat(newVerticesActual.verticesPerCell{neighbourOfBorderCell==0})];
        newVerticesNeighs2D{numSide} = [newVerticesNeighs2D{numSide}; newVerticesActual.verticesConnectCells(neighbourOfBorderCell==0, :)];
%         newVertices{numSide} = [newVertices{numSide}; vertcat(newVerticesActual.verticesPerCell{:})];
%         newVerticesNeighs2D{numSide} = [newVerticesNeighs2D{numSide}; newVerticesActual.verticesConnectCells];
    end
    
    if sum(cellfun(@isempty, newVerticesNeighs2D)) == 0
        borderCell = sum(cellfun(@(x) size(x, 1) == 1, newVerticesNeighs2D)) == 0;
        newVerticesNeighs2D(cellfun(@(x) size(x, 1) == 1, newVerticesNeighs2D)) = {[]};
    else
        borderCell = 0;
    end
    
    [newVerticesNeighs2D{1}, indicesUnique_1] = unique(newVerticesNeighs2D{1}, 'rows');
    newVertices{1} = newVertices{1}(indicesUnique_1, :);
    [newVerticesNeighs2D{2}, indicesUnique_2] = unique(newVerticesNeighs2D{2}, 'rows');
    newVertices{2} = newVertices{2}(indicesUnique_2, :);
    
    vertices2D = [newVertices{1}; newVertices{2}];
    verticesNeighs2D = [newVerticesNeighs2D{1}; newVerticesNeighs2D{2}];
end

