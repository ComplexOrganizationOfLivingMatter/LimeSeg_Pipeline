function connectVerticesOf2D(midSectionImage, neighbours2D, vertices2D, centroids, midSectionNewLabels, wholeImage, validCellsFinal, cellNumNeighbours, borderCells)
%CONNECTVERTICESOF2D Summary of this function goes here
%   Detailed explanation goes here

    for numCell = 1:max(midSectionImage(:))
        actualVertices = any(ismember(neighbours2D, numCell), 2);
        actualCellVertices = vertices2D(actualVertices, 2:-1:1);
        if ismember(numCell, borderCells)
            % Obtaining the vertices of both sides if it is a border cell
            changeOfSide = actualCellVertices(:, 1)-(size(midSectionImage, 2)/2);
            newVertices = [actualCellVertices(:, 1) - changeOfSide*2, actualCellVertices(:, 2)];
            figure; imshow(midSectionImage)
            hold on;
            actualVerticesRegion = newVertices;
            for numVertex = 1:size(actualVerticesRegion, 1)
                plot(actualVerticesRegion(numVertex, 1), actualVerticesRegion(numVertex, 2), 'x')
            end
            
            actualCellVertices = [actualCellVertices; newVertices];
        end
        cellVertices{numCell} = actualCellVertices;
    end
    
    noValidCells = ones(max(midSectionImage(:)), 1);
    noValidCells(validCellsFinal) = 0;
    
    cellVerticesNoValid = cellVertices;
    cellVerticesNoValid(noValidCells==0) = {[]};
    
    cellVerticesValid = cellVertices;
    cellVerticesValid(noValidCells==1) = {[]};
    
    ySize = size(midSectionImage, 2);
    cellInfoWithVertices = groupingVerticesPerCellSurface(midSectionImage(:, (ySize/3):(2*ySize/3)), cellVerticesValid, cellVerticesNoValid, [], 1, borderCells);
    cellInfoWithVertices(cellfun(@isempty, cellInfoWithVertices(:, 6)), :) = [];
    cellInfoWithVertices(cellfun(@(x) ismember(x, noValidCells), cellInfoWithVertices(:, 3)), :) = [];
    
    [samiraTable, cellsVoronoi] = tableWithSamiraFormat(cellInfoWithVertices, centroids, [], 1, '.', 'Voronoi_'); 
end

