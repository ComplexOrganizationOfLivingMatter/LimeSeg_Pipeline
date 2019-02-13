function samiraTable = connectVerticesOf2D(cylindre2DImage, neighbours2D, vertices2D, centroids, midSectionNewLabels, wholeImage, validCellsFinal, cellNumNeighbours, borderCells, surfaceRatio, outputDir, nameOfSimulation)
%CONNECTVERTICESOF2D Summary of this function goes here
%   Detailed explanation goes here
    
    midSectionImage_closed = imclose(cylindre2DImage>0, strel('disk', 2));
    sizeRoll = sum(midSectionImage_closed, 2);
    for numCell = 1:max(cylindre2DImage(:))
        actualVertices = any(ismember(neighbours2D, numCell), 2);
        actualCellVertices = vertices2D(actualVertices, :);
        if ismember(numCell, borderCells)
            % Obtaining the vertices of both sides if it is a border cell
            changeOfSide = actualCellVertices(:, 2)-(size(cylindre2DImage, 2)/2);
            newVertices = [actualCellVertices(:, 1), actualCellVertices(:, 2) - (sign(changeOfSide) .* sizeRoll(actualCellVertices(:, 1)))];
%             figure; imshow(cylindre2DImage)
%             hold on;
%             actualVerticesRegion = actualCellVertices;
%             for numVertex = 1:size(actualVerticesRegion, 1)
%                 plot(actualVerticesRegion(numVertex, 2), actualVerticesRegion(numVertex, 1), 'rx')
%             end
            
            actualCellVertices = [actualCellVertices; newVertices];
        end
        cellVertices{numCell} = actualCellVertices;
    end
    
    noValidCells = ones(max(cylindre2DImage(:)), 1);
    noValidCells(validCellsFinal) = 0;
    
    cellVerticesNoValid = cellVertices;
    cellVerticesNoValid(noValidCells==0) = {[]};
    
    cellVerticesValid = cellVertices;
    cellVerticesValid(noValidCells==1) = {[]};
    
    ySize = size(cylindre2DImage, 2);
    cellInfoWithVertices = groupingVerticesPerCellSurface(cylindre2DImage(:, (ySize/3):(2*ySize/3)), cellVerticesValid, cellVerticesNoValid, [], 1, borderCells);
    cellInfoWithVertices(cellfun(@isempty, cellInfoWithVertices(:, 6)), :) = [];
    cellInfoWithVertices(cellfun(@(x) ismember(x, noValidCells), cellInfoWithVertices(:, 3)), :) = [];
    
    [samiraTable, cellsVoronoi] = tableWithSamiraFormat(cellInfoWithVertices, centroids, [], surfaceRatio, outputDir, nameOfSimulation);
end

