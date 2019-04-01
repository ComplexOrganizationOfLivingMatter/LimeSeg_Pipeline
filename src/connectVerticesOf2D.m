function samiraTable = connectVerticesOf2D(cylindre2DImage, neighbours2D, vertices2D, centroids, validCellsFinal, borderCells, surfaceRatio, outputDir, nameOfSimulation, deployedImg3x, img3d)
%CONNECTVERTICESOF2D Summary of this function goes here
%   Detailed explanation goes here
    
    maxDistance = 50;

    midSectionImage_closed = imclose(cylindre2DImage>0, strel('disk', 2));
    sizeRoll = sum(midSectionImage_closed, 2);
    cellVertices = cell(1, max(cylindre2DImage(:)));
    cellNeighbours = cell(1, max(cylindre2DImage(:)));
    
    for numCell = 1:max(cylindre2DImage(:))
        actualVertices = any(ismember(neighbours2D, numCell), 2);
        actualCellVertices = vertices2D(actualVertices, :);
        actualCellNeighbours = neighbours2D(actualVertices, :);
        if ismember(numCell, borderCells)
%             % Obtaining the vertices of both sides if it is a border cell
%             changeOfSide = actualCellVertices(:, 2)-(size(cylindre2DImage, 2)/2);
%             newVertices = [actualCellVertices(:, 1), actualCellVertices(:, 2) - (sign(changeOfSide) .* sizeRoll(actualCellVertices(:, 1)))];
% %             figure; imshow(cylindre2DImage)
% %             hold on;
% %             actualVerticesRegion = actualCellVertices;
% %             for numVertex = 1:size(actualVerticesRegion, 1)
% %                 plot(actualVerticesRegion(numVertex, 2), actualVerticesRegion(numVertex, 1), 'rx')
% %             end
%             
%             actualCellVertices = [actualCellVertices; newVertices];
            [actualCellNeighbours, actualCellVertices] = obtainVerticesOfBorderCells(cylindre2DImage, deployedImg3x, img3d, numCell);
        end
        
        for numCellNeighbour = unique(actualCellNeighbours(:))'
            verticesOfActualNeighbour = cellVertices{numCellNeighbour};
            if isempty(verticesOfActualNeighbour) == 0
                %% Unify vertices that are the same triangulation
                verticesConnectingActualNeighbour = cellNeighbours{numCellNeighbour};
                [badVertices, goodVertices] = find(pdist2(actualCellVertices, verticesOfActualNeighbour) < maxDistance);
                
                sharingNeighboursVertex = all(arrayfun(@isequal, actualCellNeighbours(badVertices, :), verticesConnectingActualNeighbour(goodVertices, :)), 2);
                actualCellVertices(badVertices(sharingNeighboursVertex), :) = verticesOfActualNeighbour(goodVertices(sharingNeighboursVertex), :);
            end
        end
        
        cellVertices{numCell} = actualCellVertices;
        cellNeighbours{numCell} = actualCellNeighbours;
       
    end
    
    for numCell = 1:max(cylindre2DImage(:))
        allNeighsAtCell = cellNeighbours{numCell};
        allVerticesAtCell = cellVertices{numCell};
        neighboursPerSide = {};
        verticesPerSide = {};
        if ismember(numCell, borderCells)
            neighboursPerSide{1} = allVerticesAtCell(1:size(allVerticesAtCell, 1)/2, :);
            verticesPerSide{1} = allNeighsAtCell(1:size(allVerticesAtCell, 1)/2, :);
            
            neighboursPerSide{2} = allVerticesAtCell(size(allVerticesAtCell, 1)/2+1:size(allVerticesAtCell, 1), :);
            verticesPerSide{2} = allNeighsAtCell(size(allVerticesAtCell, 1)/2+1:size(allVerticesAtCell, 1), :);
        else
            neighboursPerSide{1} = allVerticesAtCell;
            verticesPerSide{1} = allNeighsAtCell;
        end
        
        for numParts = 1:(1 + ismember(numCell, borderCells))
            coordXY = neighboursPerSide{numParts};
            uniqueXY = unique(coordXY,'rows');
            numRepet = arrayfun(@(x,y) sum(ismember(coordXY,[x y],'rows')),uniqueXY(:,1),uniqueXY(:,2));
            crossesVert = ismember(coordXY,uniqueXY(numRepet>3,:),'rows') & ~isnan(coordXY(:,1));
            if any(crossesVert)>0
                disp('');
            end
        end
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
    cellInfoWithVertices(cellfun(@(x) ismember(x, find(noValidCells)), cellInfoWithVertices(:, 3)), :) = [];
    
    [samiraTable, cellsVoronoi] = tableWithSamiraFormat(cellInfoWithVertices, centroids, [], surfaceRatio, outputDir, nameOfSimulation);
end

