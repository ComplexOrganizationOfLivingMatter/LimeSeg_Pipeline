function samiraTable = connectVerticesOf2D(cylindre2DImage, neighbours2D, vertices2D, centroids, validCellsFinal, borderCells, surfaceRatio, outputDir, nameOfSimulation, deployedImg3x, img3d)
%CONNECTVERTICESOF2D Summary of this function goes here
%   Detailed explanation goes here
    
    maxDistance = 50;
    [~,W]=size(cylindre2DImage);
    W = W/3;
    
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
    
    uniqTriplets = unique(sort(vertcat(cellNeighbours{:}),2),'rows');
    [neighsAccumRea] = getNeighboursFromVertices(uniqTriplets);
    uniqQuartets = [];
    
    for numCell = 1:max(cylindre2DImage(:))
        allNeighsAtCell = cellNeighbours{numCell};
        allVerticesAtCell = cellVertices{numCell};
        neighboursPerSide = {};
        verticesPerSide = {};
        if ismember(numCell, borderCells)
            verticesPerSide{1} = allVerticesAtCell(1:size(allVerticesAtCell, 1)/2, :);
            neighboursPerSide{1} = allNeighsAtCell(1:size(allVerticesAtCell, 1)/2, :);
            
            verticesPerSide{2} = allVerticesAtCell(size(allVerticesAtCell, 1)/2+1:size(allVerticesAtCell, 1), :);
            neighboursPerSide{2} = allNeighsAtCell(size(allVerticesAtCell, 1)/2+1:size(allVerticesAtCell, 1), :);
        else
            verticesPerSide{1} = allVerticesAtCell;
            neighboursPerSide{1} = allNeighsAtCell;
        end
        
        for numParts = 1:(1 + ismember(numCell, borderCells))
            tripletsOfNeighs = neighboursPerSide{numParts};
            %get k4
            quartetsOfNeighs = [];
            for nTriplets = 1 : size(tripletsOfNeighs,1)
                neighsCellTriplet = arrayfun(@(x)  neighsAccumRea{x},tripletsOfNeighs(nTriplets,:),'UniformOutput',false);
                intersectionTriplet= intersect(neighsCellTriplet{1},intersect(neighsCellTriplet{2},neighsCellTriplet{3}));
                if ~isempty(intersectionTriplet)
                    for nQuartets = 1:length(intersectionTriplet)
                        quartetsOfNeighs(end+1,:) = [tripletsOfNeighs(nTriplets,:),intersectionTriplet(nQuartets)];
                    end
                end
            end
            uniqQuartets = [uniqQuartets; unique(sort(quartetsOfNeighs,2),'rows')];
        end
        
        cellVertices{numCell} = vertcat(verticesPerSide{:});
    end
    
    uniqQuartets = unique(sort(uniqQuartets,2),'rows');
    for numQuartet = 1:size(uniqQuartets, 1)
        aQuartet = uniqQuartets(numQuartet, :);
        
        neighboursPerSide{1} = [];
        neighboursPerSide{2} = [];
        verticesPerSide{1} = [];
        verticesPerSide{2} = [];
        
        verticesToChange{1} = [];
        verticesToChange{2} = [];
        
        %% Get vertices of quartet
        sizesOfCells = zeros(4, 2);
        numCell = 1;
        for numCellsToChange = aQuartet
            allNeighsAtCell = cellNeighbours{numCellsToChange};
            allVerticesAtCell = cellVertices{numCellsToChange};
            indRightBorder = (allVerticesAtCell(:,2) - W) > W/2;
            
            verticesPerSide{1} = [verticesPerSide{1}; allVerticesAtCell(indRightBorder, :)];
            verticesPerSide{2} = [verticesPerSide{2}; allVerticesAtCell(~indRightBorder, :)];
            
            neighboursPerSide{1} = [neighboursPerSide{1}; allNeighsAtCell(indRightBorder, :)];
            neighboursPerSide{2} = [neighboursPerSide{2}; allNeighsAtCell(~indRightBorder, :)];
            
            sizesOfCells(numCell, 1) = sum(indRightBorder);
            sizesOfCells(numCell, 2) = sum(indRightBorder==0);
            
            numCell = numCell + 1;
        end
        
        %% Get new vertices
        verticesPerSide_old = verticesPerSide;
        for numParts = 1:(1+ (isempty(verticesPerSide{2}) == 0))
            tripletsOfNeighs = neighboursPerSide{numParts};
            verticesToChange = all(ismember(tripletsOfNeighs, aQuartet), 2);
            newValueCentroid = round(mean(verticesPerSide{numParts}(verticesToChange, :)));
            verticesPerSide{numParts}(verticesToChange, :) = repmat(newValueCentroid, sum(verticesToChange), 1);
        end
        
        %% Replace old vertices
        for numCellsToChange = aQuartet
            allNeighsAtCell = cellNeighbours{numCellsToChange};
            allVerticesAtCell = cellVertices{numCellsToChange};
            indRightBorder = (allVerticesAtCell(:,2) - W) > W/2;
            verticesPerSide_old
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

