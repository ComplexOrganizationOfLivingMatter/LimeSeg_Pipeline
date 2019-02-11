function connectVerticesOf2D(midSectionImage, neighbours2D, vertices2D, centroids, midSectionNewLabels, wholeImage, validCellsFinal, cellNumNeighbours, borderCells)
%CONNECTVERTICESOF2D Summary of this function goes here
%   Detailed explanation goes here
    
    for numVertex = 1:size(vertices2D, 1)
        actualEdges = cellfun(@(x) length(intersect(x, neighbours2D(numVertex, :))), mat2cell(neighbours2D, ones(1, size(neighbours2D, 1)), 3));
        edgesOfVertex{numVertex, 1} = vertices2D(actualEdges>1, :);
        edgesOfVertex{numVertex, 2} = neighbours2D(actualEdges>1, :);
        edgesOfVertex{numVertex, 3} = actualEdges(actualEdges>1);
    end

    for numCell = 1:max(midSectionImage(:))
        actualVertices = any(ismember(neighbours2D, numCell), 2);
        cellVertices{numCell} = vertices2D(actualVertices, :);
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
    
    for numCentroid = 1:size(centroids, 1)
        numCell = midSectionImage(midSectionNewLabels == numCentroid);
        numCell = numCell(1);
%         if ismember(numCell, validCellsFinal) == 0
%             continue
%         end
        actualVertices = any(ismember(neighbours2D, numCell), 2);
        
        actualImg = wholeImage == numCell;
        centroids3x = regionprops(actualImg, 'Centroid');
        centroids3x = vertcat(centroids3x.Centroid);
        
        if size(centroids3x, 1) > 3
           disp('CAREEEEE') 
        end
        
        if ismember(numCell, borderCells)
            
            
            vertices2D
        else
            cellVertices = vertices2D(actualVertices, :);
            
        end
        
%         allActualVertices = [vertices2D(actualVertices, :); vertices2D_Left(actualVertices, :); vertices2D_Right(actualVertices, :)];
%         actualNeighbours2D = [neighbours2D(actualVertices, :); neighbours2D(actualVertices, :); neighbours2D(actualVertices, :)];
%         [yValidRegion, xValidRegion] = find(imdilate(actualImg, strel('disk', 10)));
%         
%         %Get only the closest vertices to the captured region
%         actualNeighbours2D(ismember(allActualVertices, [xValidRegion, yValidRegion], 'rows') == 0, :) = [];
%         allActualVertices(ismember(allActualVertices, [xValidRegion, yValidRegion], 'rows') == 0, :) = [];
%         
%         [~, closestIndices] = pdist2(centroids3x, allActualVertices, 'euclidean', 'Smallest', 1);
% 
%         [~, midCentroid] = pdist2(centroids3x, centroids(numCentroid, :), 'euclidean', 'Smallest', 1);
%         
%         sum(actualVertices)
%         cellNumNeighbours(numCell)
        
%         if sum(actualVertices) ~= cellNumNeighbours(numCell)
%             for numRegion = 1:size(centroids3x, 1)
%                 figure; imshow(actualImg)
%                 hold on;
%                 actualVerticesRegion = allActualVertices(closestIndices == numRegion, :);
%                 for numVertex = 1:size(actualVerticesRegion, 1)
%                     plot(actualVerticesRegion(numVertex, 1), actualVerticesRegion(numVertex, 2), 'x')
%                 end
%             end
%         end
        selectedNeighboursOfVertices = actualNeighbours2D(closestIndices == midCentroid, :);
        if sum(closestIndices == midCentroid) < 2
            continue;
        end
        if size(unique(selectedNeighboursOfVertices, 'rows'), 1) ~= size(allActualVertices(closestIndices == midCentroid, :), 1)
            %Wrong vertices
            
            % Remove pixels that share neighbours 
            % The ideal case would be the one in which the same vertex is
            % shared with other cells
            numCell
            [newVertOrder] = boundaryOfCell(allActualVertices(closestIndices == midCentroid, :), centroids(numCentroid, :));
        else
            [newVertOrder] = boundaryOfCell(allActualVertices(closestIndices == midCentroid, :), centroids(numCentroid, :));
        end
        
        %
        [newOrderX, newOrderY] = poly2cw(newVertOrder((1:end-1), 1), newVertOrder((1:end-1), 2));
        verticesRadius = [];
        
             
        hold on
        plot(newVertOrder(:, 1), newVertOrder(:, 2))
        plot(allActualVertices(closestIndices == midCentroid, 1), allActualVertices(closestIndices == midCentroid, 2), 'r+');
    end
    
    
end

