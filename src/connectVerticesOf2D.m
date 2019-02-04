function connectVerticesOf2D(midSectionImage, neighbours2D, vertices2D, vertices2D_Left, vertices2D_Right, centroids, midSectionNewLabels, wholeImage, validCellsFinal, cellNumNeighbours)
%CONNECTVERTICESOF2D Summary of this function goes here
%   Detailed explanation goes here
    
    for numCentroid = 1:size(centroids, 1)
        numCell = midSectionImage(midSectionNewLabels == numCentroid);
        numCell = numCell(1);
        if ismember(numCell, validCellsFinal) == 0
            continue
        end
        actualVertices = any(ismember(neighbours2D, numCell), 2);
        
        actualImg = wholeImage == numCell;
        centroids3x = regionprops(actualImg, 'Centroid');
        centroids3x = vertcat(centroids3x.Centroid);
        
        if size(centroids3x, 1) > 3
           disp('CAREEEEE') 
        end
        
        allActualVertices = [vertices2D(actualVertices, :); vertices2D_Left(actualVertices, :); vertices2D_Right(actualVertices, :)];
        
        [yValidRegion, xValidRegion] = find(imdilate(actualImg, strel('disk', 10)));
        
        %Get only the closest vertices to the captured region
        allActualVertices(ismember(allActualVertices, [xValidRegion, yValidRegion], 'rows') == 0, :) = [];
        
        [~, closestIndices] = pdist2(centroids3x, allActualVertices, 'euclidean', 'Smallest', 1);

        [~, midCentroid] = pdist2(centroids3x, centroids(numCentroid, :), 'euclidean', 'Smallest', 1);
        
        sum(actualVertices)
        cellNumNeighbours(numCell)
        
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
        
        [newVertOrder] = boundaryOfCell(allActualVertices(closestIndices == midCentroid, :), centroids(numCentroid, :));
        [newOrderX, newOrderY] = poly2cw(newVertOrder((1:end-1), 1), newVertOrder((1:end-1), 2));
        verticesRadius = [];
        
             
        hold on
        plot(newVertOrder(:, 1), newVertOrder(:, 2))
        plot(allActualVertices(closestIndices == midCentroid, 1), allActualVertices(closestIndices == midCentroid, 2), 'r+');
    end
    
    
end

