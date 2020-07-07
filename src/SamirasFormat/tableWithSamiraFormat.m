function [samiraTableVoronoi, cellsVoronoi] = tableWithSamiraFormat(cellWithVertices,cellCentroids, missingVertices, nSurfR, outputDir, nameOfSimulation)
%TABLEWITHSAMIRAFORMAT Summary of this function goes here
%   Detailed explanation goes here
%
    samiraTableVoronoi = {};
    missingVerticesCoord = [];
    for numCell = 1:size(cellWithVertices, 1)
        try
            verticesOfCellInit = cellWithVertices{numCell, end};

            numberOfVertices = (size(verticesOfCellInit, 2)/2);
            verticesOfCell = [];
            verticesOfCell(1:numberOfVertices, 1) = verticesOfCellInit(1:2:end);
            verticesOfCell(1:numberOfVertices, 2) = verticesOfCellInit(2:2:end);

            %Replace the missing cells
            for numPair = 1:size(missingVertices, 1)
                verticesToReplace = ismember(verticesOfCell, missingVertices(numPair, 1:2), 'rows');
                replaceWith = missingVertices(numPair, 3:4);

                if any(verticesToReplace) > 0
                    verticesOfCell(verticesToReplace, :) = replaceWith;
                end
            end

            verticesOfCell = double(unique(verticesOfCell, 'rows'));

            [newVertOrder] = boundaryOfCell(verticesOfCell, cellCentroids(cellWithVertices{numCell,3}, :));

            % Should be connected clockwise
            % I.e. from bigger numbers to smaller ones
            % Or the second vertex should in the left hand of the first
            [newOrderX, newOrderY] = poly2cw(newVertOrder((1:end-1), 1), newVertOrder((1:end-1), 2));
            verticesRadius = [];

    % %         figure;       
    %         hold on
    %         plot(newVertOrder(:, 1), newVertOrder(:, 2))
    %         hold on
    %         plot(verticesOfCell(:, 1), verticesOfCell(:, 2), 'r+');

            for numVertex = 1:length(newOrderX)
                verticesRadius(end+1) = newOrderX(numVertex);
                verticesRadius(end+1) = newOrderY(numVertex);
            end

            cellsVoronoi = [nSurfR, cellWithVertices{numCell, 3:5}, {verticesRadius}];
            samiraTableVoronoi = [samiraTableVoronoi; cellsVoronoi];
        catch
           ['problems with the vertices of ' num2str(numCell) ' cell'] 
        end
    end
    
    %Plot
    plotVerticesPerSurfaceRatio(samiraTableVoronoi((end-size(samiraTableVoronoi,1)+1):end,:),missingVerticesCoord, outputDir, nameOfSimulation, nSurfR)
end

