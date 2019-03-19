function [neighbours] = getNeighboursFromVertices(verticesNeighbours)
%GETNEIGHBOURSFROMVERTICES Summary of this function goes here
%   Detailed explanation goes here
            neighbours = arrayfun(@(x) unique(verticesNeighbours(any(ismember(verticesNeighbours, x), 2), :)), 1:max(verticesNeighbours(:)), 'UniformOutput', false);
            neighbours = cellfun(@(x, y) x(x ~= y), neighbours, num2cell(1:length(neighbours)), 'UniformOutput', false);
            neighbours = cellfun(@(x) x(:), neighbours, 'UniformOutput', false);
end

