function [numNeighOfNeighPerSurfacesRealization, numNeighPerSurfaceRealization] = getNumNeighsOfNeighs(neighsSurface, numberOfSurfaceRatios)
%GETNUMNEIGHSOFNEIGHS Summary of this function goes here
%   Detailed explanation goes here
    numNeighPerSurfaceRealization = cellfun(@(x) length(x),neighsSurface);

    numNeighOfNeighPerSurfacesRealization = zeros(size(neighsSurface));
    for nSR = 1:numberOfSurfaceRatios
        numNeighOfNeighPerSurfacesRealization(:,nSR) = cellfun(@(x) sum(vertcat(numNeighPerSurfaceRealization(x,nSR)))/length(x),neighsSurface(:,nSR));
    end
end

