function [areaCellsPerSurfaceRealization, volumePerSurfaceRealization, neighsSurface, neighsAccumSurfaces, percentageScutoids, apicoBasalTransitions, numLostNeighsAccum, numWonNeighsAccum] = getBasicInformationToPlot(infoPerSurfaceRatio, neighboursOfAllSurfaces, numberOfSurfaceRatios)
%GETBASICINFORMATIONTOPLOT Summary of this function goes here
%   Detailed explanation goes here
        %% Calculate variables per surface
        %Initialize all variables
        neighsSurface = cell(numberOfSurfaceRatios,1);
        neighsAccumSurfaces = cell(numberOfSurfaceRatios,1);
        percentageScutoids = cell(numberOfSurfaceRatios, 1);
        apicoBasalTransitions = cell(numberOfSurfaceRatios, 1);
        areaCells = cell(numberOfSurfaceRatios,1);
        volumes = cell(numberOfSurfaceRatios,1);
        
        infoOfCells = infoPerSurfaceRatio{1, 4};
        infoOfCells = infoOfCells{:};
        
        %Apical sides
        neighsSurface{1} = neighboursOfAllSurfaces{1};
        neighsAccumSurfaces{1} = neighboursOfAllSurfaces{1};
        percentageScutoids{1} = cellfun(@(x, y) ~isequal(x,y), neighsSurface{1}, neighsAccumSurfaces{1});
        numLostNeighsAccum{1} = cell(size(neighboursOfAllSurfaces{1}));
        numWonNeighsAccum{1} = cell(size(neighboursOfAllSurfaces{1}));
        areaCells(1) = {infoOfCells.Basal_area};
        volumes(1) = {infoOfCells.Volume};
        
        for nSR = 2:numberOfSurfaceRatios
            neighsSurface{nSR} = neighboursOfAllSurfaces{nSR};
            neighsAccumSurfaces{nSR} = cellfun(@(x,y) unique([x;y]),neighsAccumSurfaces{nSR-1},neighsSurface{nSR},'UniformOutput',false);
            percentageScutoids{nSR} = cellfun(@(x, y) ~isempty(setxor(x,y)), neighsSurface{1}, neighsSurface{nSR});

            lostNeigh = cellfun(@(x, y) setdiff(x,y), neighsAccumSurfaces{nSR-1}, neighsSurface{nSR}, 'UniformOutput',false);
            wonNeigh = cellfun(@(x, y) setdiff(y, x), neighsAccumSurfaces{nSR-1}, neighsAccumSurfaces{nSR}, 'UniformOutput',false);

            numLostNeighsAccum{nSR} = cellfun(@(x,y) unique([x;y]),lostNeigh,numLostNeighsAccum{nSR-1},'UniformOutput',false);
            numWonNeighsAccum{nSR} = cellfun(@(x,y) unique([x;y]),wonNeigh,numWonNeighsAccum{nSR-1},'UniformOutput',false);

            apicoBasalTransitions{nSR} = cellfun(@(x,y) length(([x;y])),numLostNeighsAccum{nSR},numWonNeighsAccum{nSR});  

            infoOfCells = infoPerSurfaceRatio{nSR, 4};
            infoOfCells = infoOfCells{:};
            areaCells(nSR) = {infoOfCells.Basal_area};
            volumes(nSR) = {infoOfCells.Volume};
        end

        areaCellsPerSurfaceRealization = cat(2,areaCells{:});
        volumePerSurfaceRealization = cat(2,volumes{:});
        neighsSurface = cat(1,neighsSurface{:})';
        neighsAccumSurfaces = cat(1,neighsAccumSurfaces{:})';
        percentageScutoids = cat(1,percentageScutoids{:})';
        apicoBasalTransitions = cat(1,apicoBasalTransitions{:})';
        numLostNeighsAccum = cat(1, numLostNeighsAccum{:})';
        numWonNeighsAccum = cat(1, numWonNeighsAccum{:})';

        
end

