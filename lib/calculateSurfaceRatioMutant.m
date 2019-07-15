function [realSurfaceRatio] = calculateSurfaceRatioMutant(basalLayer, apicalLayer, labelledImage, validCells)
%CALCULATESURFACERATIOMUNTANT Summary of this function goes here
%   We ought to get the surface ratio (RadiusBasal/RadiusApical). In this
%   case we want to get RadiusBasal as a cylindre amid two basal surfaces:
%   the one forming the outer gland (the maximum surface formed by the max
%   cell height of each cell) and an inner basal surface formed by the
%   minimum points of the cells intersecting.
    
    %% Get real basal layer
    [apicalLayer,basalLayer] = flattenMutantGland(apicalLayer, basalLayer, labelledImage);
    
    %% Calculate real Surface Ratio
    apicalAreaCells=regionprops(apicalLayer,'Area');
    basalAreaCells=regionprops(basalLayer,'Area');
    apicalAreaCells = apicalAreaCells(validCells);
    basalAreaCells = basalAreaCells(validCells);
    realSurfaceRatio = mean([basalAreaCells.Area]) / mean([apicalAreaCells.Area]);
    
end

