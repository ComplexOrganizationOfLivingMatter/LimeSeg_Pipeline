function [labelledImage, basalLayer, apicalLayer, colours] = postprocessGland(labelledImage,outsideGland, lumenImage, outputDir, colours, tipValue)
%POSTPROCESSGLAND Summary of this function goes here
%   Detailed explanation goes here
%     [labelledImage] = fillEmptySpacesByWatershed3D(labelledImage, outsideGland | lumenImage, 1); % error outsideGland?
%     outsideGland_NotLumen = ~outsideGland | lumenImage;

    labelledImage = fill0sWithCells(labelledImage, labelledImage, outsideGland | lumenImage);
    labelledImage(lumenImage) = 0;

    %% Get basal layer by dilating the empty space
    basalLayer = getBasalFrom3DImage(labelledImage, lumenImage, tipValue, outsideGland & imdilate(lumenImage, strel('sphere', 1)) == 0);

    %% Get apical layer by dilating the lumen
    [apicalLayer] = getApicalFrom3DImage(lumenImage, labelledImage);
    exportAsImageSequence(apicalLayer, fullfile(outputDir, 'Apical_Labelled'), colours, tipValue);

    %% Export image sequence
    [colours] = exportAsImageSequence(labelledImage, fullfile(outputDir, 'Cells', 'labelledSequence', filesep), colours, tipValue);
end
