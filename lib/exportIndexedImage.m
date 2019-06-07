function [labelledImageFinal] = exportIndexedImage(labelledImage, initialSize, tipValue, colours)
%EXPORTINDEXEDIMAGE Summary of this function goes here
%   Detailed explanation goes here

    labelledImageNoBorders = labelledImage((tipValue+2):end-(tipValue+1), (tipValue+2):end-(tipValue+1), (tipValue+2):end-(tipValue+1));
    labelledImageFinal = imresize3(labelledImageNoBorders, initialSize, 'nearest')+1;
    
    %Here 1 is the background
    for numZ = 1:initialSize(3)
        imwrite(labelledImageFinal(:, :, numZ), [1 1 1; colours], strcat('tmp/imgZCoord_', num2str(numZ),'.tif'));
    end
end

