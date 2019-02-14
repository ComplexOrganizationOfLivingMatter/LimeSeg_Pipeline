function [apicalLayer] = getApicalFrom3DImage(lumenImage, labelledImage)
%GETAPICALFRO3DIMAGE Summary of this function goes here
%   Detailed explanation goes here
    se = strel('sphere', 1);
    lumenImage = imclose(lumenImage, strel('sphere', 5));
    dilatedLumen = imdilate(lumenImage, se);
    apical1Pixel = (dilatedLumen - lumenImage);
    apicalLayer = labelledImage .* apical1Pixel;
%     [x,y,z] = ind2sub(size(apicalLayer),find(apicalLayer>0));
%     figure;
%     pcshow([x,y,z]);
end

