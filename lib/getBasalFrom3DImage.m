function [basalLayer] = getBasalFrom3DImage(labelledImage, lumenImage, tipValue)
%GETBASALFROM3DIMAGE Summary of this function goes here
%   Detailed explanation goes here
    se = strel('sphere',tipValue);
    objectDilated = imdilate(labelledImage>0, se);
    objectDilated = imfill(objectDilated, 'holes');
    finalObject = imerode(objectDilated, se);
    finalObject = bwareaopen(finalObject, 5);
%     [x,y,z] = ind2sub(size(finalObject),find(finalObject>0));
%     figure;
%     pcshow([x,y,z]);
    
    se = strel('sphere', 1);
    finalObjectEroded = imerode(finalObject, se);
    basalLayer = finalObject - finalObjectEroded;
    
    
    [~,y,~] = ind2sub(size(basalLayer),find(basalLayer>0));
    
    basalLayer(:, :, end) = finalObject(:, :, end);
    basalLayer(:, :, 1) = finalObject(:, :, 1);
    downSide = labelledImage(:, min(y), :);
    upSide = labelledImage(:, max(y), :);
    if sum(downSide(:)>0) > sum(upSide(:)>0)
        basalLayer(:, 1:min(y), :) = 0;
    else
        basalLayer(:, max(y):end, :) = 0;
    end
    
%     figure;
%     pcshow([x,y,z]);
    if isempty(lumenImage) == 0
        glandParts = bwlabeln(basalLayer>0);
        glandParts_Unique = unique(glandParts);
        
        if length(glandParts_Unique)>2
            basalLayer(imdilate(lumenImage, strel('sphere', 3))) = 0;
        end
    end
    %basalLayer = completeImageOfCells(labelledImage .* basalLayer, basalLayer == 0);
    basalLayer = labelledImage .* basalLayer;
end

