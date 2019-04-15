function [image_filled] = fillGlandWithHole(initialImage)
%FILLGLANDWITHHOLE Summary of this function goes here
%   Detailed explanation goes here
    image_closed_initial = imclose(initialImage>0, strel('sphere', 2));

    image_closed = image_closed_initial;
    [~, y, ~] = ind2sub(size(initialImage),find(initialImage>0));

    downSide = initialImage(:, min(y), :);
    upSide = initialImage(:, max(y), :);
    if sum(downSide(:)>0) > sum(upSide(:)>0)
        image_closed(:, 1:(min(y)+3), :) = 1;
    else
        image_closed(:, (max(y)-3):end, :) = 1;
    end
    %image_closed(lumenImage) = 1;
    image_filled = imfill(double(image_closed), 18, 'holes');
    if sum(downSide(:)>0) > sum(upSide(:)>0)
        image_filled(:, 1:(min(y)+3), :) = image_closed_initial(:, 1:(min(y)+3), :);
    else
        image_filled(:, (max(y)-3):end, :) = image_closed_initial(:, (max(y)-3):end, :);
    end
end

