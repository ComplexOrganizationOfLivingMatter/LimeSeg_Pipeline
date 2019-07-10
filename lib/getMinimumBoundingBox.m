function [rect] = getMinimumBoundingBox(regionPropsOfObjects)
%GETMINIMUMBOUNDINGBOX Summary of this function goes here
%   Detailed explanation goes here
%// Obtain all of the bounding box co-ordinates
    bboxCoords = reshape([regionPropsOfObjects.BoundingBox], 4, []).';

    % // Calculate top left corner
    topLeftCoords = bboxCoords(:,1:2);

    % // Calculate top right corner
    topRightCoords = [topLeftCoords(:,1) + bboxCoords(:,3) topLeftCoords(:,2)];

    % // Calculate bottom left corner
    bottomLeftCoords = [topLeftCoords(:,1) topLeftCoords(:,2) + bboxCoords(:,4)];

    % // Calculate bottom right corner
    bottomRightCoords = [topLeftCoords(:,1) + bboxCoords(:,3) ...
        topLeftCoords(:,2) + bboxCoords(:,4)];

    % // Calculating the minimum and maximum X and Y values
    finalCoords = [topLeftCoords; topRightCoords; bottomLeftCoords; bottomRightCoords];
    minX = min(finalCoords(:,1));
    maxX = max(finalCoords(:,1));
    minY = min(finalCoords(:,2));
    maxY = max(finalCoords(:,2));

    % Draw the rectangle on the screen
    width = (maxX - minX);
    height = (maxY - minY);
    rect = [minX minY width height];
end

