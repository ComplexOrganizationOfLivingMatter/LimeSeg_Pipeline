function [filledImage] = createCompleteSection(img3d, coordZ, labelledImage_realSize)
%CREATECOMPLETESECTION Summary of this function goes here
%   Detailed explanation goes here
    filledImage = imfill(double(img3d(:, :, coordZ)>0), 'holes');
    %             imshow(img3d(:, :, coordZ), colours)
    if exist('labelledImage_realSize', 'var') == 0
        filledImage = bwareafilt(filledImage>0, 1, 4);
    elseif isequal(filledImage, double(img3d(:, :, coordZ)>0)) || length(unique(bwlabel(filledImage)))>2
%                 imageLabelled = bwlabel(filledImage);
%                 infoOfObjects = regionprops(imageLabelled, 'Area');
% %                 infoOfObjects = regionprops(imageLabelled, {'Area', 'Eccentricity', 'Solidity'});
% %                 solidityOfObjects = [infoOfObjects.Solidity];
% %                 circularityOfObjects = [infoOfObjects.Eccentricity];
%                 areaOfObjects = [infoOfObjects.Area];
%                 [~, biggerArea] = max(areaOfObjects);
%                 if sum(areaOfObjects>30) > 1
%                     imgToChange = imfill(double(imclose(filledImage>0, strel('disk', 5))), 'holes');
%                     imgToChangedist = bwdist(~imgToChange);
%                     imgToChangedist = -imgToChangedist;
%                     imgToChangedist(~imgToChange) = Inf;
%                     imgToChangeWS = watershed(imgToChangedist);
%                     imgToChangeWS(~imgToChange) = 0;
%
%                     objectsToKeep = bwareafilt(imgToChangeWS>0, 1);
%                     filledImage = objectsToKeep & filledImage;
%                     %imshow(filledImageToTest);
%                     %imshow(filledImage);
%                 end
        [x, y] = find(img3d(:, :, coordZ)>0);
        coordinates = [x, y];

        userConfig = struct('xy',coordinates, 'showProg',false,'showResult',false);
        resultStruct = tspo_ga(userConfig);
        orderBoundary = [resultStruct.optRoute resultStruct.optRoute(1)];
        newVertSalesman = coordinates(orderBoundary(1:end-1), :);
        newVertSalesman = [newVertSalesman; newVertSalesman(1,:)];

        %                 for numCoord = 1:(size(newVertSalesman, 1)-1)
        %                     plot(newVertSalesman(numCoord:numCoord+1, 2), newVertSalesman(numCoord:numCoord+1, 1))
        %                 end

        [distancePxs] = pdist2(coordinates, coordinates, 'euclidean', 'Smallest', 4);
        coordinatesToUnify = coordinates(distancePxs(4, :) > 1.5, :);

        %                 figure; imshow(filledImage)
        filledImageAux = filledImage>2;
        for numCoord = 1:size(coordinatesToUnify, 1)

            midVertex = coordinatesToUnify(numCoord, :);
            midvertexIndex = find(ismember(newVertSalesman, midVertex, 'rows'));
            if length(midvertexIndex)==1
                initVertex = newVertSalesman(midvertexIndex-1, :);
                endVertex = newVertSalesman(midvertexIndex+1, :);
            else
                initVertex = newVertSalesman(2, :);
                endVertex = newVertSalesman(end-1, :);
            end

            [xnAcum1, ynAcum1] = Drawline3D(initVertex(1), initVertex(2), 0, midVertex(1), midVertex(2), 0);
            [xnAcum2, ynAcum2] = Drawline3D(midVertex(1), midVertex(2), 0, endVertex(1), endVertex(2), 0);

            xnAcum = [xnAcum1; xnAcum2];
            ynAcum = [ynAcum1; ynAcum2];
            %hold on;
            %plot(coordinatesToUnify(numCoord, 2), coordinatesToUnify(numCoord, 1), 'r*')
            indicesToSave = sub2ind(size(img3d(:, :, coordZ)), xnAcum, ynAcum);
            filledImageAux(indicesToSave) = 1;
        end
        filledImagenew = imfill(filledImageAux | filledImage, 'holes');
        if isequal(filledImage, filledImagenew) == 0
            filledImage = filledImagenew;
            filledImage = bwareafilt(filledImage>0, 1, 4);
        else
            disp('Warning');
        end
    else
        filledImage = bwareafilt(filledImage>0, 1);
    end
end

