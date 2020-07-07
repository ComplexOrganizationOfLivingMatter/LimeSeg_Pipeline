function [orderedLabels] = perim2line(filledImage, img3d, img3dComplete, coordZ)
%PERIM2LINE Summary of this function goes here
%   Detailed explanation goes here
    finalPerimImage = bwperim(filledImage);
    solidityOfObjects = regionprops(filledImage>0, 'Solidity');
    %imshow(filledImage)
    solidityThreshold = 0.6;
    %Check if there is a hole
    if length(solidityOfObjects) == 1
        goCalculatePerim = solidityOfObjects.Solidity < solidityThreshold;
    else
        goCalculatePerim = 1;
    end
    if goCalculatePerim
%                 convexPerimImage = regionprops(imclose(finalPerimImage, strel('disk', 5)), 'convexHull');
%                 convexPerimImage = convexPerimImage.ConvexHull;
%
%                 validRegion = zeros(size(finalPerimImage));
%                 [xq, yq] = find(validRegion==0);
%                 in = inpolygon(xq,yq, round(convexPerimImage(:, 2)), round(convexPerimImage(:, 1)));
%
%                 indicesInsideCell = sub2ind(size(finalPerimImage), xq, yq);
%
%                 validRegion(indicesInsideCell(in)) = 1;
%
%                 finalPerimImage = bwperim(validRegion);


        finalPerimImage = uint16(bwskel(filledImage>0));
        %fill0sWithCells(img3d(:, :, coordZ) ,validRegion);

        [X,Y] = meshgrid(1:size(filledImage,2), 1:size(filledImage,1));

        s = regionprops(filledImage>0, 'BoundingBox');
        [rect] = getMinimumBoundingBox(s);

        %bb = floor(s.BoundingBox); %// Could be floating point, so floor it
        bb = rect;
        cenx = bb(1) + (bb(3) / 2.0); %// Get the centre of the bounding box
        ceny = bb(2) + (bb(4) / 2.0);

        radi = max(bb(3), bb(4)) / 2; %// Find the best radius
        perimeterNew = ((X - cenx).^2 + (Y - ceny).^2) <= radi^2; %// Draw our circle and place in a temp. image
        perimeterNew = bwperim(perimeterNew); %// Add this circle on top of our output image
        %figure; imshow(perimeterNew)
    end

    %imshow(finalPerimImage)

    %% Obtaining the center of the cylinder
    [x, y] = find(finalPerimImage > 0);
    centroidCoordZ = mean([x, y], 1); % Centroid of each real Y of the cylinder
    centroidX = centroidCoordZ(1);
    centroidY = centroidCoordZ(2);

    [x, y] = find(finalPerimImage > 0);

    %% labelled mask
    if goCalculatePerim
        centroidX = ceny;
        centroidY = cenx;

        [x_PerimeterNew, y_PerimeterNew] = find(perimeterNew > 0);
        angleLabelCoord_NewPerimeter = atan2(y_PerimeterNew - centroidY, x_PerimeterNew - centroidX);
        angleLabelCoord_NewPerimeter_Sorted = sort(angleLabelCoord_NewPerimeter);
        minDistance = abs(angleLabelCoord_NewPerimeter_Sorted(1) - angleLabelCoord_NewPerimeter_Sorted(2));
        img3dPerimFilled = fill0sWithCells(img3d(:, :, coordZ), img3dComplete(:, :, coordZ),filledImage==0);
        maskLabel=uint16(finalPerimImage).*img3dPerimFilled;
    else
        if isequal(filledImage, double(img3d(:, :, coordZ)>0))==0
            img3dPerimFilled = fill0sWithCells(img3d(:, :, coordZ), img3dComplete(:, :, coordZ),filledImage==0);
        else
            img3dPerimFilled = img3d(:, :, coordZ);
        end
        maskLabel=uint16(finalPerimImage).*img3dPerimFilled;
    end
    %angles label coord regarding centroid
    angleLabelCoord = atan2(y - centroidY, x - centroidX);
    [angleLabelCoordSort, orderedIndices] = sort(angleLabelCoord);

    %% Completing the missing parts of the circle perim
    if goCalculatePerim
        distanceToNextPoint = angleLabelCoordSort([2:end 1]) - angleLabelCoordSort;
        distanceToNextPoint(end) = distanceToNextPoint(end) + 6;
        if max(distanceToNextPoint) > minDistance*3
            [~, positionsToFill] = max(distanceToNextPoint);
            if positionsToFill+1 > length(distanceToNextPoint)
                newAngles = angleLabelCoord_NewPerimeter_Sorted(angleLabelCoordSort(1) > angleLabelCoord_NewPerimeter_Sorted);

                angleLabelCoordSort = [newAngles; angleLabelCoordSort];
                orderedIndices = [zeros(size(newAngles)); orderedIndices];
            else
                newAngles = angleLabelCoord_NewPerimeter_Sorted(angleLabelCoordSort(positionsToFill) < angleLabelCoord_NewPerimeter_Sorted & angleLabelCoordSort(positionsToFill+1) > angleLabelCoord_NewPerimeter_Sorted);
                angleLabelCoordSort = [angleLabelCoordSort(1:positionsToFill); newAngles; angleLabelCoordSort(positionsToFill+1:end)];
                orderedIndices = [orderedIndices(1:positionsToFill); zeros(size(newAngles)); orderedIndices(positionsToFill+1:end)];
            end
        end
    end

    %% Assing label to pixels of perimeters
    %If a perimeter coordinate have no label pixels in a range of pi/45 radians, it label is 0
    orderedLabels = zeros(1,length(angleLabelCoordSort));
    for nCoord = 1:length(angleLabelCoordSort)
        if orderedIndices(nCoord) ~= 0
            indicesClosest = sub2ind(size(maskLabel), x(orderedIndices(nCoord)), y(orderedIndices(nCoord)));
            pixelLabel = maskLabel(indicesClosest);
            orderedLabels(nCoord) = pixelLabel;
        else
            orderedLabels(nCoord) = 0;
        end
    end
end

