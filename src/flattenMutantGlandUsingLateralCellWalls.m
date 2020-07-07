function labelledImageRealSR = flattenMutantGlandUsingLateralCellWalls(labelledImage, lumenImage)
    
    pixel2dilate=2;
    %%dilate each cell and save indexes
    se = strel('sphere', pixel2dilate);
    cells = unique(labelledImage);
    cells = cells(cells~=0)';
    pixelsDilatedPerCells = cell(1,max(cells));
    for numCell = cells
        perimCell = bwperim(labelledImage==numCell);
        perimCellDilated = imdilate(perimCell,se);
        idsPerimDilated = find(perimCellDilated>0);
        pixelsDilatedPerCells{numCell} = idsPerimDilated;
    end
    
    totalPixelsIds = sort([vertcat(pixelsDilatedPerCells{:})]);
    uniqPixels = sort(unique(totalPixelsIds)); % which will give you the unique elements of A in array B
    Ncount = histc(totalPixelsIds,uniqPixels); % this willgive the number of occurences of each unique element
    %overlaping pixels between several cells
    idRepeated = uniqPixels(Ncount>1);
    
    labelledImageToCalculateSR = zeros(size(labelledImage));  
    
    for numCell = cells
        try 
            indicesCell = find(labelledImage == numCell);
            [xC,yC,zC]=ind2sub(size(labelledImage),indicesCell);

            indicesCellWall = indicesCell(ismember(indicesCell,idRepeated));
            [xCW,yCW,zCW]=ind2sub(size(labelledImage),indicesCellWall);
            plotCell = alphaShape([xCW,yCW,zCW], 1);
            %enormous alpha radius. It doesnt matter the number while the
            %cell si close
            plotCell.Alpha = 500;

            pixelsIn = inShape(plotCell, xC,yC,zC);
            labelledImageToCalculateSR(indicesCell(pixelsIn)) = numCell;
        catch
            %try/catch for non-valid cells with strange shape (avoid flatten)
            ['Cell ' num2str(numCell) ' without flatten']
            labelledImageToCalculateSR(labelledImage==numCell)=numCell;
        end
    end
    labelledImageRealSR = uint16(imclose(labelledImageToCalculateSR > 0,strel('sphere',2)));
    labelledImageRealSR = labelledImageRealSR .* labelledImage;
    labelledImageRealSR(lumenImage>0) = 0;
    
end

