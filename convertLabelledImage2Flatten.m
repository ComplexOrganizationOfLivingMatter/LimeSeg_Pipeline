function flattenImage = convertLabelledImage2Flatten(labelledImage, lateralLayer)

    
    [row,col,z]=ind2sub(size(lateralLayer),find(lateralLayer>0));

    shp = alphaShape(col,row,z);   
%     pc = criticalAlpha(shp,'one-region');

    shp.Alpha = 120;
    
    idsCells = find(labelledImage>0);
    [qRow,qCol,qZ]=ind2sub(size(labelledImage),idsCells);
    
    tf = inShape(shp,qCol,qRow,qZ);
    idsCellIn = idsCells(tf);
    
    flattenImage=uint8(zeros(size(labelledImage)));
    flattenImage(idsCellIn)=labelledImage(idsCellIn);
    
end