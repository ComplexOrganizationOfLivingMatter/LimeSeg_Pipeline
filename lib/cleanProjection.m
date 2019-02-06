function maskClean= cleanProjection(imgProj)

    maskClean = zeros(size(imgProj));
    areaCells = regionprops(imgProj,'Area');
    [ord,ind]=sort(cat(1,areaCells.Area),'descend');
    cellsImg = ind(ord>0)';
    for nCell = cellsImg
        binMask = imgProj==nCell;
        binMask = imfill(binMask,'holes');
        binMask = bwareaopen(binMask,2,8).*nCell;
        maskClean(maskClean==0)=binMask(maskClean==0);
    end

end

