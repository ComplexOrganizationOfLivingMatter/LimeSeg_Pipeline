function getFourProjectedPlanesFrom3Dgland(img3d,colours)

    %crop img3d
    [allX,allY,allZ]=ind2sub(size(img3d),find(img3d>0));
    img3dCrop = img3d(min(allX):max(allX),min(allY):max(allY),min(allZ):max(allZ));
    
    %divide image in 2 parts and rotate each subset to minimize the projection
    %error due to the gland orientation
    [xSize,ySize,~] = size(img3dCrop);
    if xSize > ySize
        img3d_section1 = img3dCrop(1:round(xSize/2 + xSize/8),:,:);
        img3d_section2 = img3dCrop(round(xSize/2 - xSize/8):end,:,:);
    else
        img3d_section1 = permute(img3dCrop(:,1:round(ySize/2 + ySize/8),:),[2,1,3]);
        img3d_section2 = permute(img3dCrop(:,round(ySize/2 - ySize/8):end,:),[2,1,3]);        
    end

    img3d_rotA = rotateImg3(img3d_section1);
    [allX,allY,allZ]=ind2sub(size(img3d_rotA),find(img3d_rotA>0));
    img3d_rotACrop = img3d_rotA(min(allX):max(allX),min(allY):max(allY),min(allZ):max(allZ));
    
    img3d_rotB = rotateImg3(img3d_section2);
    [allX,allY,allZ]=ind2sub(size(img3d_rotB),find(img3d_rotB>0));
    img3d_rotBCrop = img3d_rotB(min(allX):max(allX),min(allY):max(allY),min(allZ):max(allZ));

    get4ProjectionsAlongMajorAxis(img3d_rotACrop,colours)
    get4ProjectionsAlongMajorAxis(img3d_rotBCrop,colours)

    

end