function [maskUpClean,maskDownClean,maskLeftClean,maskRightClean] = get4ProjectionsAlongMajorAxis(img3d,colours)
    
    [xLength,yLength,zLength] = size(img3d);

    closedImg3D = imdilate(img3d>0,strel('sphere',5));
    maskFill3d = imerode(bwmorph3(closedImg3D,'fill'),strel('sphere',5));
    
    setCentroids = zeros(zLength,3);

    %force Y axis as major axis
    if xLength > yLength
       img3d = permute(img3d,[2,1,3]); 
       maskFill3d = permute(maskFill3d,[2,1,3]); 
    end
    img3dUp = img3d;
    img3dDown = img3d;
    img3dLeft = img3d;
    img3dRight = img3d;
    
    for nY = 1 : size(maskFill3d,2)
        [x,~,z] = ind2sub(size(maskFill3d(:,nY,:)),find(maskFill3d(:,nY,:)));
        setCentroids(nY,:) = [round(mean(x)),nY,round(mean(z))];
        
        img3dUp (:,nY,1:setCentroids(nY,3)) = 0;
        img3dDown (:,nY,setCentroids(nY,3)+1:end) = 0;
        img3dLeft(1:setCentroids(nY,1),nY,:) = 0;
        img3dRight(setCentroids(nY,1)+1:end,nY,:) = 0;
    end
    
    %project Up
    maskUp = zeros(size(img3dUp,1),size(img3dUp,2));
    for nZ = size(img3dUp,3):-1:1
        sliceZ = img3dUp(:,:,nZ);
        maskUp(maskUp==0) = sliceZ(maskUp==0);
    end
    
    %project Down
    maskDown = zeros(size(img3dDown,1),size(img3dDown,2));
    for nZ = 1:size(img3dDown,3)
        sliceZ = img3dDown(:,:,nZ);
        maskDown(maskDown==0) = sliceZ(maskDown==0);
    end
    
    %project Left
    maskLeft = zeros(size(img3dLeft,2),size(img3dLeft,3));
    for nX = 1:size(img3dLeft,1)
        sliceX = img3dLeft(nX,:,:);
        maskLeft(maskLeft==0) = sliceX(maskLeft==0);
    end
    
    %project Right
    maskRight = zeros(size(img3dRight,2),size(img3dRight,3)); 
    for nX = size(img3dRight,1):-1:1
        sliceX = img3dRight(nX,:,:);
        maskRight(maskRight==0) = sliceX(maskRight==0);
    end
    
    

    %%clean projections   
    maskUpClean = cleanProjection(maskUp);
    maskDownClean = cleanProjection(maskDown);
    maskLeftClean = cleanProjection(maskLeft);
    maskRightClean = cleanProjection(maskRight);
    
    
%     colours = [0 0 0;colours];

%     figure,imshow(maskUp,colours)
%     figure,imshow(maskUpClean,colours)
%     
%     figure,imshow(maskDown,colours)
%     figure,imshow(maskDownClean,colours)
%     
%     figure,imshow(maskLeft,colours)
%     figure,imshow(maskLeftClean,colours)
%     
%     figure,imshow(maskRight,colours)
%     figure,imshow(maskRightClean,colours)
    
    
    
end

