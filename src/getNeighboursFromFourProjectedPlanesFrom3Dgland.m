function uniqNeigh = getNeighboursFromFourProjectedPlanesFrom3Dgland(img3d,colours)

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

    [upProjA,downProjA,leftProjA,rightProjA] = get4ProjectionsAlongMajorAxis(img3d_rotACrop,colours);
    radiusNeigh = 1;
    
    maxCell = max(img3d(:));
    
    neighUpA = calculateNeighbours(upProjA,radiusNeigh);
    neighUpA = [neighUpA,cell(1,maxCell-length(neighUpA))];
    neighDownA = calculateNeighbours(downProjA,radiusNeigh);
    neighDownA = [neighDownA,cell(1,maxCell-length(neighDownA))];
    neighLeftA = calculateNeighbours(leftProjA,radiusNeigh);
    neighLeftA = [neighLeftA,cell(1,maxCell-length(neighLeftA))];
    neighRightA = calculateNeighbours(rightProjA,radiusNeigh);
    neighRightA = [neighRightA,cell(1,maxCell-length(neighRightA))];

    [upProjB,downProjB,leftProjB,rightProjB] = get4ProjectionsAlongMajorAxis(img3d_rotBCrop,colours);
    neighUpB = calculateNeighbours(upProjB,radiusNeigh);
    neighUpB = [neighUpB,cell(1,maxCell-length(neighUpB))];
    neighDownB = calculateNeighbours(downProjB,radiusNeigh);
    neighDownB = [neighDownB,cell(1,maxCell-length(neighDownB))];
    neighLeftB = calculateNeighbours(leftProjB,radiusNeigh);
    neighLeftB = [neighLeftB,cell(1,maxCell-length(neighLeftB))];
    neighRightB = calculateNeighbours(rightProjB,radiusNeigh);
    neighRightB = [neighRightB,cell(1,maxCell-length(neighRightB))];
    
    %unify neighs
    uniqNeigh = cellfun(@(a,b,c,d,e,f,g,h) unique([a;b;c;d;e;f;g;h]),neighUpA,neighDownA,...
        neighLeftA,neighRightA,neighUpB,neighDownB,neighLeftB,neighRightB,'UniformOutput',false);
    
    

end