function tableDividedImages = interpolateImagesBySR(labelledImage,apicalLayer, basalLayer,finalSR, desiredSRs,path2save)
%% function designed to interpolate images at different heights between 2 surfaces

    %get closest apical coordinate from each basal coordinate
    allCells = unique(labelledImage(labelledImage>0));
    
    basalCoordCell = cell(max(allCells),1);
    closestApiCoordCell = cell(max(allCells),1);
    if exist(fullfile(path2save,'closestPaired_basalApical_coords.mat'),'file')==0
        for nCell = allCells'
            basalCell = basalLayer==nCell;
            [rowBas,colBas,zBas]=ind2sub(size(basalLayer),find(basalCell));

            apicalCell = apicalLayer==nCell;
            [rowApi,colApi,zApi]=ind2sub(size(apicalLayer),find(apicalCell));

            apicalIdsClosest=zeros(size(rowBas));
            for nId = 1:length(rowBas)
                distCoord = pdist2([colBas(nId),rowBas(nId), zBas(nId)],[colApi,rowApi, zApi]);
                [~,idClosest]=min(distCoord);
                apicalIdsClosest(nId) = idClosest;
            end

            basalCoordCell{nCell} = [rowBas,colBas,zBas];
            closestApiCoordCell{nCell} = [rowApi(apicalIdsClosest),colApi(apicalIdsClosest), zApi(apicalIdsClosest)];

        end

        save(fullfile(path2save,'closestPaired_basalApical_coords.mat'),'basalCoordCell','closestApiCoordCell')
    else
        load(fullfile(path2save,'closestPaired_basalApical_coords.mat'),'basalCoordCell','closestApiCoordCell')
    end
    
    %go over every desiredSR to get the interpolated coordinates  
    
    idsCells = find(labelledImage>0);
    [qRow,qCol,qZ]=ind2sub(size(basalLayer),idsCells);
    
    mkdir(fullfile(path2save,'dividedGlandBySr'))
    for selSR = 1:length(desiredSRs)
        
        distanceFactor = (desiredSRs(selSR)-1)/(finalSR-1);
        interpolatedCoord = cellfun(@(x,y) round(((x-y).*distanceFactor)+y) ,basalCoordCell,closestApiCoordCell,'UniformOutput',false);
        coordInterp = unique(vertcat(interpolatedCoord{:}),'rows');
        
        shp = alphaShape(coordInterp(:,2),coordInterp(:,1),coordInterp(:,3));
        pc = criticalAlpha(shp,'one-region');
        
        shp.Alpha = pc*2;
        tf = inShape(shp,qCol,qRow,qZ);
        
        idsCellIn = idsCells(tf);
        interpImageCells=uint8(zeros(size(labelledImage)));
        interpImageCells(idsCellIn)=labelledImage(idsCellIn);
        
        save(fullfile(path2save,'dividedGlandBySr',['sr_' num2str(desiredSRs(selSR)) '.mat']),'interpImageCells')
    end
end



