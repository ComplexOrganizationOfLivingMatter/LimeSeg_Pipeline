function realisticSR = calculateRealisticSR(meanCellProperties, apicalLayer,lumenSkeleton, validCells,pixelScale,folderFeatures)


    centroidsApical = table2array(regionprops3(apicalLayer,'Centroid'));
    centroidsApicalAux = centroidsApical;
    
    lumenRadius2cell= nan(size(centroidsApical,1),1);
    
    %get all lumen skel voxels
    idsSkel = find(lumenSkeleton);
    [rowSkel, colSkel, zSkel] = ind2sub(size(lumenSkeleton),idsSkel);
    
    for idCell = 1:length(lumenRadius2cell)
        if ~isnan(centroidsApical(idCell,1))
            %get all apical cell voxels
            idsApical = find(apicalLayer==idCell);
            [rowApical, colApical, zApical] = ind2sub(size(apicalLayer),idsApical);
            
            %calculate the closest apical pixel to the calculated apical centroid
            distCoord = pdist2([colApical,rowApical, zApical],centroidsApical(idCell,:));
            [~,idSeedMin]=min(distCoord);
            centroidsApicalAux(idCell,:)= [colApical(idSeedMin),rowApical(idSeedMin), zApical(idSeedMin)];
            
            %measure distance between the apical centroid and its closest skeleton voxel
            distCoord = pdist2(centroidsApicalAux(idCell,:),[colSkel, rowSkel, zSkel]);
            lumenRadius2cell(idCell) = min(distCoord);
        end
    end
    
    meanLumenDiameter = mean(lumenRadius2cell(validCells))*pixelScale;
    realisticSR = (meanCellProperties.Fun_cell_height+meanLumenDiameter)/meanLumenDiameter;
    
    save(fullfile(folderFeatures,'realisticSR.mat'),'realisticSR')

end
