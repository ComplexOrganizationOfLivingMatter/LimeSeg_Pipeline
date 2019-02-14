function unrollTube_parallel(selpath)
%UNROLLTUBE_PARALLEL Summary of this function goes here
%   Detailed explanation goes here

    idName_splitted = strsplit(selpath, filesep);
    idName = strjoin(idName_splitted(end-3:end-1), '_');
    
    if exist(strcat(selpath, '\', idName ,'_samirasFormat.xls'), 'file') == 0
        
        load(fullfile(selpath, '3d_layers_info.mat'), 'labelledImage', 'apicalLayer', 'basalLayer', 'colours', 'glandOrientation', 'lumenImage');
        load(fullfile(selpath, 'valid_cells.mat'), 'validCells', 'noValidCells');

        apicalAreaValidCells = 100;
        disp('Apical');
        [apicalSamiraTable, apicalAreaValidCells] = unrollTube(apicalLayer, fullfile(selpath, 'apical'), noValidCells, colours);

        disp('Basal');
        basalSamiraTable = unrollTube(basalLayer, fullfile(selpath, 'basal'), noValidCells, colours, apicalAreaValidCells);

        samiraTable = [apicalSamiraTable; basalSamiraTable];


        samiraTableT = cell2table(samiraTable, 'VariableNames',{'Radius', 'CellIDs', 'TipCells', 'BorderCell','verticesValues_x_y'});

        newCrossesTable = lookFor4cellsJunctionsAndExportTheExcel(samiraTableT);

        writetable(samiraTableT, strcat(selpath, '\', idName ,'_samirasFormat.xls'));
        writetable(newCrossesTable, strcat(selpath, '\', idName ,'_VertCrosses.xls'));
    end
end

