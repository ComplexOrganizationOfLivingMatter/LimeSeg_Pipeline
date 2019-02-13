function unrollTube_parallel(selpath)
%UNROLLTUBE_PARALLEL Summary of this function goes here
%   Detailed explanation goes here
    
    load(fullfile(selpath, '3d_layers_info.mat'), 'labelledImage', 'apicalLayer', 'basalLayer', 'colours', 'glandOrientation', 'lumenImage');
    load(fullfile(selpath, 'valid_cells.mat'), 'validCells', 'noValidCells');

    apicalAreaValidCells = 100;
    disp('Apical');
    [apicalSamiraTable, apicalAreaValidCells] = unrollTube(apicalLayer, fullfile(selpath, 'apical'), noValidCells, colours);
    
    disp('Basal');
    basalSamiraTable = unrollTube(basalLayer, fullfile(selpath, 'basal'), noValidCells, colours, apicalAreaValidCells);
    
    samiraTable = [apicalSamiraTable; basalSamiraTable];
    
    idName_splitted = strsplit(selpath, filesep);
    
    idName = strjoin(idName_splitted(end-3:end-1), '_');
    samiraTableT = cell2table(samiraTable, 'VariableNames',{'Radius', 'CellIDs', 'TipCells', 'BorderCell','verticesValues_x_y'});

    newCrossesTable = lookFor4cellsJunctionsAndExportTheExcel(samiraTableT);

    writetable(samiraTableT, strcat(dir2save, '\', idName ,'_samirasFormat_', date, '.xls'));
    writetable(newCrossesTable, strcat(dir2save, '\', idName ,'_VertCrosses_', date, '.xls'));
end

