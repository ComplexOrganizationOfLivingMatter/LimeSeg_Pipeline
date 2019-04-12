function unroll_OnlyApicalAndBasal(selpath)
%UNROLL_ONLYAPICALANDBASAL Summary of this function goes here
%   Detailed explanation goes here

    %% ONLY APICAL AND BASAL
    load(fullfile(selpath, '3d_layers_info.mat'), 'labelledImage', 'apicalLayer', 'basalLayer', 'colours', 'glandOrientation', 'lumenImage');
    load(fullfile(selpath, 'valid_cells.mat'), 'validCells', 'noValidCells');
    
    apicalAreaValidCells = 100;
    disp('Apical');
    mkdir(fullfile(selpath, 'unrolledGlands', 'gland_SR_1'));
    [~, apicalAreaValidCells] = unrollTube(apicalLayer, fullfile(selpath, 'unrolledGlands', 'gland_SR_1'), noValidCells, colours);

    disp('Basal');
    mkdir(fullfile(selpath, 'unrolledGlands', 'gland_SR_basal'));
    unrollTube(basalLayer, fullfile(selpath, 'unrolledGlands', 'gland_SR_basal'), noValidCells, colours, apicalAreaValidCells);
end

