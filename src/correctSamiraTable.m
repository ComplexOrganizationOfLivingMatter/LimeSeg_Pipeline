addpath(genpath('src'))
addpath(genpath('lib'))
addpath(genpath('gui'))

close all
clear all

files = dir('**/data/Salivary gland_ToDivideInConstantPieces/**/Results/glandDividedInSurfaceRatios_AllUnrollFeatures.mat');

nonDiscardedFiles = cellfun(@(x) contains(lower(x), 'discarded') == 0 && contains(lower(x), 'wildtype'), {files.folder});
files = files(nonDiscardedFiles);
for numFile = 1:length(files)
    files(numFile).folder
    gland_SRs= dir(fullfile(files(numFile).folder, '**', 'samiraTable.mat'));
    
    for numSR= 1:length(gland_SRs)
        load(fullfile(gland_SRs(numSR).folder, gland_SRs(numSR).name));
        load(fullfile(gland_SRs(numSR).folder, 'allInfo.mat'), 'deployedImg3x');
        [allSamiraTable{numSR}, glandMinValue(numSR)] = normalizeVerticesGland(samiraTable,round(size(deployedImg3x,1)/2));
    end
    
    for numSR= 1:length(gland_SRs) 
        [normalizedSamiraTable{numSR}] = normalizeVerticesGland(allSamiraTable{numSR}, min(glandMinValue));
    end
    save(fullfile(gland_SRs(numSR).folder, gland_SRs(numSR).name), 'normalizedSamiraTable', '-append');
    
    %% Creating samira table
    samiraTableT = cell2table(normalizedSamiraTable, 'VariableNames',{'Radius', 'CellIDs', 'TipCells', 'BorderCell','verticesValues_x_y'});
    
    writetable(samiraTableT, strcat(selpath, '\', idName ,'_samirasFormat.xls'));
end