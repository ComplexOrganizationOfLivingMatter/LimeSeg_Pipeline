addpath(genpath('src'))
addpath(genpath('lib'))
addpath(genpath('gui'))

close all
clear all

files = dir('**/data/Salivary gland_ToDivideInConstantPieces/**/Results/glandDividedInSurfaceRatios_AllUnrollFeatures.mat');

nonDiscardedFiles = cellfun(@(x) contains(lower(x), 'discarded') == 0 && contains(lower(x), 'wildtype'), {files.folder});
files = files(nonDiscardedFiles);

minNumberOfSurfaceRatios = 7;
steps = 2.5/(minNumberOfSurfaceRatios-1);
surfaceRatiosExtrapolatedFrom3D = 1:steps:((steps*(minNumberOfSurfaceRatios-1))+1);

for numFile = 1:length(files)
    files(numFile).folder
    gland_SRs= dir(fullfile(files(numFile).folder, '**', 'samiraTable.mat'));
    
    allSamiraTable=cell(minNumberOfSurfaceRatios,1);
    glandMinValue=zeros(minNumberOfSurfaceRatios,1);
    
    for numSR= 1:minNumberOfSurfaceRatios
        load(fullfile(gland_SRs(numSR).folder, gland_SRs(numSR).name));
        load(fullfile(gland_SRs(numSR).folder, 'allInfo.mat'), 'deployedImg3x');
        %round(size(deployedImg3x,2)/2)
        [allSamiraTable{numSR}, glandMinValue(numSR)] = normalizeVerticesGland(samiraTable,round(size(deployedImg3x,2)/2));
    end
    
    %plotVerticesPerSurfaceRatio(vertcat(allSamiraTable{:}));
    
    load(fullfile(files(numFile).folder,'glandDividedInSurfaceRatios_AllUnrollFeatures.mat'));
    normalizedSamiraTable=cell(minNumberOfSurfaceRatios, 1);
    
    for numSR= 1:minNumberOfSurfaceRatios
        [normalizedSamiraTable{numSR}] = normalizeVerticesGland(allSamiraTable{numSR}, min(glandMinValue));
        actualSamiraTable=normalizedSamiraTable{numSR};
        actualSamiraTable(:,1)={surfaceRatiosExtrapolatedFrom3D(numSR)};
        normalizedSamiraTable{numSR} = actualSamiraTable;
    end
    save(fullfile(gland_SRs(numSR).folder, gland_SRs(numSR).name), 'normalizedSamiraTable', '-append');
       
    %plotVerticesPerSurfaceRatio(vertcat(normalizedSamiraTable{:}));
    
    %% Creating samira table
    normalizedSamiraTable = vertcat(normalizedSamiraTable{:});
    samiraTableT = cell2table(normalizedSamiraTable, 'VariableNames',{'Radius', 'CellIDs', 'TipCells', 'BorderCell','verticesValues_x_y'});
    idName_splitted = strsplit(files(numFile).folder, filesep);
    idName = strjoin(idName_splitted(end-3:end-1), '_');
    writetable(samiraTableT, fullfile(files(numFile).folder, strcat(idName ,'_samirasFormat.xls')));
end