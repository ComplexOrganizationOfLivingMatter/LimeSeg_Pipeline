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
        [samiraTable] = normalizeVerticesGland(samiraTable,deployedImg3x);
    end
    
    
end