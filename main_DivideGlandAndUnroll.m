
addpath(genpath('src'))
addpath(genpath('lib'))
addpath(genpath('gui'))

close all
clear all

files = dir('**/data/Salivary gland/**/Results/3d_layers_info.mat');
nonDiscardedFiles = cellfun(@(x) contains(lower(x), 'discarded') == 0 && contains(lower(x), 'wildtype'), {files.folder});
files = files(nonDiscardedFiles);

disp('----------- DIVIDING GLANDS -------------')
parfor numFile = 1:length(files)
    files(numFile).folder
    selpath = files(numFile).folder;
    unroll_OnlyApicalAndBasal(selpath)
    divideObjectInSurfaceRatios(selpath);
end

disp('----------- UNROLLING TUBES -------------')
parfor numFile = 1:13%length(files)
    files(numFile).folder
    selpath = files(numFile).folder;
    
    unrollTube_parallel(selpath);
end

