%% Pipeline to unroll the divide the glands in SR and unroll the tubes to obtain them in 2D.
addpath(genpath('src'))
addpath(genpath('lib'))
addpath(genpath('gui'))

close all
clear all

files = dir('**/data/Salivary gland/**/Results/3d_layers_info.mat');
nonDiscardedFiles = cellfun(@(x) contains(lower(x), 'discarded') == 0 && contains(lower(x), 'e-cadh'), {files.folder});
files = files(nonDiscardedFiles);

disp('----------- UNROLLING TUBES DIVIDING GLANDS -------------')
for numFile = 1:length(files)
    files(numFile).folder
    selpath = files(numFile).folder;
    
    unroll_OnlyApicalAndBasal(selpath);
    calculate3DMorphologicalFeatures(files(numFile).folder);
    if contains(lower(files(numFile).folder), 'e-cadh') == 0 || contains(lower(files(numFile).folder), 'flatten')
        divideObjectInSurfaceRatios(selpath);
        unrollTube_parallel(selpath);
    end
end

