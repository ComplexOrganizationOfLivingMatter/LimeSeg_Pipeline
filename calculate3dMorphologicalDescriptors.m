close all
clear all

files = dir('**/data/Salivary gland/**/Results/3d_layers_info.mat');

totalMeanFeatures = [];
totalStdFeatures = [];

for numFiles=1:length(files)
load(fullfile(files(numFiles).folder, '3d_layers_info.mat'), 'labelledImage'); 
load(fullfile(files(numFiles).folder, 'valid_cells.mat'), 'validCells'); 

%% Extract each cell and calculate 3D features        
cells3dFeatures = regionprops3(labelledImage, 'PrincipalAxisLength', 'Volume', 'ConvexVolume', 'Solidity', 'SurfaceArea', 'EquivDiameter');
aspectRatio = max(cells3dFeatures.PrincipalAxisLength,[],2) ./ min(cells3dFeatures.PrincipalAxisLength,[],2);
sphereArea = 4 * pi .* ((cells3dFeatures.EquivDiameter) ./ 2) .^ 2;
sphericity = sphereArea ./ cells3dFeatures.SurfaceArea;

%% Save variables and export to excel
cells3dFeatures = [table(validCells', 'VariableNames', {'ID_Cell'}), cells3dFeatures(validCells,:), table(aspectRatio(validCells,:), sphericity(validCells,:), 'VariableNames', {'AspectRatio','Sphericity'})];
writetable(cells3dFeatures,fullfile(files(numFiles).folder,'3dFeatures_LimeSeg3DSegmentation.xls'), 'Range','B2');
save(fullfile(files(numFiles).folder, 'morphological3dFeatures.mat'), 'cells3dFeatures') 

%% Calculate mean and std of 3D features
meanFeatures = varfun(@(x) mean(x),cells3dFeatures);
stdFeatures = varfun(@(x) std(x),cells3dFeatures);

totalMeanFeatures = [totalMeanFeatures; meanFeatures];
totalStdFeatures = [totalStdFeatures; stdFeatures];

end

selpath = dir('**/data/Salivary gland/');
save(fullfile(selpath(1).folder, 'WT_3dFeatures.mat'), 'totalMeanFeatures','totalStdFeatures')
writetable([totalMeanFeatures; totalStdFeatures],fullfile(selpath(1).folder,'WT_3dFeatures.xls'), 'Range','B2');
