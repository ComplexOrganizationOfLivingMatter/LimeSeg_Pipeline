close all
clear all

files = dir('**/data/Salivary gland/**/Results/3d_layers_info.mat');

for numFiles=1:length(files)
load(fullfile(files(numFiles).folder, '3d_layers_info.mat'), 'labelledImage'); 
load(fullfile(files(numFiles).folder, 'valid_cells.mat'), 'noValidCells'); 

morphological3dFeatures={};
totalCellPrincipleAxis={};

%% Extract each cell and calculate 3D features
for numCell=1:max(max(max(labelledImage)))
        labelledCell=zeros(size(labelledImage,1),size(labelledImage,2),size(labelledImage,3));
        labelledCell=(labelledImage==numCell);
        
        cellPrincipleAxis = table2array(regionprops3(labelledCell, 'PrincipalAxisLength'));
        aspectRatio = max(cellPrincipleAxis)/min(cellPrincipleAxis);

        volume=table2array(regionprops3(labelledCell, 'Volume'));
        convexVolume = table2array(regionprops3(labelledCell, 'ConvexVolume'));
        solidity = table2array(regionprops3(labelledCell, 'Solidity'));

        surfaceArea = table2array(regionprops3(labelledCell, 'SurfaceArea'));
        sphereDiameter = table2array(regionprops3(labelledCell, 'EquivDiameter'));
        sphereArea = 4*pi*((sphereDiameter)/2)^2;
        sphericity = sphereArea/surfaceArea;
        
        morphological3dFeatures{numCell,1} = numCell;
        morphological3dFeatures{numCell,2} = aspectRatio;
        morphological3dFeatures{numCell,3} = volume;
        morphological3dFeatures{numCell,4} = convexVolume;
        morphological3dFeatures{numCell,5} = solidity;
        morphological3dFeatures{numCell,6} = surfaceArea;
        morphological3dFeatures{numCell,7} = sphericity;

        totalCellPrincipleAxis{numCell,1} = cellPrincipleAxis;
        
        cellPrincipleAxis = [];
        aspectRatio = [];
        volume= [];
        convexVolume = [];
        solidity = [];
        surfaceArea = [];
        sphereDiameter = [];
        sphereArea = [];
        sphericity = [];
end

morphological3dFeatures(noValidCells,:)=[];
totalCellPrincipleAxis(noValidCells,:)=[];

%% Save variables
totalAspectRatio = {morphological3dFeatures{:,2}};
totalVolume = {morphological3dFeatures{:,3}};
totalConvexVolume = {morphological3dFeatures{:,4}};
totalSolidity = {morphological3dFeatures{:,5}};
totalSurfaceArea = {morphological3dFeatures{:,6}};
totalSphericity = {morphological3dFeatures{:,7}};

save(fullfile(files(numFiles).folder, 'morphological3dFeatures.mat'), 'totalAspectRatio','totalVolume','totalConvexVolume','totalSolidity','totalSurfaceArea','totalSphericity','totalCellPrincipleAxis'); 

%% Export to excel
xls3dFeatures=cell2table(morphological3dFeatures);
xls3dFeatures.Properties.VariableNames = {'ID_Cell','AspectRatio','Volume','ConvexVolume','Solidity','SurfaceArea','Sphericity'};
writetable(xls3dFeatures,fullfile(files(numFiles).folder,'3dFeatures_LimeSeg3DSegmentation.xls'), 'Range','B2');

end