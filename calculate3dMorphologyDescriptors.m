close all
clear all

files = dir('**/data/Salivary gland/**/Results/3d_layers_info.mat');

for numFiles=1:length(files)
load(fullfile(files(numFiles).folder, '3d_layers_info.mat'), 'labelledImage'); 
load(fullfile(files(numFiles).folder, 'valid_cells.mat'), 'noValidCells'); 

morphology3dFeatures={};
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
        
        morphology3dFeatures{numCell,1} = numCell;
        morphology3dFeatures{numCell,2} = aspectRatio;
        morphology3dFeatures{numCell,3} = volume;
        morphology3dFeatures{numCell,4} = convexVolume;
        morphology3dFeatures{numCell,5} = solidity;
        morphology3dFeatures{numCell,6} = surfaceArea;
        morphology3dFeatures{numCell,7} = sphericity;

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

morphology3dFeatures(noValidCells,:)=[];
totalCellPrincipleAxis(noValidCells,:)=[];

%% Save variables
totalAspectRatio = {morphology3dFeatures{:,2}};
totalVolume = {morphology3dFeatures{:,3}};
totalConvexVolume = {morphology3dFeatures{:,4}};
totalSolidity = {morphology3dFeatures{:,5}};
totalSurfaceArea = {morphology3dFeatures{:,6}};
totalSphericity = {morphology3dFeatures{:,7}};

save(fullfile(files(numFiles).folder, 'morphology3dFeatures.mat'), 'totalAspectRatio','totalVolume','totalConvexVolume','totalSolidity','totalSurfaceArea','totalSphericity','totalCellPrincipleAxis'); 

%% Export to excel
xls3dFeatures=cell2table(morphology3dFeatures);
xls3dFeatures.Properties.VariableNames = {'ID_Cell','AspectRatio','Volume','ConvexVolume','Solidity','SurfaceArea','Sphericity'};
writetable(xls3dFeatures,fullfile(files(numFiles).folder,'3dFeatures_LimeSeg3DSegmentation.xls'), 'Range','B2');

end