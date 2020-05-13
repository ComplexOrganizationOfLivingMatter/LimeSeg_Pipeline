selectedPath= 'E:\Pedro\LimeSeg_Pipeline\data\E-cadh Inhibited (flatten)\2019-02-07\3a1539\';
unrollPath = [selectedPath 'Results\unrolledGlands\'];
selectedSR={'gland_SR_1'};

%% CorrectIndividual SR unrolled.
% %1st - 
% %remove -> verticesInfo.mat, allInfo.mat, samiraTable.mat and
% %polygon_distribution.xls
% for nSr = 1:length(selectedSR)
%     delete([unrollPath selectedSR{nSr} '\samiraTable.mat'])
%     delete([unrollPath selectedSR{nSr} '\verticesInfo.mat'])
%     delete([unrollPath selectedSR{nSr} '\allInfo.mat'])
%     delete([unrollPath selectedSR{nSr} '\polygon_distribution.xls'])
% end

%2nd - 
%run mappCylindricalCoordinatesInto2D(img3d, img3dComplete, closingPxAreas2D, noValidCells, colours, outputDir)
%run connectVerticesOf2D(deployedImg, newVerticesNeighs2D, newVertices2D, centroids, validCells, borderCells, surfaceRatio, outputDir, nameOfSimulation, deployedImg3x, img3d);
load([selectedPath 'Results\3d_layers_info.mat'],'colours')
load([selectedPath 'Results\valid_cells.mat'],'validCells','noValidCells')
closingPxAreas2D = 10;
colours = [1,1,1; colours];
for nSr = 1:length(selectedSR)
    outputDirActual =[unrollPath selectedSR{nSr} '\'];
    load([outputDirActual 'final3DImg.mat'],'img3d','img3dComplete')
%     mappCylindricalCoordinatesInto2D(img3d, img3dComplete, closingPxAreas2D, noValidCells, colours, outputDir)
    load(fullfile(outputDirActual, 'allInfo.mat'));
    load(fullfile(outputDirActual, 'verticesInfo.mat'));
    samiraTable = connectVerticesOf2D(deployedImg, newVerticesNeighs2D, newVertices2D, centroids, validCells, borderCells, surfaceRatio, outputDirActual, nameOfSimulation, deployedImg3x, img3d);
    save([outputDirActual 'samiraTable.mat'], 'samiraTable');
end

%3rd -
%remove -> _samirasFormat.xls ; _VertCrosses.xls and
%glandDividedInSurfacesRatios_AllUnrollFeatures
splittedPath = strsplit(selectedPath,'\');
delete([selectedPath 'Results\' splittedPath{end-3} '_' splittedPath{end-2} '_' splittedPath{end-1} '_samirasFormat.xls']);
delete([selectedPath 'Results\' splittedPath{end-3} '_' splittedPath{end-2} '_' splittedPath{end-1} '_VertCrosses.xls']);
delete([selectedPath 'Results\glandDividedInSurfacesRatios_AllUnrollFeatures.mat']);

%once everything corrected
unrollTube_parallel([selectedPath 'Results\'])
