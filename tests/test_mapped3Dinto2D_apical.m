function test_mapped3Dinto2D_apical(testCase)
%TEST_MAPPED3DINTO2D_APICAL Summary of this function goes here
%   Detailed explanation goes here
    warning('off','all')
    outputDir = 'D:\Pablo\LimeSeg_Pipeline\data\Salivary gland_ExtractedVertices_Correct\Wildtype\2017-12-04\1a\Results\unrolledGlands\gland_SR_1\';
    
    %% Actual solution
    load(fullfile(outputDir, 'verticesInfo.mat'));
    [cylindre2DImage, newVerticesNeighs2D, newVertices2D, centroids, ...
        validCellsFinal, borderCells, surfaceRatio, nameOfSimulation, ...
        areaOfValidCells, deployedImg, deployedImg3x, wholeImage, ...
        polygon_distribution, newNeighbours2D, newNeighbours2D_Checked] = mappCylindricalCoordinatesInto2D(img3d, img3dComplete, testing);
    

    actSolution = {cylindre2DImage, newVerticesNeighs2D, newVertices2D, centroids, ...
        validCellsFinal, borderCells, surfaceRatio, nameOfSimulation, ...
        areaOfValidCells, deployedImg, deployedImg3x, wholeImage, ...
        polygon_distribution, newNeighbours2D, newNeighbours2D_Checked};
    
    
    %% Expected solution
    load(fullfile(outputDir, 'allInfo.mat'));
    load(fullfile(outputDir, 'verticesInfo.mat'));
    
    expSolution = {cylindre2DImage, newVerticesNeighs2D, newVertices2D, centroids, ...
        validCellsFinal, borderCells, surfaceRatio, nameOfSimulation, ...
        areaOfValidCells, deployedImg, deployedImg3x, wholeImage, ...
        polygon_distribution, newNeighbours2D, newNeighbours2D_Checked};
    
    verifyEqual(testCase,actSolution,expSolution);

    warning('on','all')
end

