function test_mapped3Dinto2D(testCase, selpath, surfaceName)
%TEST_MAPPED3DINTO2D_APICAL Summary of this function goes here
%   Detailed explanation goes here
    warning('off','all')
    outputDir = fullfile(selpath, strcat('unrolledGlands\gland_SR_', surfaceName,'\'));
    
    
    %% Actual solution
    load(fullfile(outputDir, 'final3DImg.mat'));
    
    if sum(size(img3d) > 1000) >= 2 %Real size image
        closingPxAreas2D = 1;
        closingPxAreas3D = 0;
    else %Reduced image
        closingPxAreas3D = 10;
        closingPxAreas2D = closingPxAreas3D;
    end
    
    load(fullfile(selpath, 'valid_cells.mat'));
    load(fullfile(selpath, '3d_layers_info.mat'), 'colours');
    
    [cylindre2DImage, newVerticesNeighs2D, newVertices2D, centroids, ...
        validCellsFinal, borderCells, surfaceRatio, nameOfSimulation, ...
        areaOfValidCells, deployedImg, deployedImg3x, wholeImage, ...
        polygon_distribution, newNeighbours2D, newNeighbours2D_Checked] = mappCylindricalCoordinatesInto2D(img3d, img3dComplete, closingPxAreas2D, noValidCells, colours);
    

    actSolution = {cylindre2DImage, newVerticesNeighs2D, newVertices2D, ...
        validCellsFinal, borderCells, surfaceRatio, nameOfSimulation, ...
        areaOfValidCells, deployedImg, deployedImg3x, wholeImage, ...
        polygon_distribution, newNeighbours2D, newNeighbours2D_Checked};
    
    actualNeighbours = {newVerticesNeighs2D, newVertices2D};
    
    %% Expected solution
    load(fullfile(outputDir, 'allInfo.mat'));
    load(fullfile(outputDir, 'verticesInfo.mat'));
    
    expSolution = {cylindre2DImage, newVerticesNeighs2D, newVertices2D, ...
        validCellsFinal, borderCells, surfaceRatio, nameOfSimulation, ...
        areaOfValidCells, deployedImg, deployedImg3x, wholeImage, ...
        polygon_distribution, newNeighbours2D, newNeighbours2D_Checked};
    
    verifyEqual(testCase,actSolution,expSolution);
    
    disp('%----- Checking neighbours -----%')
    expectedNeighbours = {newVerticesNeighs2D, newVertices2D};
    verifyEqual(testCase, actualNeighbours, expectedNeighbours);

    warning('on','all')
end

