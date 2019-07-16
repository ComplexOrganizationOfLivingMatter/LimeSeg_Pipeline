function tests = executesTests
    addpath(genpath('tests'));
    mode = matlab.unittest.TestCase.forInteractiveUse;
    
    perfectSample = '';
    
    %Testing function 'divideObjectInSurfaceRatios'
    %test_divideObjectInSurfaceRatios(mode);
    
    %Testing function 'mappCylindricalCoordinatesInto2D' from 'unrollTube'
    test_mapped3Dinto2D_apical(mode)
end