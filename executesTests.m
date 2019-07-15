function tests = executesTests
    addpath(genpath('tests'));
    mode = matlab.unittest.TestCase.forInteractiveUse;
    tests = test_divideObjectInSurfaceRatios(mode);
end