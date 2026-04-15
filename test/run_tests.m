function results = run_tests()
%RUN_TESTS Run the FMAM_code MATLAB regression tests.

    testDir = fileparts(mfilename('fullpath'));
    rootDir = fileparts(testDir);
    addpath(rootDir);

    suite = matlab.unittest.TestSuite.fromFolder(testDir, ...
        'IncludingSubfolders', false);
    runner = matlab.unittest.TestRunner.withTextOutput('Verbosity', 2);
    results = runner.run(suite);
end
