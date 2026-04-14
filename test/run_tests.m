function results = run_tests()
%RUN_TESTS Run the FMAM_code MATLAB regression tests.

    suite = matlab.unittest.TestSuite.fromFolder(fileparts(mfilename('fullpath')), ...
        'IncludingSubfolders', false);
    runner = matlab.unittest.TestRunner.withTextOutput('Verbosity', 2);
    results = runner.run(suite);
end
