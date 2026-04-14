function tests = test_legacy_boundary_contracts
tests = functiontests(localfunctions);
end

function testFmamStateOpsDoesNotReferenceUtilsImpl(testCase)
source = fileread(fullfile(repo_root(), 'fmam_state_ops.m'));
verifyFalse(testCase, contains(source, 'utils.state_ops'));
end

function testUtilsStateOpsIsShimToFmamStateOps(testCase)
source = fileread(fullfile(repo_root(), '+utils', 'state_ops.m'));
verifyTrue(testCase, contains(source, 'fmam_state_ops.'));
end

function testProductionCodeDoesNotReadStateDefaultConstantsOutsideCompatLayer(testCase)
files = { ...
    fullfile(repo_root(), 'FMAM_ODE.m'), ...
    fullfile(repo_root(), 'assemble_newton_linear_system.m'), ...
    fullfile(repo_root(), 'fmam_state_ops.m'), ...
    fullfile(repo_root(), 'Circadian', '+circadian', 'run_fmam_task.m'), ...
    fullfile(repo_root(), 'flexible_modulators', '+flexmod', 'run_fmam_task.m')};

patterns = {'state.Lconst', 'state.LphiConst', 'state.countMax', 'state.errMax'};
for i = 1:numel(files)
    source = fileread(files{i});
    for j = 1:numel(patterns)
        verifyFalse(testCase, contains(source, patterns{j}), ...
            sprintf('Unexpected legacy default constant reference in %s: %s', files{i}, patterns{j}));
    end
end
end

function root = repo_root()
root = fileparts(fileparts(mfilename('fullpath')));
end
