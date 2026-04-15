function data = build_fig5b_data(cfg)
if nargin < 1
    cfg = struct();
end
circadian.ensure_paths();
cfg = circadian.merge_config(default_config(), cfg);

periods = cfg.fig5b.periods;
seedKd_vec = cfg.fig5b.seedKd;
lowerSeedAt_vec = cfg.fig5b.lowerSeedAt;
upperSeedAt_vec = cfg.fig5b.upperSeedAt;

assert(numel(seedKd_vec) == numel(periods), ...
    'cfg.fig5b.seedKd must align with cfg.fig5b.periods.');
assert(numel(lowerSeedAt_vec) == numel(periods), ...
    'cfg.fig5b.lowerSeedAt must align with cfg.fig5b.periods.');
assert(numel(upperSeedAt_vec) == numel(periods), ...
    'cfg.fig5b.upperSeedAt must align with cfg.fig5b.periods.');

curves = cell(1, numel(periods));
curveFiles = cell(1, numel(periods));
for i = 1:numel(periods)
    curveFiles{i} = circadian.curve_cache_file(cfg, 'fig5b', 'period', periods(i));
    curves{i} = circadian.cache_get_or_create(curveFiles{i}, @() build_iso_period_curve( ...
        cfg, periods(i), seedKd_vec(i), lowerSeedAt_vec(i), upperSeedAt_vec(i)));
end

data = struct();
data.curves = curves;
data.curveFiles = curveFiles;
end

function curve = build_iso_period_curve(cfg, targetPeriod, seedKd, lowerSeedAt, upperSeedAt)
lowerSeed = solve_seed(cfg, targetPeriod, seedKd, lowerSeedAt);
upperSeed = solve_seed(cfg, targetPeriod, seedKd, upperSeedAt);

lowerLeft = solve_kd_branch(cfg, targetPeriod, lowerSeed, cfg.fig5b.KdAxis(1));
lowerRight = solve_kd_branch(cfg, targetPeriod, lowerSeed, cfg.fig5b.KdAxis(2));
upperLeft = solve_kd_branch(cfg, targetPeriod, upperSeed, cfg.fig5b.KdAxis(1));
upperRight = solve_kd_branch(cfg, targetPeriod, upperSeed, cfg.fig5b.KdAxis(2));

curve = combine_branches(lowerLeft, lowerRight, upperLeft, upperRight, targetPeriod);
curve.targetPeriod = targetPeriod;
curve.figureId = 'fig5b';
curve.featureKind = 'period';
curve.targetValue = targetPeriod;
end

function seed = solve_seed(~, targetPeriod, seedKd, seedAt)
seed = run_circadian_fmam_task(struct( ...
    'startParams', [seedKd, seedAt], ...
    'goalOrder', {{'period', 'AT'}}, ...
    'goals', struct('period', targetPeriod, 'AT', seedAt), ...
    'controlledParams', {{'Kd', 'AT'}}, ...
    'fmam', struct('dlambdaCap', 0.05)));
end

function branch = solve_kd_branch(cfg, targetPeriod, seed, targetKd)
result = run_circadian_fmam_task(struct( ...
    'initialSolverView', seed.solverView, ...
    'initialDerivedView', seed.derivedView, ...
    'goalOrder', {{'period', 'Kd'}}, ...
    'goals', struct('period', targetPeriod, 'Kd', targetKd), ...
    'controlledParams', {{'Kd', 'AT'}}, ...
    'fmam', struct('dlambdaCap', cfg.fmam.dlambdaCap)));

branch = branch_from_path(circadian.trim_path(result.path, 'iso_period', cfg));
end

function curve = combine_branches(lowerLeft, lowerRight, upperLeft, upperRight, targetPeriod)
branches = struct( ...
    'lowerLeft', lowerLeft, ...
    'lowerRight', lowerRight, ...
    'upperLeft', upperLeft, ...
    'upperRight', upperRight);

curve = struct();
curve.targetPeriod = targetPeriod;
curve.paramNames = lowerLeft.paramNames;
curve.branches = branches;
curve.allPoints = flatten_branches(branches);
curve.paramMatrix = curve.allPoints.paramMatrix;
curve.Kd = curve.allPoints.Kd;
curve.AT = curve.allPoints.AT;
curve.period = curve.allPoints.period;
curve.obsAmplitude = curve.allPoints.obsAmplitude;
curve.obsMax = curve.allPoints.obsMax;
curve.directConditionEstimate = curve.allPoints.directConditionEstimate;
end

function branch = branch_from_path(path)
iKd = circadian.parameter_index(path, 'Kd');
iAT = circadian.parameter_index(path, 'AT');

branch = struct();
branch.paramNames = path.paramNames;
branch.paramMatrix = zeros(0, size(path.paramMatrix, 2));
branch.Kd = zeros(0, 1);
branch.AT = zeros(0, 1);
branch.period = zeros(0, 1);
branch.obsAmplitude = zeros(0, 1);
branch.obsMax = zeros(0, 1);
branch.directConditionEstimate = zeros(0, 1);
branch.stopReason = path.stopReason;

if isempty(path.paramMatrix)
    return
end

order = sortrows([(1:size(path.paramMatrix, 1)).', path.paramMatrix(:, iKd)], 2);
idx = order(:, 1);

branch.paramMatrix = path.paramMatrix(idx, :);
branch.Kd = branch.paramMatrix(:, iKd);
branch.AT = branch.paramMatrix(:, iAT);
branch.period = path.period(idx);
branch.obsAmplitude = path.obsAmplitude(idx);
branch.obsMax = path.obsMax(idx);
branch.directConditionEstimate = path.directConditionEstimate(idx);

keep = dedupe_adjacent_rows(branch.paramMatrix, 1e-10);
branch.paramMatrix = branch.paramMatrix(keep, :);
branch.Kd = branch.Kd(keep);
branch.AT = branch.AT(keep);
branch.period = branch.period(keep);
branch.obsAmplitude = branch.obsAmplitude(keep);
branch.obsMax = branch.obsMax(keep);
branch.directConditionEstimate = branch.directConditionEstimate(keep);
end

function allPoints = flatten_branches(branches)
branchNames = {'lowerLeft', 'lowerRight', 'upperLeft', 'upperRight'};
allPoints = struct();
allPoints.paramMatrix = zeros(0, size(branches.lowerLeft.paramMatrix, 2));
allPoints.Kd = zeros(0, 1);
allPoints.AT = zeros(0, 1);
allPoints.period = zeros(0, 1);
allPoints.obsAmplitude = zeros(0, 1);
allPoints.obsMax = zeros(0, 1);
allPoints.directConditionEstimate = zeros(0, 1);
allPoints.branchId = cell(0, 1);

for i = 1:numel(branchNames)
    branchName = branchNames{i};
    branch = branches.(branchName);
    n = size(branch.paramMatrix, 1);
    if n == 0
        continue
    end

    allPoints.paramMatrix = [allPoints.paramMatrix; branch.paramMatrix];
    allPoints.Kd = [allPoints.Kd; branch.Kd(:)];
    allPoints.AT = [allPoints.AT; branch.AT(:)];
    allPoints.period = [allPoints.period; branch.period(:)];
    allPoints.obsAmplitude = [allPoints.obsAmplitude; branch.obsAmplitude(:)];
    allPoints.obsMax = [allPoints.obsMax; branch.obsMax(:)];
    allPoints.directConditionEstimate = [allPoints.directConditionEstimate; branch.directConditionEstimate(:)];
    allPoints.branchId = [allPoints.branchId; repmat({branchName}, n, 1)];
end
end

function keep = dedupe_adjacent_rows(matrix, tol)
n = size(matrix, 1);
keep = true(n, 1);
if n < 2
    return
end

for i = 2:n
    if norm(matrix(i, :) - matrix(i - 1, :), inf) < tol
        keep(i) = false;
    end
end
end
