function bifurcation = base_bifurcation_curve(cfg)
if nargin < 1
    cfg = struct();
end

flexmod.ensure_paths();
cfg = flexmod.merge_config(default_config(), cfg);

xRange = [min(cfg.grid.I1Values), max(cfg.grid.I1Values)];
yRange = [min(cfg.grid.ETValues), max(cfg.grid.ETValues)];
raw = load_or_build_raw_bifurcation(cfg);

hopfCurves = clip_curve_set(raw.hopfCurves, xRange, yRange);
lpCurves = clip_curve_set(raw.lpCurves, xRange, yRange);

lowerHopf = select_primary_component(build_boundary_envelope(hopfCurves, xRange, yRange, 'lower'));
lowerLp = select_primary_component(build_boundary_envelope(lpCurves, xRange, yRange, 'lower'));
upperBoundary = lowerLp;
allCurves = join_curve_sets(hopfCurves, lpCurves);
visibleHopfCurves = {};
visibleLpCurves = {};
if ~isempty(lowerHopf.I1)
    visibleHopfCurves = {lowerHopf};
end
if ~isempty(lowerLp.I1)
    visibleLpCurves = {lowerLp};
end

[lowerRegionX, lowerRegionY] = build_fill_region(lowerHopf, yRange, 'bottom');
[upperRegionX, upperRegionY] = build_fill_region(upperBoundary, yRange, 'top');

bifurcation = struct();
bifurcation.raw = raw;
bifurcation.hopfCurves = hopfCurves;
bifurcation.lpCurves = lpCurves;
bifurcation.visibleHopfCurves = visibleHopfCurves;
bifurcation.visibleLpCurves = visibleLpCurves;
bifurcation.allCurves = allCurves;
bifurcation.lowerCurve = lowerHopf;
bifurcation.upperCurve = upperBoundary;
bifurcation.lowerLpCurve = lowerLp;
bifurcation.lowerRegionX = lowerRegionX;
bifurcation.lowerRegionY = lowerRegionY;
bifurcation.upperRegionX = upperRegionX;
bifurcation.upperRegionY = upperRegionY;
bifurcation.xRange = xRange;
bifurcation.yRange = yRange;
end

function raw = load_or_build_raw_bifurcation(cfg)
cacheFile = fullfile(cfg.io.cacheDir, 'base_bifurcation_matcont_v1.mat');
if isfile(cacheFile)
    loaded = load(cacheFile);
    if isfield(loaded, 'data')
        raw = loaded.data;
        return
    end
end

raw = build_raw_bifurcation(cfg);

cacheDir = fileparts(cacheFile);
if ~isfolder(cacheDir)
    mkdir(cacheDir);
end
data = raw;
save(cacheFile, 'data', '-v7.3');
end

function raw = build_raw_bifurcation(cfg)
c = cfg.base.constants;
matcont_base_system('set_constants', c);

seedParams = cfg.base.initialParams(:);
seedState = find_equilibrium(seedParams, c);

activeParam1D = 2;
[epX, epS] = continue_equilibrium(seedState, seedParams, activeParam1D);

hopfSeeds = collect_singularity_seeds(epX, epS, 'H', seedParams, 2);
lpSeeds = collect_singularity_seeds(epX, epS, 'LP', seedParams, 2);

hopfCurves = continue_curve_family(@init_H_H, @hopf, hopfSeeds, [1 2], 5);
lpCurves = continue_curve_family(@init_LP_LP, @limitpoint, lpSeeds, [1 2], 4);

raw = struct();
raw.seedState = seedState;
raw.seedParams = seedParams;
raw.equilibrium = struct('x', epX, 's', epS);
raw.hopfCurves = hopfCurves;
raw.lpCurves = lpCurves;
end

function state = find_equilibrium(params, c)
I1 = params(1);
ET = params(2);
drive = regulated_drive(I1, c);
yCandidates = linspace(1e-6, 4, 4000);
values = equilibrium_residual(yCandidates, ET, drive, c);
signChange = find(values(1:end-1) .* values(2:end) <= 0, 1, 'first');

if isempty(signChange)
    [~, idx] = min(abs(values));
    yEq = yCandidates(idx);
else
    left = yCandidates(signChange);
    right = yCandidates(signChange + 1);
    yEq = fzero(@(y) equilibrium_residual(y, ET, drive, c), [left, right]);
end

denominator = c.Km + yEq + c.KI * yEq^2;
xEq = (c.kdy * yEq + c.k2 * ET * yEq / denominator) / c.ksy;
state = [xEq; yEq];
end

function value = equilibrium_residual(y, ET, drive, c)
production = c.k1 * c.S * c.Kd^c.p ./ (c.Kd^c.p + y.^c.p);
denominator = c.Km + y + c.KI * y.^2;
value = c.ksy * production * drive / c.kdx - c.kdy * y - c.k2 * ET * y ./ denominator;
end

function [x, s] = continue_equilibrium(seedState, seedParams, activeParam)
opt = contset;
opt = contset(opt, 'Singularities', 1);
opt = contset(opt, 'InitStepsize', 1e-3);
opt = contset(opt, 'MinStepsize', 1e-6);
opt = contset(opt, 'MaxStepsize', 0.05);
opt = contset(opt, 'MaxNumPoints', 400);

[x0, v0] = init_EP_EP(@matcont_base_system, seedState, seedParams, activeParam);
[forwardX, ~, forwardS] = cont(@equilibrium, x0, v0, opt);

opt = contset(opt, 'Backward', 1);
[backwardX, ~, backwardS] = cont(@equilibrium, x0, v0, opt);

[x, s] = merge_continuation_runs( ...
    extract_equilibrium_curve(forwardX, forwardS, activeParam), ...
    extract_equilibrium_curve(backwardX, backwardS, activeParam), ...
    activeParam);
end

function curve = extract_equilibrium_curve(x, s, activeParam)
curve = struct();
curve.points = x;
curve.s = reindex_singularities(s, size(x, 2));
curve.activeParam = activeParam;
end

function seeds = collect_singularity_seeds(epX, epS, label, baseParams, activeParam)
labels = arrayfun(@(item) strtrim(item.label), epS, 'UniformOutput', false);
matching = strcmp(labels, label);
selected = epS(matching);
seeds = repmat(struct('state', [], 'params', []), 1, numel(selected));
for i = 1:numel(selected)
    index = selected(i).index;
    params = baseParams;
    params(activeParam) = epX(end, index);
    seeds(i).state = epX(1:2, index);
    seeds(i).params = params;
end
end

function curves = continue_curve_family(initFunction, curveFunction, seeds, activeParams, minColumns)
curves = {};
for i = 1:numel(seeds)
    seed = seeds(i);
    directions = {'forward', 'backward'};
    for j = 1:numel(directions)
        direction = directions{j};
        opt = contset;
        opt = contset(opt, 'Singularities', 1);
        opt = contset(opt, 'InitStepsize', 1e-3);
        opt = contset(opt, 'MinStepsize', 1e-6);
        opt = contset(opt, 'MaxStepsize', 0.03);
        opt = contset(opt, 'MaxNumPoints', 600);
        opt = contset(opt, 'FunTolerance', 1e-8);
        opt = contset(opt, 'VarTolerance', 1e-8);
        if strcmp(direction, 'backward')
            opt = contset(opt, 'Backward', 1);
        end

        try
            [x0, v0] = initFunction(@matcont_base_system, seed.state, seed.params, activeParams);
            [x, ~] = cont(curveFunction, x0, v0, opt);
        catch
            continue
        end

        if size(x, 2) < minColumns
            continue
        end

        curve = struct();
        curve.I1 = x(3, :).';
        curve.ET = x(4, :).';
        curve.state = x(1:2, :).';
        curve.seedIndex = i;
        curve.direction = direction;
        curves{end + 1} = sanitize_curve(curve); %#ok<AGROW>
    end
end
end

function [mergedX, mergedS] = merge_continuation_runs(forward, backward, ~)
forwardPoints = forward.points;
backwardPoints = fliplr(backward.points);
nBackwardOriginal = size(backward.points, 2);
backwardS = reindex_singularities(backward.s, nBackwardOriginal);
for i = 1:numel(backwardS)
    backwardS(i).index = nBackwardOriginal - backwardS(i).index + 1;
end

if ~isempty(backwardPoints) && ~isempty(forwardPoints)
    duplicate = max(abs(backwardPoints(:, end) - forwardPoints(:, 1))) < 1e-9;
    if duplicate
        backwardPoints(:, end) = [];
        keep = [backwardS.index] <= size(backwardPoints, 2);
        backwardS = backwardS(keep);
    end
end

mergedX = [backwardPoints, forwardPoints];
forwardS = reindex_singularities(forward.s, size(forward.points, 2));
if ~isempty(backwardPoints)
    offset = size(backwardPoints, 2);
    for i = 1:numel(forwardS)
        forwardS(i).index = forwardS(i).index + offset;
    end
end
mergedS = concat_singularity_arrays(backwardS, forwardS);

if isempty(mergedX)
    mergedX = zeros(3, 0);
end
if isempty(mergedS)
    mergedS = struct('label', {}, 'index', {});
end
end

function singularities = reindex_singularities(singularities, count)
if isempty(singularities)
    if isstruct(singularities)
        singularities = reshape(singularities, 1, 0);
    else
        singularities = struct('label', {}, 'index', {});
    end
    return
end

keep = [singularities.index] >= 1 & [singularities.index] <= count;
singularities = singularities(keep);
singularities = reshape(singularities, 1, []);
end

function merged = concat_singularity_arrays(first, second)
first = normalize_singularity_array(first);
second = normalize_singularity_array(second);

if isempty(first)
    merged = second;
    return
end
if isempty(second)
    merged = first;
    return
end

fieldOrder = union(fieldnames(first), fieldnames(second), 'stable');
first = align_singularity_fields(first, fieldOrder);
second = align_singularity_fields(second, fieldOrder);
merged = [first, second];
end

function singularities = normalize_singularity_array(singularities)
if isempty(singularities)
    if isstruct(singularities)
        singularities = reshape(singularities, 1, 0);
    end
    return
end
singularities = reshape(singularities, 1, []);
end

function singularities = align_singularity_fields(singularities, fieldOrder)
for i = 1:numel(fieldOrder)
    fieldName = fieldOrder{i};
    if ~isfield(singularities, fieldName)
        [singularities.(fieldName)] = deal([]);
    end
end
singularities = orderfields(singularities, fieldOrder);
end

function normalized = sanitize_curve(curve)
normalized = curve;
if isempty(curve.I1)
    return
end

keep = [true; abs(diff(curve.I1)) > 1e-10 | abs(diff(curve.ET)) > 1e-10];
normalized.I1 = curve.I1(keep);
normalized.ET = curve.ET(keep);
if isfield(curve, 'state')
    normalized.state = curve.state(keep, :);
end
end

function clippedCurves = clip_curve_set(curves, xRange, yRange)
clippedCurves = {};
for i = 1:numel(curves)
    curve = curves{i};
    visible = curve.I1 >= xRange(1) & curve.I1 <= xRange(2) & ...
        curve.ET >= yRange(1) & curve.ET <= yRange(2);
    runs = contiguous_true_runs(visible);
    for j = 1:numel(runs)
        idx = runs{j};
        if numel(idx) < 2
            continue
        end

        clipped = curve;
        clipped.I1 = curve.I1(idx);
        clipped.ET = curve.ET(idx);
        if isfield(curve, 'state')
            clipped.state = curve.state(idx, :);
        end
        clippedCurves{end + 1} = sanitize_curve(clipped); %#ok<AGROW>
    end
end
end

function runs = contiguous_true_runs(mask)
runs = {};
if ~any(mask)
    return
end

indices = find(mask);
splitPoints = [0; find(diff(indices) > 1); numel(indices)];
for i = 1:(numel(splitPoints) - 1)
    runs{end + 1} = indices((splitPoints(i) + 1):splitPoints(i + 1)); %#ok<AGROW>
end
end

function combined = join_curve_sets(firstSet, secondSet)
combined = firstSet;
for i = 1:numel(secondSet)
    combined{end + 1} = secondSet{i}; %#ok<AGROW>
end
end

function boundary = build_boundary_envelope(curves, xRange, yRange, mode)
boundary = empty_curve();
if isempty(curves)
    return
end

xGrid = linspace(xRange(1), xRange(2), 600).';
boundaryY = NaN(size(xGrid));
for i = 1:numel(curves)
    candidate = interpolate_curve(curves{i}, xGrid);
    boundaryY = combine_boundary_values(boundaryY, candidate, mode);
end

valid = isfinite(boundaryY);
boundary.I1 = xGrid(valid);
boundary.ET = min(max(boundaryY(valid), yRange(1)), yRange(2));
end

function curve = select_primary_component(curveRaw)
curve = curveRaw;
if numel(curveRaw.I1) < 2
    return
end

dx = median(diff(curveRaw.I1));
if ~(isfinite(dx) && dx > 0)
    return
end

gapMask = diff(curveRaw.I1) > 1.5 * dx;
if ~any(gapMask)
    return
end

breaks = [0; find(gapMask); numel(curveRaw.I1)];
bestIdx = [];
bestSpan = -Inf;
bestLevel = Inf;
for i = 1:(numel(breaks) - 1)
    idx = (breaks(i) + 1):breaks(i + 1);
    if numel(idx) < 2
        continue
    end

    span = curveRaw.I1(idx(end)) - curveRaw.I1(idx(1));
    level = median(curveRaw.ET(idx));
    if span > bestSpan + 1e-12 || (abs(span - bestSpan) <= 1e-12 && level < bestLevel)
        bestIdx = idx;
        bestSpan = span;
        bestLevel = level;
    end
end

if isempty(bestIdx)
    curve = empty_curve();
    return
end

curve.I1 = curveRaw.I1(bestIdx);
curve.ET = curveRaw.ET(bestIdx);
end

function values = interpolate_curve(curve, xGrid)
values = NaN(size(xGrid));
if numel(curve.I1) < 2
    return
end

segments = split_monotone_segments(curve.I1);
for i = 1:numel(segments)
    idx = segments{i};
    xSegment = curve.I1(idx);
    ySegment = curve.ET(idx);
    [xUnique, uniqueIdx] = unique(xSegment, 'stable');
    yUnique = ySegment(uniqueIdx);
    if numel(xUnique) < 2
        continue
    end

    if xUnique(1) > xUnique(end)
        xUnique = flipud(xUnique);
        yUnique = flipud(yUnique);
    end

    inside = xGrid >= xUnique(1) & xGrid <= xUnique(end);
    candidate = NaN(size(xGrid));
    candidate(inside) = interp1(xUnique, yUnique, xGrid(inside), 'pchip');
    values = combine_boundary_values(values, candidate, 'lower');
end
end

function segments = split_monotone_segments(x)
segments = {};
if numel(x) < 2
    segments{1} = 1:numel(x);
    return
end

dx = diff(x(:));
startIdx = 1;
lastSign = 0;
for i = 1:numel(dx)
    currentSign = sign(dx(i));
    if currentSign == 0
        continue
    end
    if lastSign == 0
        lastSign = currentSign;
        continue
    end
    if currentSign ~= lastSign
        segments{end + 1} = startIdx:i; %#ok<AGROW>
        startIdx = i;
        lastSign = currentSign;
    end
end
segments{end + 1} = startIdx:numel(x);
end

function merged = combine_boundary_values(base, candidate, mode)
merged = base;
candidateOnly = isnan(base) & isfinite(candidate);
merged(candidateOnly) = candidate(candidateOnly);

both = isfinite(base) & isfinite(candidate);
switch lower(mode)
    case 'lower'
        merged(both) = min(base(both), candidate(both));
    case 'upper'
        merged(both) = max(base(both), candidate(both));
    otherwise
        error('flexmod:UnknownBoundaryMode', 'Unsupported boundary mode ''%s''.', mode);
end
end

function [regionX, regionY] = build_fill_region(curve, yRange, sideName)
if isempty(curve.I1)
    regionX = zeros(0, 1);
    regionY = zeros(0, 1);
    return
end

curveX = curve.I1(:);
curveY = min(max(curve.ET(:), yRange(1)), yRange(2));
switch lower(sideName)
    case 'bottom'
        boundaryY = yRange(1) * ones(size(curveY));
    case 'top'
        boundaryY = yRange(2) * ones(size(curveY));
    otherwise
        error('flexmod:UnknownFillSide', 'Unsupported fill side ''%s''.', sideName);
end

regionX = [curveX; flipud(curveX)];
regionY = [curveY; flipud(boundaryY)];
end

function curve = empty_curve()
curve = struct('I1', zeros(0, 1), 'ET', zeros(0, 1));
end

function value = regulated_drive(I1, c)
hU = c.U * I1^2 / (c.K1 + I1^2);
value = c.bU * hU^2 / (c.KU + hU^2);
end
