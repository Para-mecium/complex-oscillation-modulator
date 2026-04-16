function f2Curve = solve_curve_continuation(f1Grid, residualFn, seed, opts)
%SOLVE_CURVE_CONTINUATION Solve residual(f1,f2)=0 along a prescribed f1 grid.

if nargin < 4
    opts = struct();
end

cfg = normalform.merge_options(default_options(), opts);

f1Grid = f1Grid(:);
f2Curve = nan(size(f1Grid));

if isempty(seed)
    seed = 0;
end

if cfg.use_fsolve && exist('fsolve', 'file') == 2
    f2Curve = solve_with_fsolve(f1Grid, residualFn, seed);
    return
end

f2Curve = solve_with_bracket(f1Grid, residualFn, seed, cfg.search_interval);
end

function cfg = default_options()
cfg = struct();
cfg.use_fsolve = true;
cfg.search_interval = [-4, 4];
end

function f2Curve = solve_with_fsolve(f1Grid, residualFn, seed)
f2Curve = nan(size(f1Grid));
solverOptions = optimoptions('fsolve', 'Display', 'off');

state = seed;
positiveIdx = find(f1Grid >= 0);
negativeIdx = flipud(find(f1Grid < 0));

for idx = positiveIdx(:).'
    state = fsolve(@(f2) residualFn(f1Grid(idx), f2), state, solverOptions);
    f2Curve(idx) = state;
end

for idx = negativeIdx(:).'
    state = fsolve(@(f2) residualFn(f1Grid(idx), f2), state, solverOptions);
    f2Curve(idx) = state;
end
end

function f2Curve = solve_with_bracket(f1Grid, residualFn, seed, searchInterval)
f2Curve = nan(size(f1Grid));
anchorIdx = find(abs(f1Grid) == min(abs(f1Grid)), 1, 'first');
f2Curve(anchorIdx) = solve_root_for_point(f1Grid(anchorIdx), residualFn, seed, searchInterval);

for idx = (anchorIdx + 1):numel(f1Grid)
    f2Curve(idx) = solve_root_for_point(f1Grid(idx), residualFn, f2Curve(idx - 1), searchInterval);
end

for idx = (anchorIdx - 1):-1:1
    f2Curve(idx) = solve_root_for_point(f1Grid(idx), residualFn, f2Curve(idx + 1), searchInterval);
end
end

function root = solve_root_for_point(f1, residualFn, guess, searchInterval)
residualAtGuess = residualFn(f1, guess);
if isfinite(residualAtGuess) && abs(residualAtGuess) < 1e-12
    root = guess;
    return
end

[left, right] = find_bracket(f1, residualFn, guess, searchInterval);
if left == right
    root = left;
    return
end

root = fzero(@(f2) residualFn(f1, f2), [left, right]);
end

function [left, right] = find_bracket(f1, residualFn, guess, searchInterval)
span = 0.15;

for iter = 1:8
    grid = linspace(max(searchInterval(1), guess - span), ...
        min(searchInterval(2), guess + span), 200);
    [left, right, found] = bracket_from_grid(f1, residualFn, grid, guess);
    if found
        return
    end
    span = span * 2;
end

grid = linspace(searchInterval(1), searchInterval(2), 2000);
[left, right, found] = bracket_from_grid(f1, residualFn, grid, guess);
if ~found
    error('normalform:NoBracket', ...
        'Failed to bracket a root for continuation point f1=%g.', f1);
end
end

function [left, right, found] = bracket_from_grid(f1, residualFn, grid, guess)
values = nan(size(grid));
for idx = 1:numel(grid)
    values(idx) = residualFn(f1, grid(idx));
end

finiteMask = isfinite(values);
grid = grid(finiteMask);
values = values(finiteMask);

left = nan;
right = nan;
found = false;
candidates = zeros(0, 2);

for idx = 1:(numel(grid) - 1)
    if values(idx) == 0
        candidates(end + 1, :) = [grid(idx), grid(idx)]; %#ok<AGROW>
    end
    if values(idx) * values(idx + 1) < 0
        candidates(end + 1, :) = [grid(idx), grid(idx + 1)]; %#ok<AGROW>
    end
end

if isempty(candidates)
    return
end

midpoints = mean(candidates, 2);
[~, bestIdx] = min(abs(midpoints - guess));
left = candidates(bestIdx, 1);
right = candidates(bestIdx, 2);
found = true;
end
