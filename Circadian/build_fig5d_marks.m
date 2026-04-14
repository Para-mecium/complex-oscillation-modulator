function data = build_fig5d_marks(cfg)
if nargin < 1
    cfg = struct();
end
circadian.ensure_paths();
cfg = circadian.merge_config(default_config(), cfg);

markPeriods = cfg.fig5d.markPeriods;
markSeeds = cfg.fig5d.markSeeds;
assert(size(markSeeds, 1) == numel(markPeriods) && size(markSeeds, 2) == 2, ...
    'cfg.fig5d.markSeeds must be an N x 2 matrix aligned with cfg.fig5d.markPeriods.');

markResults = cell(1, numel(markPeriods));
markFiles = cell(1, numel(markPeriods));
for i = 1:numel(markPeriods)
    markFiles{i} = circadian.mark_cache_file(cfg, 'fig5d', sprintf('period_%02d', i));
    targetPeriod = markPeriods(i);
    seedParams = markSeeds(i, :);
    markResults{i} = circadian.cache_get_or_create(markFiles{i}, ...
        @() solve_mark_point(cfg, cfg.fig5d.markMaximum, targetPeriod, seedParams, i));
end

data = struct();
data.markResults = markResults;
data.markFiles = markFiles;
data.markPeriods = markPeriods;
end

function mark = solve_mark_point(cfg, targetMaximum, targetPeriod, startParams, markIndex)
mark = run_modulation_task(struct( ...
    'startParams', startParams, ...
    'goalOrder', {{'maximum', 'period'}}, ...
    'goals', struct('maximum', targetMaximum, 'period', targetPeriod), ...
    'controlledParams', {{'Kd', 'AT'}}, ...
    'fmam', struct('dlambdaCap', 1)));
mark.orbit = circadian.shift_cycle_to_max(mark.orbit);
mark.figureId = 'fig5d';
mark.markIndex = markIndex;
mark.markType = 'period';
mark.targetMaximum = targetMaximum;
mark.targetPeriod = targetPeriod;
mark.seedParams = reshape(startParams, 1, []);
end
