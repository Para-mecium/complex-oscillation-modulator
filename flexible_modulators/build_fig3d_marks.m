function data = build_fig3d_marks(cfg)
if nargin < 1
    cfg = struct();
end
flexmod.ensure_paths();
cfg = flexmod.merge_config(default_config(), cfg);

markPeriods = cfg.fig3d.markPeriods;
markSeeds = cfg.fig3d.markSeeds;
assert(size(markSeeds, 1) == numel(markPeriods) && size(markSeeds, 2) == 2, ...
    'cfg.fig3d.markSeeds must be an N x 2 matrix aligned with cfg.fig3d.markPeriods.');

markResults = cell(1, numel(markPeriods));
markFiles = cell(1, numel(markPeriods));
for i = 1:numel(markPeriods)
    markFiles{i} = flexmod.mark_cache_file(cfg, 'fig3d', sprintf('period_%02d', i));
    seedParams = markSeeds(i, :);
    targetPeriod = markPeriods(i);
    markResults{i} = flexmod.cache_get_or_create(markFiles{i}, ...
        @() solve_mark_point(cfg, cfg.fig3d.markAmplitude, targetPeriod, seedParams, i));
end

data = struct();
data.markResults = markResults;
data.markFiles = markFiles;
data.markPeriods = markPeriods;
end

function mark = solve_mark_point(cfg, targetAmplitude, targetPeriod, startParams, markIndex)
mark = run_modulation_task(struct( ...
    'modelType', 'base', ...
    'startParams', startParams, ...
    'goalOrder', {{'period', 'amplitude'}}, ...
    'goals', struct('period', targetPeriod, 'amplitude', targetAmplitude), ...
    'controlledParams', {{'I1', 'ET'}}, ...
    'fmam', struct('dlambdaCap', 1)));
mark.orbit = flexmod.shift_cycle_to_max(mark.orbit);
mark.figureId = 'fig3d';
mark.markIndex = markIndex;
mark.markType = 'period';
mark.targetPeriod = targetPeriod;
mark.targetAmplitude = targetAmplitude;
mark.seedParams = reshape(startParams, 1, []);
end
