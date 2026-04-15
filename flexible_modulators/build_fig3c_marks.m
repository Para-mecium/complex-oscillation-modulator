function data = build_fig3c_marks(cfg)
if nargin < 1
    cfg = struct();
end
flexmod.ensure_paths();
cfg = flexmod.merge_config(default_config(), cfg);

markAmplitudes = cfg.fig3c.markAmplitudes;
markSeeds = cfg.fig3c.markSeeds;
assert(size(markSeeds, 1) == numel(markAmplitudes) && size(markSeeds, 2) == 2, ...
    'cfg.fig3c.markSeeds must be an N x 2 matrix aligned with cfg.fig3c.markAmplitudes.');

markResults = cell(1, numel(markAmplitudes));
markFiles = cell(1, numel(markAmplitudes));
for i = 1:numel(markAmplitudes)
    markFiles{i} = flexmod.mark_cache_file(cfg, 'fig3c', sprintf('amplitude_%02d', i));
    seedParams = markSeeds(i, :);
    targetAmplitude = markAmplitudes(i);
    markResults{i} = flexmod.cache_get_or_create(markFiles{i}, ...
        @() solve_mark_point(cfg, cfg.fig3c.markPeriod, targetAmplitude, seedParams, i));
end

data = struct();
data.markResults = markResults;
data.markFiles = markFiles;
data.markAmplitudes = markAmplitudes;
end

function mark = solve_mark_point(cfg, targetPeriod, targetAmplitude, startParams, markIndex)
mark = run_flexmod_fmam_task(struct( ...
    'modelType', 'base', ...
    'startParams', startParams, ...
    'goalOrder', {{'period', 'amplitude'}}, ...
    'goals', struct('period', targetPeriod, 'amplitude', targetAmplitude), ...
    'controlledParams', {{'I1', 'ET'}}, ...
    'fmam', struct('dlambdaCap', 1)));
mark.orbit = flexmod.shift_cycle_to_max(mark.orbit);
mark.figureId = 'fig3c';
mark.markIndex = markIndex;
mark.markType = 'amplitude';
mark.targetPeriod = targetPeriod;
mark.targetAmplitude = targetAmplitude;
mark.seedParams = reshape(startParams, 1, []);
end
