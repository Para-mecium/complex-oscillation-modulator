function data = build_fig5b_marks(cfg)
if nargin < 1
    cfg = struct();
end
circadian.ensure_paths();
cfg = circadian.merge_config(default_config(), cfg);

markAmplitudes = cfg.fig5b.markAmplitudes;
markSeeds = cfg.fig5b.markSeeds;
assert(size(markSeeds, 1) == numel(markAmplitudes) && size(markSeeds, 2) == 2, ...
    'cfg.fig5b.markSeeds must be an N x 2 matrix aligned with cfg.fig5b.markAmplitudes.');

markResults = cell(1, numel(markAmplitudes));
markFiles = cell(1, numel(markAmplitudes));
for i = 1:numel(markAmplitudes)
    markFiles{i} = circadian.mark_cache_file(cfg, 'fig5b', sprintf('amplitude_%02d', i));
    targetAmplitude = markAmplitudes(i);
    seedParams = markSeeds(i, :);
    markResults{i} = circadian.cache_get_or_create(markFiles{i}, ...
        @() solve_mark_point(cfg, cfg.fig5b.markPeriod, targetAmplitude, seedParams, i));
end

data = struct();
data.markResults = markResults;
data.markFiles = markFiles;
data.markAmplitudes = markAmplitudes;
end

function mark = solve_mark_point(cfg, targetPeriod, targetAmplitude, startParams, markIndex)
mark = run_modulation_task(struct( ...
    'startParams', startParams, ...
    'goalOrder', {{'period', 'amplitude'}}, ...
    'goals', struct('period', targetPeriod, 'amplitude', targetAmplitude), ...
    'controlledParams', {{'Kd', 'AT'}}, ...
    'fmam', struct('dlambdaCap', 1)));
mark.orbit = circadian.shift_cycle_to_max(mark.orbit);
mark.figureId = 'fig5b';
mark.markIndex = markIndex;
mark.markType = 'amplitude';
mark.targetPeriod = targetPeriod;
mark.targetAmplitude = targetAmplitude;
mark.seedParams = reshape(startParams, 1, []);
end
