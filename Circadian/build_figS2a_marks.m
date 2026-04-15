function data = build_figS2a_marks(cfg)
if nargin < 1
    cfg = struct();
end
circadian.ensure_paths();
cfg = circadian.merge_config(default_config(), cfg);

targetAmplitude = cfg.figS2a.markAmplitude;
markPeriods = cfg.figS2a.markPeriods;
markSeeds = cfg.figS2a.markSeeds;
assert(size(markSeeds, 1) == numel(markPeriods) && size(markSeeds, 2) == 2, ...
    'cfg.figS2a.markSeeds must be an N x 2 matrix aligned with cfg.figS2a.markPeriods.');

markResults = cell(1, numel(markPeriods));
markFiles = cell(1, numel(markPeriods));
for i = 1:numel(markPeriods)
    markFiles{i} = circadian.mark_cache_file(cfg, 'figS2a', sprintf('period_%02d', i));
    targetPeriod = markPeriods(i);
    seedParams = markSeeds(i, :);
    markResults{i} = circadian.cache_get_or_create(markFiles{i}, ...
        @() solve_mark_point(cfg, targetAmplitude, targetPeriod, seedParams, i));
end

data = struct();
data.markResults = markResults;
data.markFiles = markFiles;
data.markAmplitude = targetAmplitude;
data.markPeriods = markPeriods;
end

function mark = solve_mark_point(cfg, targetAmplitude, targetPeriod, startParams, markIndex)
mark = run_circadian_fmam_task(struct( ...
    'startParams', startParams, ...
    'goalOrder', {{'amplitude', 'period'}}, ...
    'goals', struct('amplitude', targetAmplitude, 'period', targetPeriod), ...
    'controlledParams', {{'Kd', 'AT'}}, ...
    'fmam', struct('dlambdaCap', 1)));
mark.orbit = circadian.shift_cycle_to_max(mark.orbit);
mark.figureId = 'figS2a';
mark.markIndex = markIndex;
mark.markType = 'period';
mark.targetAmplitude = targetAmplitude;
mark.targetPeriod = targetPeriod;
mark.seedParams = reshape(startParams, 1, []);
end
