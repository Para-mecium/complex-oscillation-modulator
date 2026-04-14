function data = build_fig5c_data(cfg)
if nargin < 1
    cfg = struct();
end
circadian.ensure_paths();
cfg = circadian.merge_config(default_config(), cfg);

markData = build_fig5b_marks(cfg);
markResults = markData.markResults;

series = cell(1, numel(markResults));
for i = 1:numel(markResults)
    orbit = circadian.shift_cycle_to_max(markResults{i}.orbit);
    series{i} = struct( ...
        't', orbit.t, ...
        'Ptot', orbit.obs, ...
        'Pc', orbit.y(:, 2), ...
        'Pn', orbit.y(:, 3), ...
        'amplitude', markResults{i}.measures.amplitude, ...
        'params', markResults{i}.params);
end

data = struct();
data.markResults = markResults;
data.markAmplitudes = markData.markAmplitudes;
data.series = series;
end
