function data = load_ergodic_results(filePath)
%LOAD_ERGODIC_RESULTS Normalize legacy and refactored ergodic saves.

raw = load(filePath);
if isfield(raw, 'schemaVersion') && isfield(raw, 'baseline') && isfield(raw, 'repeats')
    data = normalize_v2(raw, filePath);
else
    data = normalize_legacy(raw, filePath);
end
end

function data = normalize_v2(raw, filePath)
raw.baseline = normalize_orbit_record(raw.baseline);
nodeCount = infer_node_count(raw.baseline.observables.amplitude);
repeats = reshape(raw.repeats, 1, []);
for i = 1:numel(repeats)
    repeats(i) = normalize_repeat_group(repeats(i));
    if ~isfield(repeats(i), 'samples')
        repeats(i).samples = struct([]);
    end
    repeats(i).successfulSamples = collect_successful_samples(repeats(i).samples);
    repeats(i).plottableSamples = collect_plottable_samples(repeats(i).samples);
end

data = struct( ...
    'filePath', filePath, ...
    'schemaVersion', raw.schemaVersion, ...
    'modelName', raw.modelName, ...
    'netName', raw.netName, ...
    'weightPer', raw.weightPer, ...
    'nodeCount', nodeCount, ...
    'requestedRepeatCount', raw.numRepeat, ...
    'baseline', raw.baseline, ...
    'repeats', repeats, ...
    'summary', raw.summary, ...
    'config', raw.config);
end

function data = normalize_legacy(raw, filePath)
baselineObs = normalize_legacy_observables(raw.props(1).prop);
nodeCount = infer_node_count(baselineObs.amplitude);
requestedRepeatCount = raw.numRepeat;

repeats = repmat(struct( ...
    'nPerturbedEdges', 0, ...
    'requestedRepeatCount', requestedRepeatCount, ...
    'attemptCount', 0, ...
    'successCount', 0, ...
    'failureCount', 0, ...
    'samples', struct([]), ...
    'successfulSamples', struct([]), ...
    'plottableSamples', struct([])), 1, max(numel(raw.props) - 1, 0));

for i = 2:numel(raw.props)
    legacySamples = reshape(raw.props(i).prop, 1, []);
    samples = repmat(struct( ...
        'attemptIndex', 0, ...
        'success', true, ...
        'hasOrbit', true, ...
        'isCandidate', false, ...
        'status', 'success', ...
        'orbitStatus', 'converged_periodic_orbit', ...
        'orbitCode', 2, ...
        'message', '', ...
        'sampleSeed', [], ...
        'networkSeed', [], ...
        'perturbationSeed', [], ...
        'perturbation', struct([]), ...
        'networkSummary', struct([]), ...
        'observables', struct([]), ...
        'diagnostics', struct()), 1, numel(legacySamples));

    for j = 1:numel(legacySamples)
        samples(j).attemptIndex = j;
        samples(j).observables = normalize_legacy_observables(legacySamples(j));
    end

    repeats(i - 1).nPerturbedEdges = raw.props(i).n_per;
    repeats(i - 1).attemptCount = numel(samples);
    repeats(i - 1).successCount = numel(samples);
    repeats(i - 1).samples = samples;
    repeats(i - 1).successfulSamples = samples;
    repeats(i - 1).plottableSamples = samples;
end

data = struct( ...
    'filePath', filePath, ...
    'schemaVersion', 1, ...
    'modelName', infer_model_name(filePath), ...
    'netName', '', ...
    'weightPer', [], ...
    'nodeCount', nodeCount, ...
    'requestedRepeatCount', requestedRepeatCount, ...
    'baseline', struct( ...
        'success', true, ...
        'hasOrbit', true, ...
        'isCandidate', false, ...
        'status', 'success', ...
        'orbitStatus', 'converged_periodic_orbit', ...
        'orbitCode', 2, ...
        'message', '', ...
        'TS', {{[], []}}, ...
        'observables', baselineObs, ...
        'diagnostics', struct()), ...
    'repeats', repeats, ...
    'summary', struct([]), ...
    'config', struct([]));
end

function repeat = normalize_repeat_group(repeat)
if isfield(repeat, 'samples')
    samples = reshape(repeat.samples, 1, []);
    for j = 1:numel(samples)
        samples(j) = normalize_orbit_record(samples(j));
    end
    repeat.samples = samples;
end
end

function record = normalize_orbit_record(record)
if ~isfield(record, 'hasOrbit')
    record.hasOrbit = getfield_with_default(record, 'success', false);
end
if ~isfield(record, 'isCandidate')
    record.isCandidate = false;
end
if ~isfield(record, 'orbitCode')
    if record.isCandidate
        record.orbitCode = 1;
    elseif getfield_with_default(record, 'hasOrbit', false)
        record.orbitCode = 2;
    else
        record.orbitCode = 0;
    end
end
if ~isfield(record, 'orbitStatus')
    if record.isCandidate
        record.orbitStatus = 'candidate_periodic_orbit_not_converged';
    elseif getfield_with_default(record, 'hasOrbit', false)
        record.orbitStatus = 'converged_periodic_orbit';
    else
        record.orbitStatus = 'no_periodic_orbit_detected_on_tspan';
    end
end
if ~isfield(record, 'diagnostics')
    record.diagnostics = struct();
end
end

function value = getfield_with_default(s, name, defaultValue)
if isfield(s, name)
    value = s.(name);
else
    value = defaultValue;
end
end

function successfulSamples = collect_successful_samples(samples)
if isempty(samples)
    successfulSamples = struct([]);
    return
end

mask = arrayfun(@(sample) isfield(sample, 'success') && sample.success, samples);
successfulSamples = samples(mask);
end

function plottableSamples = collect_plottable_samples(samples)
if isempty(samples)
    plottableSamples = struct([]);
    return
end

mask = arrayfun(@is_sample_plottable, samples);
plottableSamples = samples(mask);
end

function tf = is_sample_plottable(sample)
success = isfield(sample, 'success') && sample.success;
candidate = isfield(sample, 'isCandidate') && sample.isCandidate;
tf = success || candidate;
end

function observables = normalize_legacy_observables(sample)
observables = struct( ...
    'amplitude', get_field(sample, {'amp', 'amplitude'}), ...
    'period', get_field(sample, {'period', 'peroid'}), ...
    'maxVariable', get_field(sample, {'max', 'maxVariable'}), ...
    'minVariable', get_field(sample, {'min', 'minVariable'}));
end

function value = get_field(sample, names)
value = [];
for i = 1:numel(names)
    if isfield(sample, names{i})
        value = sample.(names{i});
        return
    end
end
end

function nodeCount = infer_node_count(amplitude)
nodeCount = floor(numel(amplitude) / 2);
end

function modelName = infer_model_name(filePath)
[~, name] = fileparts(filePath);
if contains(name, 'GRN')
    modelName = 'GRN';
else
    modelName = 'FHN';
end
end
