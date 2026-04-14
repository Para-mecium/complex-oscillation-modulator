function log_orbit_result(label, result, indent)
%LOG_ORBIT_RESULT Print a human-readable orbit-detection summary.

if nargin < 3 || isempty(indent)
    indent = '  ';
end

fprintf('%s%s: %s\n', indent, label, result.status);
detailLines = build_detail_lines(result);
for i = 1:numel(detailLines)
    fprintf('%s  %s\n', indent, detailLines{i});
end
end

function lines = build_detail_lines(result)
lines = {};

if isfield(result, 'orbitStatus') && ~isempty(result.orbitStatus) && ...
        ~strcmp(result.orbitStatus, result.status)
    if isfield(result, 'orbitCode') && ~isnan(result.orbitCode)
        lines{end + 1} = sprintf('orbitStatus=%s (code=%g)', result.orbitStatus, result.orbitCode); %#ok<AGROW>
    else
        lines{end + 1} = sprintf('orbitStatus=%s', result.orbitStatus); %#ok<AGROW>
    end
end

if isfield(result, 'message') && ~isempty(result.message)
    lines{end + 1} = sprintf('message=%s', result.message); %#ok<AGROW>
end

diagnostics = struct();
if isfield(result, 'diagnostics') && isstruct(result.diagnostics)
    diagnostics = result.diagnostics;
end

lines = [lines, build_tail_lines(diagnostics)]; %#ok<AGROW>
lines = [lines, build_strict_check_lines(diagnostics)]; %#ok<AGROW>
lines = [lines, build_candidate_lines(result, diagnostics)]; %#ok<AGROW>
end

function lines = build_tail_lines(diagnostics)
lines = {};
thresholds = get_thresholds(diagnostics);
if isempty(thresholds)
    return
end

if isfield(diagnostics, 'numCrossings') && diagnostics.numCrossings < thresholds.minCrossings
    lines{end + 1} = sprintf('crossings after transient=%d < minCrossings=%d', ... %#ok<AGROW>
        diagnostics.numCrossings, thresholds.minCrossings);
end

if isfield(diagnostics, 'tailRange') && isfield(diagnostics, 'tailStd')
    smallRange = diagnostics.tailRange <= thresholds.nonoscAmpTol;
    smallStd = diagnostics.tailStd <= thresholds.nonoscStdTol;
    if smallRange || smallStd
        lines{end + 1} = sprintf( ... %#ok<AGROW>
            'tail oscillation: range=%s (tol=%s), std=%s (tol=%s)', ...
            fmt_num(diagnostics.tailRange), fmt_num(thresholds.nonoscAmpTol), ...
            fmt_num(diagnostics.tailStd), fmt_num(thresholds.nonoscStdTol));
    end
end
end

function lines = build_strict_check_lines(diagnostics)
lines = {};
thresholds = get_thresholds(diagnostics);
if isempty(thresholds)
    return
end

periods = get_numeric_vector(diagnostics, 'periods');
poincareResidual = get_numeric_vector(diagnostics, 'poincareResidual');
periodRelChange = get_numeric_vector(diagnostics, 'periodRelChange');
amplitudeRelChange = get_numeric_vector(diagnostics, 'amplitudeRelChange');
nCycles = numel(periods);
m = thresholds.consecutiveCycles;

if nCycles < m
    if nCycles > 0
        lines{end + 1} = sprintf('strict convergence window unavailable: nCycles=%d < required=%d', nCycles, m); %#ok<AGROW>
    end
    return
end

idx = (nCycles - m + 1):(nCycles - 1);
lines = append_check_line(lines, ...
    'strict poincare residual', poincareResidual(idx), thresholds.strict.poincareTol);
lines = append_check_line(lines, ...
    'strict period relative change', periodRelChange(idx), thresholds.strict.periodTol);
lines = append_check_line(lines, ...
    'strict amplitude relative change', amplitudeRelChange(idx), thresholds.strict.amplitudeTol);
end

function lines = build_candidate_lines(result, diagnostics)
lines = {};
if ~(isfield(result, 'status') && strcmp(result.status, 'candidate'))
    return
end

if isfield(diagnostics, 'poincareResidual') && ~isempty(diagnostics.poincareResidual)
    lastResidual = diagnostics.poincareResidual(end);
    if isfield(diagnostics, 'candidateResidualTol')
        lines{end + 1} = sprintf('candidate residual check: last=%s, tolerance=%s', ... %#ok<AGROW>
            fmt_num(lastResidual), fmt_num(diagnostics.candidateResidualTol));
    else
        lines{end + 1} = sprintf('candidate residual check: last=%s', fmt_num(lastResidual)); %#ok<AGROW>
    end
end

if isfield(diagnostics, 'cvPeriod_recent') || isfield(diagnostics, 'cvAmplitude_recent')
    cvPeriod = NaN;
    cvAmplitude = NaN;
    if isfield(diagnostics, 'cvPeriod_recent')
        cvPeriod = diagnostics.cvPeriod_recent;
    end
    if isfield(diagnostics, 'cvAmplitude_recent')
        cvAmplitude = diagnostics.cvAmplitude_recent;
    end
    lines{end + 1} = sprintf('candidate variability: cvPeriod=%s, cvAmplitude=%s', ... %#ok<AGROW>
        fmt_num(cvPeriod), fmt_num(cvAmplitude));
end
end

function lines = append_check_line(lines, label, values, tolerance)
if isempty(values)
    return
end

maxValue = max(values);
if maxValue > tolerance
    lines{end + 1} = sprintf('%s failed: recent max=%s > tol=%s', ... %#ok<AGROW>
        label, fmt_num(maxValue), fmt_num(tolerance));
end
end

function thresholds = get_thresholds(diagnostics)
thresholds = [];
if isfield(diagnostics, 'thresholds') && isstruct(diagnostics.thresholds)
    thresholds = diagnostics.thresholds;
end
end

function values = get_numeric_vector(s, fieldName)
values = [];
if ~isfield(s, fieldName) || isempty(s.(fieldName))
    return
end

values = s.(fieldName)(:);
values = values(isfinite(values));
end

function text = fmt_num(value)
if isempty(value) || ~isscalar(value) || ~isfinite(value)
    text = 'NaN';
else
    text = sprintf('%.3e', value);
end
end
