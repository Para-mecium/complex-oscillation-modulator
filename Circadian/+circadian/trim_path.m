function trimmed = trim_path(path, mode, cfg)
keep = true(size(path.period));

for i = 2:numel(keep)
    if any(~isfinite(path.paramMatrix(i, :))) || ...
            ~isfinite(path.period(i)) || ~isfinite(path.obsAmplitude(i)) || ...
            ~isfinite(path.obsMax(i))
        keep(i:end) = false;
        break
    end

    paramJump = abs(path.paramMatrix(i, :) - path.paramMatrix(i - 1, :));
    switch lower(mode)
        case 'iso_period'
            if any(paramJump > cfg.fig5b.maxParamJump) || ...
                    abs(path.obsAmplitude(i) - path.obsAmplitude(i - 1)) > cfg.fig5b.maxAmplitudeJump
                keep(i:end) = false;
                break
            end
        case 'iso_maximum'
            if any(paramJump > cfg.fig5d.maxParamJump) || ...
                    abs(path.period(i) - path.period(i - 1)) > cfg.fig5d.maxPeriodJump
                keep(i:end) = false;
                break
            end
        case 'iso_amplitude'
            if any(paramJump > cfg.figS2a.maxParamJump) || ...
                    abs(path.period(i) - path.period(i - 1)) > cfg.figS2a.maxPeriodJump
                keep(i:end) = false;
                break
            end
    end
end

trimmed = path;
fields = {'paramMatrix', 'period', 'obsAmplitude', 'obsMax', 'obsMin', 'lambda', ...
    'directConditionEstimate', 'finalConditionEstimate', 'success'};
for i = 1:numel(fields)
    trimmed.(fields{i}) = path.(fields{i})(keep, :);
end
end
