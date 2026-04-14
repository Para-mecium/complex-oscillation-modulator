function trimmed = trim_path(path, mode, cfg)
keep = true(size(path.period));

for i = 2:numel(keep)
    if any(~isfinite(path.paramMatrix(i, :))) || ...
            ~isfinite(path.period(i)) || ~isfinite(path.amplitude(i))
        keep(i:end) = false;
        break
    end

    paramJump = abs(path.paramMatrix(i, :) - path.paramMatrix(i - 1, :));
    if strcmpi(mode, 'iso_amplitude')
        if any(paramJump > cfg.fig3d.maxParamJump) || ...
                abs(path.period(i) - path.period(i - 1)) > cfg.fig3d.maxPeriodJump
            keep(i:end) = false;
            break
        end
    else
        if any(paramJump > [0.2, 0.25]) || abs(path.amplitude(i) - path.amplitude(i - 1)) > 0.8
            keep(i:end) = false;
            break
        end
    end
end

trimmed = path;
fields = {'paramMatrix', 'period', 'amplitude', 'yMax', 'yMin', 'lambda', ...
    'directConditionEstimate', 'finalConditionEstimate', 'success'};
for i = 1:numel(fields)
    trimmed.(fields{i}) = path.(fields{i})(keep, :);
end
end
