function data = build_fig3b_data(cfg)
if nargin < 1
    cfg = struct();
end
flexmod.ensure_paths();
cfg = flexmod.merge_config(default_config(), cfg);

model = flexmod.build_model(cfg);
I1Values = cfg.grid.I1Values;
ETValues = cfg.grid.ETValues;

periodGrid = NaN(numel(ETValues), numel(I1Values));
amplitudeGrid = NaN(numel(ETValues), numel(I1Values));
for i = 1:numel(I1Values)
    for j = 1:numel(ETValues)
        disp(['computing: (' num2str(i) ',' num2str(j) ')...'])
        orbit = flexmod.find_orbit(model, [I1Values(i), ETValues(j)], cfg, cfg.initialState);
        if orbit.success
            periodGrid(j, i) = orbit.period;
            amplitudeGrid(j, i) = orbit.amplitude;
        end
    end
end

fixedET = model.defaultParams(2);
varyI1 = build_single_parameter_series(model, cfg, [I1Values(:), fixedET * ones(numel(I1Values), 1)], 'I1');
fixedI1 = cfg.grid.fixedI1;
varyET = build_single_parameter_series(model, cfg, [fixedI1 * ones(numel(ETValues), 1), ETValues(:)], 'ET');

seriesOrbits = cell(1, numel(cfg.grid.seriesET));
for i = 1:numel(cfg.grid.seriesET)
    seriesOrbits{i} = flexmod.find_orbit(model, [cfg.grid.fixedI1, cfg.grid.seriesET(i)], cfg, cfg.initialState);
    if seriesOrbits{i}.success
        seriesOrbits{i} = flexmod.shift_cycle_to_max(seriesOrbits{i});
    end
end

data = struct();
data.I1Values = I1Values;
data.ETValues = ETValues;
data.periodGrid = periodGrid;
data.amplitudeGrid = amplitudeGrid;
data.varyI1 = varyI1;
data.varyET = varyET;
data.seriesET = cfg.grid.seriesET;
data.seriesOrbits = seriesOrbits;
end

function series = build_single_parameter_series(model, cfg, paramPairs, varyingName)
n = size(paramPairs, 1);
period = NaN(n, 1);
amplitude = NaN(n, 1);
for i = 1:n
    orbit = flexmod.find_orbit(model, paramPairs(i, :), cfg, cfg.initialState);
    if orbit.success
        period(i) = orbit.period;
        amplitude(i) = orbit.amplitude;
    end
end

series = struct('params', paramPairs, 'period', period, 'amplitude', amplitude, 'varyingName', varyingName);
end
