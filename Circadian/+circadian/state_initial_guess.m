function y0 = state_initial_guess(stateLike, cfg)
series = stateLike.TS_var;
nSamples = size(series, 1);
phaseFraction = cfg.orbit.stateInitPhaseFraction;
sampleIndex = 1 + floor(phaseFraction * max(nSamples - 1, 0));
sampleIndex = min(max(sampleIndex, 1), nSamples);
y0 = series(sampleIndex, :).';
end
