function shifted = shift_cycle_to_max(orbit)
t = orbit.t(:);
y = orbit.y;
obs = orbit.obs(:);
[~, idx] = max(obs);

period = orbit.period;
order = [idx:numel(t), 1:idx-1];
shiftedT = t(order) - t(idx);
wrapMask = shiftedT < 0;
shiftedT(wrapMask) = shiftedT(wrapMask) + period;

[shiftedT, sortIdx] = sort(shiftedT);
shiftedY = y(order(sortIdx), :);
shiftedObs = obs(order(sortIdx));

shifted = orbit;
shifted.t = shiftedT;
shifted.y = shiftedY;
shifted.obs = shiftedObs;
end
