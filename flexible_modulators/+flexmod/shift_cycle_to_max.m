function shifted = shift_cycle_to_max(orbit)
t = orbit.t(:);
y = orbit.y;
protein = y(:, 2);
[~, idx] = max(protein);

period = orbit.period;
order = [idx:numel(t), 1:idx-1];
shiftedT = t(order) - t(idx);
wrapMask = shiftedT < 0;
shiftedT(wrapMask) = shiftedT(wrapMask) + period;

[shiftedT, sortIdx] = sort(shiftedT);
shiftedY = y(order(sortIdx), :);

shifted = orbit;
shifted.t = shiftedT;
shifted.y = shiftedY;
end
