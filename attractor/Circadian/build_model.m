function dydt = build_model(~, y, Parameters)
%BUILD_MODEL RHS for the circadian oscillator.

validateattributes(Parameters, {'numeric'}, {'vector', 'numel', 2}, mfilename, 'Parameters');
validateattributes(y, {'numeric'}, {'vector', 'numel', 3}, mfilename, 'y');

Parameters = reshape(Parameters, 1, []);
y = reshape(y, 1, []);

beta = 0.1572;
Kd = Parameters(1);
AT = Parameters(2);

M = y(:, 1);
Pc = y(:, 2);
Pn = y(:, 3);

Af = free_active(AT, Kd, Pn);

dMdt = beta * (Af ./ AT - M);
dPc_dt = beta * (M - Pc);
dPn_dt = beta * (Pc - Pn);

dydt = [dMdt; dPc_dt; dPn_dt];
end

function value = free_active(AT, Kd, Pn)
delta = AT - Pn - Kd;
value = 0.5 * (delta + sqrt(delta.^2 + 4 * Kd * AT));
end
