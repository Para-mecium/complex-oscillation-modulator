function dydt = build_model(~, y, Parameters)
%BUILD_MODEL RHS for the repressilator-like transistor circuit.

validateattributes(Parameters, {'numeric'}, {'vector', 'numel', 7}, mfilename, 'Parameters');
validateattributes(y, {'numeric'}, {'vector', 'numel', 3}, mfilename, 'y');

Parameters = reshape(Parameters, 1, []);
y = reshape(y, 1, []);

R_C = Parameters(1);
n = Parameters(2);
V_th = Parameters(3);
beta_1 = Parameters(4);
beta_2 = Parameters(5);
beta_3 = Parameters(6);
I_max = Parameters(7);

V_1 = y(:, 1);
V_2 = y(:, 2);
V_3 = y(:, 3);

I_t1 = I_max * V_th^n ./ (V_th^n + V_3.^n);
I_t2 = I_max * V_th^n ./ (V_th^n + V_1.^n);
I_t3 = I_max * V_th^n ./ (V_th^n + V_2.^n);

dp1dt = beta_1 * 1e3 * (-V_1 / R_C + I_t1);
dp2dt = beta_2 * 1e3 * (-V_2 / R_C + I_t2);
dp3dt = beta_3 * 1e3 * (-V_3 / R_C + I_t3);

dydt = [dp1dt; dp2dt; dp3dt];
end
