function dydt = build_model(~, y, activeParams)
%BUILD_MODEL RHS for the base flexible-modulator model.

validateattributes(activeParams, {'numeric'}, {'vector', 'numel', 2}, mfilename, 'activeParams');
validateattributes(y, {'numeric'}, {'vector', 'numel', 2}, mfilename, 'y');

activeParams = reshape(activeParams, 1, []);
y = reshape(y, 1, []);

p = 4;
Kd = 1;
k1 = 0.05;
kdx = 0.05;
bU = 2;
U = 2;
K1 = 1;
KU = 1;
S = 1;
Km = 0.1;
KI = 2;
ksy = 1;
kdy = 0.05;
k2 = 1;

X = y(:, 1);
Y = y(:, 2);
I1 = activeParams(1);
ET = activeParams(2);

hU = U * I1^2 / (K1 + I1^2);
fU = bU * hU^2 / (KU + hU^2);

dXdt = k1 * S * Kd^p ./ (Kd^p + Y.^p) .* fU - kdx * X;
dYdt = ksy * X - kdy * Y - k2 * ET .* Y ./ (Km + Y + KI * Y.^2);

dydt = [dXdt; dYdt];
end
