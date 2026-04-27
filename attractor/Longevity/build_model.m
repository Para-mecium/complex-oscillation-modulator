function dydt = build_model(~, y, Parameters)
%BUILD_MODEL RHS for the longevity oscillator.

validateattributes(Parameters, {'numeric'}, {'vector', 'numel', 13}, mfilename, 'Parameters');
validateattributes(y, {'numeric'}, {'vector', 'numel', 4}, mfilename, 'y');

Parameters = reshape(Parameters, 1, []);
y = reshape(y, 1, []);

alphaS = Parameters(1);
alphaH = Parameters(2);
alphaS0 = Parameters(3);
alphaH0 = Parameters(4);
betaS = Parameters(5);
betaH = Parameters(6);
deltam = Parameters(7);
deltaS = Parameters(8);
deltaH = Parameters(9);
KH = Parameters(10);
KS = Parameters(11);
n1 = Parameters(12);
n2 = Parameters(13);

mS = y(:, 1);
mH = y(:, 2);
S = y(:, 3);
H = y(:, 4);

dmSdt = alphaS0 + alphaS * H.^n1 ./ (KH^n1 + H.^n1) - deltam * mS;
dmHdt = alphaH0 + alphaH * KS^n2 ./ (KS^n2 + S.^n2) - deltam * mH;
dSdt = betaS * mS - deltaS * S;
dHdt = betaH * mH - deltaH * H;

dydt = [dmSdt; dmHdt; dSdt; dHdt];
end
