sys = cell(1,2);
sys{1} = str2func('dXdt');
sys{2} = str2func('dYdt');

function output = dXdt(TS_variable,parameter)
    p = 4;
    Kd = 1;
    k1 = 0.05;
    kdx = 0.05;
    bU = 2;
    U = 2;
    K1 = 1;
    KU = 1;
    S = 1;

    X = TS_variable(:,1);
    Y = TS_variable(:,2);

    I1 = parameter(1);
    hU = U * I1^2 / (K1 + I1^2);
    fU = bU * hU^2 / (KU + hU^2);

    output = k1 * S * Kd^p ./ (Kd^p + Y.^p) .* fU - kdx * X;
end

function output = dYdt(TS_variable,parameter)
    Km = 0.1;
    KI = 2;
    ksy = 1;
    kdy = 0.05;
    k2 = 1;

    Y = TS_variable(:,2);

    X = TS_variable(:,1);
    ET = parameter(2);
    output = ksy * X - kdy * Y - k2 * ET .* Y ./ (Km + Y + KI * Y.^2);
end
