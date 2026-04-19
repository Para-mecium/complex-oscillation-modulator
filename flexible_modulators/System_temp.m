sys = cell(1,2);
sys{1} = str2func('dXdt');
sys{2} = str2func('dYdt');

function output = dXdt(TS_variable,parameter)
    p = 4;
    Kd = 1;
    k1_ref = 0.05;
    kdx_ref = 0.05;
    bU = 2;
    U = 2;
    K1 = 1;
    KU = 1;
    I1 = 1;
    S = 1.5;
    R = 8.314;
    Tref = 298;
    Ek1 = 87.73e3;
    Ekdx = 36.43e3;

    X = TS_variable(:,1);
    Y = TS_variable(:,2);

    Temperature = parameter(2);
    hU = U * I1^2 / (K1 + I1^2);
    fU = bU * hU^2 / (KU + hU^2);

    k1 = k1_ref * exp(Ek1 / R * (1 / Tref - 1 ./ Temperature));
    kdx = kdx_ref * exp(Ekdx / R * (1 / Tref - 1 ./ Temperature));

    output = k1 * S * Kd^p ./ (Kd^p + Y.^p) .* fU - kdx * X;
end

function output = dYdt(TS_variable,parameter)
    Km = 0.1;
    KI = 2;
    ksy_ref = 1;
    kdy_ref = 0.05;
    k2_ref = 1;
    R = 8.314;
    Tref = 298;
    Eksy = 61.12e3;
    Ekdy = 47.84e3;
    Ek2 = 63.31e3;

    Y = TS_variable(:,2);

    X = TS_variable(:,1);
    ET = parameter(1);
    Temperature = parameter(2);

    ksy = ksy_ref * exp(Eksy / R * (1 / Tref - 1 ./ Temperature));
    kdy = kdy_ref * exp(Ekdy / R * (1 / Tref - 1 ./ Temperature));
    k2 = k2_ref * exp(Ek2 / R * (1 / Tref - 1 ./ Temperature));

    output = ksy * X - kdy * Y - k2 * ET .* Y ./ (Km + Y + KI * Y.^2);
end
