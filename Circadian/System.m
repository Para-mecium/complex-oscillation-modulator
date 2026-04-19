sys = cell(1, 3);
sys{1} = str2func('dMdt');
sys{2} = str2func('dPc_dt');
sys{3} = str2func('dPn_dt');

function output = dMdt(TS_variable, parameter)
    beta = 0.1572;

    M = TS_variable(:, 1);
    Pn = TS_variable(:, 3);

    Kd = parameter(1);
    AT = parameter(2);
    Af = free_active(AT, Kd, Pn);

    output = beta * (Af ./ AT - M);
end

function output = dPc_dt(TS_variable, ~)
    beta = 0.1572;

    M = TS_variable(:, 1);
    Pc = TS_variable(:, 2);

    output = beta * (M - Pc);
end

function output = dPn_dt(TS_variable, ~)
    beta = 0.1572;

    Pc = TS_variable(:, 2);
    Pn = TS_variable(:, 3);

    output = beta * (Pc - Pn);
end

function value = free_active(AT, Kd, Pn)
    delta = AT - Pn - Kd;
    value = 0.5 * (delta + sqrt(delta.^2 + 4 * Kd * AT));
end
