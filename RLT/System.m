sys = cell(1, 6);
sys{1} = str2func('dm1dt');
sys{2} = str2func('dm2dt');
sys{3} = str2func('dm3dt');
sys{4} = str2func('dp1dt');
sys{5} = str2func('dp2dt');
sys{6} = str2func('dp3dt');

function output = dm1dt(y, Parameters)
    n = Parameters(1);
    alpha = Parameters(2);
    alpha0 = Parameters(3);

    m1 = y(:, 1);
    p3 = y(:, 6);

    output = -m1 + alpha ./ (1 + p3 .^ n) + alpha0;
end

function output = dm2dt(y, Parameters)
    n = Parameters(1);
    alpha = Parameters(2);
    alpha0 = Parameters(3);

    m2 = y(:, 2);
    p1 = y(:, 4);

    output = -m2 + alpha ./ (1 + p1 .^ n) + alpha0;
end

function output = dm3dt(y, Parameters)
    n = Parameters(1);
    alpha = Parameters(2);
    alpha0 = Parameters(3);

    m3 = y(:, 3);
    p2 = y(:, 5);

    output = -m3 + alpha ./ (1 + p2 .^ n) + alpha0;
end

function output = dp1dt(y, Parameters)
    beta = Parameters(4);
    phi1 = Parameters(5);

    m1 = y(:, 1);
    p1 = y(:, 4);

    output = beta * (m1 - phi1 * p1);
end

function output = dp2dt(y, Parameters)
    beta = Parameters(4);
    phi2 = Parameters(6);

    m2 = y(:, 2);
    p2 = y(:, 5);

    output = beta * (m2 - phi2 * p2);
end

function output = dp3dt(y, Parameters)
    beta = Parameters(4);

    m3 = y(:, 3);
    p3 = y(:, 6);

    output = beta * (m3 - p3);
end
