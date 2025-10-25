sys = cell(1,3);
sys{1} = str2func('dp1dt');
sys{2} = str2func('dp2dt');
sys{3} = str2func('dp3dt');

function output = dp1dt(y,Parameters)
    R_C = Parameters(1);
    n = Parameters(2);
    V_th = Parameters(3);
    beta_1 = Parameters(4);
    I_max = Parameters(7);

    V_1 = y(:,1);
    V_3 = y(:,3);

    I_t1 = I_max*V_th^n./(V_th^n+V_3.^n);

    output = beta_1*1e3*(-V_1/R_C + I_t1);
end

function output = dp2dt(y,Parameters)
    R_C = Parameters(1);
    n = Parameters(2);
    V_th = Parameters(3);
    beta_2 = Parameters(5);
    I_max = Parameters(7);

    V_1 = y(:,1);
    V_2 = y(:,2);

    I_t2 = I_max*V_th^n./(V_th^n+V_1.^n);

    output = beta_2*1e3*(-V_2/R_C + I_t2);
end

function output = dp3dt(y,Parameters)
    R_C = Parameters(1);
    n = Parameters(2);
    V_th = Parameters(3);
    beta_3 = Parameters(6);
    I_max = Parameters(7);

    V_2 = y(:,2);
    V_3 = y(:,3);

    I_t3 = I_max*V_th^n./(V_th^n+V_2.^n);

    output = beta_3*1e3*(-V_3/R_C + I_t3);
end

