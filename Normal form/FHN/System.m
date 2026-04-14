sys = cell(1,2);
sys{1} = str2func('dVdt');
sys{2} = str2func('dWdt');

function output = dVdt(TS_variable,parameter)
    theta = parameter(1);
    I = parameter(3);
    f11 = parameter(5);
    f12 = parameter(6);

    V = TS_variable(:,1);
    W = TS_variable(:,2);

    output = V.*(V-theta).*(1-V)-W+I+f11*(V-0.3068)+f12*(W-0.1227);
end

function output = dWdt(TS_variable,parameter)
    epsilon = parameter(4);
    gamma = parameter(2);
    f21 = parameter(7);
    f22 = -parameter(5);

    V = TS_variable(:,1);
    W = TS_variable(:,2);

    output = epsilon*(V-gamma*W)+f21*(V-0.3068)+f22*(W-0.1227);
end