sys = cell(1,2);
sys{1} = str2func('dxdt');
sys{2} = str2func('dydt');

function output = dxdt(TS_variable,parameter)
    alpha = parameter(1);
    kappa = parameter(2);
    gamma1 = parameter(3);
    gamma2 = parameter(4); 
    epsilon = parameter(5);
    f11 = parameter(6);
    f12 = parameter(7);

    x = TS_variable(:,1);
    y = TS_variable(:,2);

    output = (alpha+kappa*x.^2./(gamma1+x.^2+gamma2*y)-x)/epsilon+f11*(x-7/4)./(x-7/4+5)+f12*(y-11/4)./(y-11/4+5);
end

function output = dydt(TS_variable,parameter)
    f21 = parameter(8);
    f22 = -parameter(6);
    
    x = TS_variable(:,1);
    y = TS_variable(:,2);
    
    output = 1+x-y+f21*(x-7/4)./(x-7/4+5)+f22*(y-11/4)./(y-11/4+5);
end                                                    