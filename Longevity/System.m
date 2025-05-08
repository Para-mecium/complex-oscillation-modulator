sys = cell(1,4);
sys{1} = str2func('dmSdt');
sys{2} = str2func('dmHdt');
sys{3} = str2func('dSdt');
sys{4} = str2func('dHdt');

function output = dmSdt(TS_variable,parameter)
    alphaS = parameter(1);
    alphaS0 = parameter(3);
    deltam = parameter(7);
    KH = parameter(10);
    n1 = parameter(12);

    mS = TS_variable(:,1);
    H = TS_variable(:,4);

    output = alphaS0+alphaS*H.^n1./(KH^n1+H.^n1)-deltam*mS;
end

function output = dmHdt(TS_variable,parameter)
    alphaH = parameter(2);
    alphaH0 = parameter(4);
    deltam = parameter(7);
    KS = parameter(11);
    n2 = parameter(13);

    mH = TS_variable(:,2);
    S = TS_variable(:,3);

    output = alphaH0+alphaH*KS^n2./(KS^n2+S.^n2)-deltam*mH;
end

function output = dSdt(TS_variable,parameter)
    beta1 = parameter(5);
    deltaS = parameter(8);

    mS = TS_variable(:,1);
    S = TS_variable(:,3);

    output = beta1*mS-deltaS*S;
end

function output = dHdt(TS_variable,parameter)
    beta2 = parameter(6);
    deltaH = parameter(9);

    mH = TS_variable(:,2);
    H = TS_variable(:,4);

    output = beta2*mH-deltaH*H;
end