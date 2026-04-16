function dydt = po_call_odefun(odefun, t, y, parameter)
narginHandle = nargin(odefun);
if narginHandle == 2
    dydt = odefun(t, y);
else
    dydt = odefun(t, y, parameter);
end

dydt = dydt(:);
end
