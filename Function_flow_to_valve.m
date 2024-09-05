function Avalve = Function_flow_to_valve(H1, d, Z0, Qin_new)
% ------------------------------------------------------------------------------------------------------------
% Actuator models: transfer valve's setpoints to flow
% Input: inflow rate to the storage tank, level before the orifice
% Output: orifice's setpoint
% ------------------------------------------------------------------------------------------------------------
Cd = 0.65;
g = 9.81;
enta = min(round((H1 - Z0) / d,2),1);

if enta > 0
    x01 = 0:1/100:enta;
    N1 = size(x01, 2);
    A01 = d^2 / 4 * (acos(1 - 2*x01) - (1 - 2*x01) .* sin(acos(1 - 2*x01)));
    He = H1 - Z0 - d / 2 * x01;
    y01(1) = 0;
    y01(2:N1) = Cd * A01(2:N1) .* (2 * 9.81 * He(2:N1)).^0.5 * 1000;
    % yi1 = Qin_new;
    delta1 = abs(y01-Qin_new);
    id = find(delta1 == min(delta1), 1, 'first');
    x1 = x01(id);
    
    % getenv('BLAS_VERSION')
    % setenv('BLAS_VERSION','')
    
    x02 = enta:1/100:1;
    A02 = d^2 / 4 * (acos(1 - 2*x02) - (1 - 2*x02) .* sin(acos(1 - 2*x02)));
    y02 = Cd * (H1 - Z0)^1.5 * g^0.5 * A02 ./ (x02 * d) * 1000;
    delta2 = abs(y02-Qin_new);
    id = find(delta2 == min(delta2), 1, 'last');
    x2 = x02(id);
    
    
    if min(delta1) <= min(delta2)
        Avalve = x1;
    else
        Avalve = x2;
    end
    
else
    Avalve = 1;
end

end