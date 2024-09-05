function Qin_new = Function_valve_to_flow(H1, d, Z0, Av)
% ------------------------------------------------------------------------------------------------------------
% Actuator models: transfer valve's setpoints to flow
% Input: orifice's setpoint, water level before the orifice
% Output: inflow rate to the storage tank
% ------------------------------------------------------------------------------------------------------------
%H1 = 0.3005;
%d = 0.3;
%Av = 1;
%Z0 = 0.2;
Cd = 0.65;
g = 9.81;
enta = (H1 - Z0) / d;

if enta > 0
    x0 = Av;
    A0 = d^2 / 4 * (acos(1 - 2*x0) - (1 - 2*x0) .* sin(acos(1 - 2*x0)));
    
    if x0  <  enta
        He = H1 - Z0 - d / 2 * x0;
        y0 = Cd * A0 * (2 * 9.81 * He)^0.5 * 1000;
    else
        y0 = Cd * (H1 - Z0)^1.5 * g^0.5 * A0/ (x0 * d) * 1000;
    end
else
    y0 = 0;
end

% getenv('BLAS_VERSION')
% setenv('BLAS_VERSION','')

Qin_new = y0;
end