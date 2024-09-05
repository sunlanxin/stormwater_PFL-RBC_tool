function [Ao, Ap1, Ap2] = control_model_wet(RTCmodel, StatePara, TankPara, ControlPara, PredictPara, SetpointPara, timedate)
% ------------------------------------------------------------------------------------------------------------
% This function is used to generate a control strategy during wet period
% Inputs: system physical properties, real-time system states, rainfall forecast data
% Outputs: setpoints for the inlet orifice
% RTCmodel 1: RBC
% RTCmodel 2: PFL-RBC
% RTCmodel 3: MPC
% ------------------------------------------------------------------------------------------------------------

% ------------------------------------------------------------------------------------------------------------
% Geting the inputs of system physical properties 
% ------------------------------------------------------------------------------------------------------------
Hmin = TankPara.Hmin; % 闸门初始状态
Hmax = TankPara.Hmax;
As = TankPara.As;
d = TankPara.d;
Z0 = TankPara.Z0;
constep = ControlPara.constep;
% ------------------------------------------------------------------------------------------------------------
% Geting the inputs of state variables 
% ------------------------------------------------------------------------------------------------------------
Qrunoff = StatePara.Qrunoff; 
H1 = StatePara.H1;
h = StatePara.h; 
Qin = StatePara.Qin; 
Qout = StatePara.Qout; 
Qover = StatePara.Qover; 
Quse = StatePara.Quse; 
Q = StatePara.Q; 
Aorifice = SetpointPara.Aorifice;
Apump1 = SetpointPara.Apump1; 
Apump2 = SetpointPara.Apump2; 
Qmax = ControlPara.Qmax;
% ------------------------------------------------------------------------------------------------------------
% Geting the inputs of rainfall forecast data
% ------------------------------------------------------------------------------------------------------------
QP = PredictPara.QP; 

switch RTCmodel
    case 1 % RBC
        Rule = 2; 
        
    case 2 % PFL-RBC
        Rule = 1;
        
    case 3 % MPC
        Tc = ControlPara.Tc;
        Tp = ControlPara.Tp;
        Nc = Tc * 60 / constep;
        Np = Tp * 60 / constep;
        Qin_new = Function_mpc(h, Qrunoff, Q, Hmax, Hmin, Qmax, constep, As, Nc, QP, timedate);
        if Qin_new >= 0.9 * Qrunoff
            Ao = 1;
        else
            Ao = Function_flow_to_valve(H1, d, Z0, Qin_new); % transfer the flow to the orifice's setpoint
        end
        Rule = 0;
end

switch Rule
    case 1
        if h < Hmax
            if Qrunoff >= Qmax
                Qin_new = Qrunoff - Qmax;
                Ao = Function_flow_to_valve(H1, d, Z0, Qin_new);
                if Qin_new == Qrunoff
                    Ao = 1;
                end
            else
                Ao = 0;
            end
        else
            Ao = 0;
        end
        
    case 2
        if h < Hmax
            Ao = 1;
        else
            Ao = 0;
        end
        
    case 0
        
end

Ap1 = 0;
Ap2 = 0;

end


function Qin_new = Function_mpc(H0, Qrunoff, Q, Hmax, Hmin, Qmax, constep, As, Nc, QP, timedate)
% ------------------------------------------------------------------------------------------------------------
% This function is used to generate a decision sequence for model predictive control during wet periods
% Inputs: predictive system states, rainfall forecast data
% Outputs: the controlled inflow rate (Qin)
% ------------------------------------------------------------------------------------------------------------
In = QP;
InitPara = struct('H0', H0, 'Q0', Q, 'x0', Q/Qrunoff);
LimitPara = struct('Hmax', Hmax, 'Hmin', Hmin, 'Qmax', Qmax);
% ------------------------------------------------------------------------------------------------------------
% Online optimization based on the Genetic Algorithm (GA)
% ------------------------------------------------------------------------------------------------------------
for k = 1 : Nc
    params{k} = char(strcat("x", string(k)));
    lb(k) = 0;
    ub(k) = 1;
end
nvars = length(params);
popSize = 100;
crossRate = 0.8;
maxGen = 20;
gaOpts = optimoptions('ga', 'PopulationSize', popSize, 'CrossoverFraction', crossRate, 'MaxGenerations', maxGen);
functn = @(x)ObjFunc(x, In, InitPara, LimitPara, constep, As, Nc);
[x, J] = ga(functn, nvars, [], [], [], [], lb, ub, [], gaOpts);

% ------------------------------------------------------------------------------------------------------------
% 输出最优的决策序列
% ------------------------------------------------------------------------------------------------------------
disp([timedate, 'Optimal decision sequence: ', num2str(x)]);
disp([timedate, 'Optimal fitness: ', num2str(J)]);
Qin_new = (1 - roundn(x(1),-2)) * Qrunoff;
end

function J = ObjFunc(x, In, InitPara, LimitPara, constep, As, Nc)
% ------------------------------------------------------------------------------------------------------------
% This function is used to calculate the fitness function (objective function) for the optimisation problem in the mpc algorithm
% The objective function is set as to minimize the total outflow, the flow variation, and the punishment of exceeding constraints
% J = Σ|Q| + Σ|x-x0| + Punishment
% ------------------------------------------------------------------------------------------------------------
x0 = InitPara.x0; 
Q(1) = InitPara.Q0;
h(1) = InitPara.H0;
Hmax = LimitPara.Hmax;
Hmin = LimitPara.Hmin;
Qmax = LimitPara.Qmax;
M1 = 100;
M2 = 100;
M3 = 1;

for t = 1 : Nc   
    % ------------------------------------------------------------------------------------------------------------
    % Calculating the lower boundary of flow constraint - Qallow, based on a peak inflow-based prediction method
    % ------------------------------------------------------------------------------------------------------------
    hreq(t) = Function_Vreq(t, In, InitPara.H0, LimitPara, constep, As); 
    Vreq(t) = min(hreq(t), Hmax - Hmin) * As;
    Vleft(t) = max((Hmax - h(t)), 0) * As;
    if Vleft(t) > Vreq(t)
        Qalow(t) = max(In(t) - (Vleft(t) - Vreq(t))/(constep * 60 / 1000), 0); % Qallow < Q < Qmax
        Qalow(t) = min(Qalow(t), Qmax);
    else
        Qalow(t) = In(t);
    end
    % ------------------------------------------------------------------------------------------------------------
    % System state transition: linear assumption
    % ------------------------------------------------------------------------------------------------------------
    Q(t) = roundn(x(t),-2) * In(t);
    h(t+1) = h(t) + (In(t) - Q(t)) * (constep * 60 / 1000) / As;
    
end

% ------------------------------------------------------------------------------------------------------------
% Calculating the terms in objective function
% ------------------------------------------------------------------------------------------------------------
for t = 1 : Nc   
    varia_Q(t) = Q(t);    
    
    varia_A(t) = M3 * abs(roundn(x(t),-2) - x0);
    x0 = roundn(x(t),-2);    
    
    if h(t+1) < Hmin
        P1(t) = M1 * abs(Hmin - h(t+1)) * (As / (constep * 60 / 1000));
    elseif  h(t+1) > Hmax
        P1(t) = M1 * abs(h(t+1) - Hmax) * (As / (constep * 60 / 1000));
        Qover(t) = max(h(t+1) - Hmax, 0) * (As / (constep * 60 / 1000)); % 溢流量
        Q(t) = Q(t) + Qover(t);
    else
        P1(t) = 0;
    end
    if  Q(t) < Qalow(t) || Q(t) > Qmax
        P2(t) = M2 * max(Qalow(t) - Q(t), 0) + M2 * max(Q(t) - Qmax, 0);
    else
        P2(t) = 0;
    end    
    PM(t) = P1(t) + P2(t);
end

J = sum(varia_Q) + sum(varia_A) + sum(PM);

end


function hreq = Function_Vreq(t, In, H0, LimitPara, constep, As)
% ------------------------------------------------------------------------------------------------------------
% This function is used to calculate the required storage volume and water% level (hreq) at the t-th control period.
% Within the prediction horizon, hreq is used to limit the Qin before the flood peaks , mainly at the early stage of a rainfall, to reserve sufficient space for large flood peaks
% t: the current time step
% In: Incoming total inflow process in the prediction horizon
% LimtPara: constraints
% Nc, Np: control level, prediction level
% constep: simulation time step
% As: bottom area of the tank
% ------------------------------------------------------------------------------------------------------------
hreq = 0;
Qmax = LimitPara.Qmax;
Hmax = LimitPara.Hmax;
Wmax = (Hmax-H0)*As;
Inp(:) = In(t:end);
m = find(Inp(:) == max(Inp(:)), 1);

if sum(Inp) * (constep * 60 / 1000) < Wmax
    hreq = 0;
else
    if Inp(m) > Qmax
        if t ~= m 
            for i = 1:length(Inp)
                if Inp(i) > Qmax
                    hreq = hreq + (Inp(i) - Qmax) * (constep * 60 / 1000) / As;
                end
            end
        else
            hreq = 0;
        end
        
    elseif Inp(m) > Qmax * 0.5 
        if t ~= m 
            for i = 1:length(Inp)
                if Inp(i) > Qmax
                    n = n + 1;
                    hreq = hreq + (Inp(i) - Qmax * 0.5) * (constep * 60 / 1000) / As;
                end
            end
        else
            hreq = 0;
        end
        
    else 
        if t ~= m
            hreq = sum(Inp(:)) / length(Inp) * constep * 60 / 1000 / As;
        else
            hreq = 0;
        end
    end
end

end