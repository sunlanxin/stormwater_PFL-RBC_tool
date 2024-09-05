function [Ao, Ap1, Ap2] = control_model_dry(RTCmodel, StatePara, TankPara, ControlPara, PredictPara, SetpointPara, td, ~)
% ------------------------------------------------------------------------------------------------------------
% This function is used to generate a control strategy during dry period
% Inputs: system physical properties, real-time system states, rainfall forecast data
% Outputs: setpoints for the inlet orifice
% RTCmodel 1: RBC
% RTCmodel 2: PFL-RBC
% RTCmodel 3: MPC
% ------------------------------------------------------------------------------------------------------------

% ------------------------------------------------------------------------------------------------------------
% Geting the inputs of system physical properties 
% ------------------------------------------------------------------------------------------------------------
Hmin = TankPara.Hmin;
Hmax = TankPara.Hmax;
As = TankPara.As;
A = TankPara.A;
constep = ControlPara.constep;
Tmin = ControlPara.Tmin;
RA = ControlPara.RA;
Td = ControlPara.Td; 
Tmax = ControlPara.Tmax;
Wuse = ControlPara.Wuse;
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

% ------------------------------------------------------------------------------------------------------------
% Geting the inputs of rainfall forecast data
% ------------------------------------------------------------------------------------------------------------
PP24 = PredictPara.PP24;
PP6 = PredictPara.PP6;

switch RTCmodel
    case 1 % RBC
        if td >= Td * 60 / constep
            Rule = 4;
        elseif td >= Tmin * 60 / constep
            Rule = 5;
        else
            Rule = 6;
        end
        
    case 2 % PFL-RBC
        if PP24 > 0
            Rule = 3;
        elseif td >= Tmax * 60 / constep 
            Rule = 4;
        elseif td >= Tmin * 60 / constep
            Rule = 5;
        else
            Rule = 6;
        end      
    case 3 % MPC
        if PP24 > 0
            Rule = 3;
        elseif td >= Tmax * 60 / constep 
            Rule = 4;
        elseif td >= Tmin * 60 / constep
            Rule = 5;
        else
            Rule = 6;
        end    
end

switch Rule
    case 3 % pre-storm release rule
        if PP6 > 0
            Hp = Function_runoffpredict(PP24, A, As, Hmax, Hmin, RA);
            if h > Hp + 0.0001 
                Ap1 = 1;
            else
                Ap1 = 0;
            end            
        
        else
            Ap1 = 0;            
        end

        if h > Hmin && Wuse > 0
            Ap2 = 1;
        else
            Ap2 = 0;
        end
        Ao = 0;        
        
    case 4 % empty rule
        if h > Hmin
            Ap1 = 1;
        else
            Ap1 = 0;
        end
        Ap2 = 0;
        Ao = 0;
        
    case 5 % reuse rule
        if h > Hmin && Wuse > 0
            Ap2 = 1;
        else
            Ap2 = 0;
        end
        Ap1 = 0;
        Ao = 0;
        
    case 6 % detention rule
        Ap1 = 0;
        Ap2 = 0;
        Ao = 0;
end 

end



