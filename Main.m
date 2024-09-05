%% Statements
% ------------------------------------------------------------------------------------------------------------
% This .m file is the main function to implement the PFL-RBC approach, an RTC algorithm, 
% which incorporate the fuzzy logic control and the rule-based control.

% Demonstration case: a  storage tank with inlet orifice and outlet pump

% Control objectives: 
% i) water quantity (peak flow reduction, runoff volume reduction), 
% ii) water quality (detention time), and 
% iii) economy (stormwater reuse rate)

% System model: 
% i) TVGM-Urban urban hydrological model and
% ii) SWMM model

% Control model(RTCmodel=): 
% 1) Rule-Based Control (RBC), 
% 2) Predictive Fuzzy Logic and Rule-Based Control (PFL-RBC) 
% 3) Model Prediction Control (MPC)

% Created time: 17 October 2023 by Lanxin Sun
% ------------------------------------------------------------------------------------------------------------
clc;clear all;
%% Control model selection
RTCmodel = 1; 
%% Reading input data
path = 'F:\孙蓝心\成果\SCI论文3_2024.3~2024\返修1-WR\stormwater_PFL-RBC_tool\Data\';
%path = 'D:\PFL-RBC\Data\';
load(strcat(path, 'paras_TVGM.mat')); % TVGM-Urban urban hydrological model parameters
load(strcat(path, 'PPx.mat'));% forecast rainfall data
load(strcat(path, 'Pcs.mat'));% monitored rainfall data
load (strcat(path, 'time.mat'));% time series
swmminp = strcat(path, 'tankswmm-zynf_17.5.inp');% input files of SWMM model
%% Parameters and variables
% ------------------------------------------------------------------------------------------------------------
% Definition and initialtion of variables
% ------------------------------------------------------------------------------------------------------------
P = []; % monitored rainfall series, mm
PP6 = []; % forecasted rainfall depth for the next 6 hours, mm
PP24 = []; % forecasted rainfall depth for the next 24 hours, mm
h = []; % water level inside the storage tank, m
H1 = [];% water level before the orifice, m
It = []; % total inflow to the storage system, L/s
QIN = []; % inflow rate to the tank (controlled by the inlet orifice), L/s
QOUT = []; % outflow rate from the tank (controlled by the outlet pump1), L/s
Ot = []; % total outflow to the downstream, L/s
QOVER = []; % overflow rate, L/s
QUSE = []; % reuse flow rate (controlled by the outlet pump2), L/s
Hp = []; % the prepared water level for pre-storm release, m
orifice(1) = 0; % setpoint of the inlet orifice, 0~1
pump1(1) = 0; % setpoint of the outlet pump1, 0 or 1
pump2(1) = 0;% setpoint of the outlet pump2, 0 or 1
% ------------------------------------------------------------------------------------------------------------
% Definition of system model parameters
% ------------------------------------------------------------------------------------------------------------
% area, TVGM-Urban model - total catchment area, (m2)
% area_sub, TVGM-Urban model - subsurface components, [permeable ratio, LID ratio, impervious ratio]
% dt, TVGM-Urban model - simulation time step, (min)
% ------------------------------------------------------------------------------------------------------------
% Definition of control model parameters
% ------------------------------------------------------------------------------------------------------------
% Hmin, minimum water level of the storage tank, (m)
% Hmax, maximum water level of the storage tank, (m)
% As, bottom area of storage pool tank, (m2)
% d, orifice diameter, (m)
% Z0, orifice inlet drop, (m)
% Vmax, maximum storage capacity, (m3)
% constep, control step length, (min)
% Td, the fixed detention time  in RBC (h)
% Tmin, the minimum detention time before reuse in PFL-RBC/MPC, (h) 
% Tmax, the maximum detention time before empty in PFL-RBC/MPC, (h) 
% RA, the rainfall-runoff coefficient for calculating pre-storm release volume
% Pmin, wet/dry period discrimination: minimum rainfall threshold, (mm)
% Qmin, wet/dry period discrimination: minimum flow threshold, (L/s)
% Tuse1, start time of water reuse per day, hh:mm:ss
% Tuse2, end time of water reuse per day, hh:mm:ss
% usetime, length of reuse time per day, (min)
% Qmax, target flow in PFL-RBC, (L/s)
% Tc, the length of control horizon in MPC, (h) 
% Tp,  the length of prediction horizon in MPC, (h)
% ------------------------------------------------------------------------------------------------------------
% Assignment of system model parameters
% ------------------------------------------------------------------------------------------------------------
area = 4447.1;
area_sub = [0.33,0.5,0.17]; 
dt = 5; 
RunoffPara = struct('Id', 1, 'g1', paras_TVGM(1,1), 'g2', paras_TVGM(2,1), 'g3', paras_TVGM(3,1),...
    'g4', paras_TVGM(4,1), 'g5', paras_TVGM(5,1), 'g6', paras_TVGM(6,1), 'W0', paras_TVGM(7,1),...
    'Wp', paras_TVGM(8,1), 'Wlid', paras_TVGM(9,1), 'Dm', paras_TVGM(10,1), 'IMD', paras_TVGM(11,1), 'TP', paras_TVGM(12,1));
FlowPara = struct('Id', 2, 'n', paras_TVGM(13,1), 'k', paras_TVGM(14,1));
FlowPara2 = struct('Id', 2, 'n', paras_TVGM(15,1), 'k', paras_TVGM(16,1));
% ------------------------------------------------------------------------------------------------------------
% Assignment of control model parameters
% ------------------------------------------------------------------------------------------------------------
TankPara = struct('Hmin', 1, 'Hmax', 3.1, 'As', 8.75, 'd', 0.3, 'Z0', 0.2, 'A',area);
ControlPara = struct('constep', 5, 'Tmin', 19, 'Td', 72, 'Tmax', 168, 'RA', 0.3, 'Qmax', 0, 'Wuse', 1, 'Tc', 0.5, 'Tp', 24);
constep = ControlPara.constep;
Vmax = TankPara.As * (TankPara.Hmax - TankPara.Hmin);
Pmin = 0.2;
Qmin = 0.1;
Tuse1 = "10:00:00";
Tuse2 = "15:00:00";
usetime = 60;

%% Simmulation of the total inflow (It)
% ------------------------------------------------------------------------------------------------------------
% TVGM-Urban model
% ------------------------------------------------------------------------------------------------------------
P = Pcs;
N = size(P,1);
PI = P/dt;
It = TVGM_URBAN(P, PI, area/1000000, area_sub, RunoffPara, FlowPara, FlowPara2, dt/60) * 1000;
It(It<0.1) = 0;
% ------------------------------------------------------------------------------------------------------------
% Loose coupling of TVGM-Urban and SWMM model
% ------------------------------------------------------------------------------------------------------------
name(1:N,1) = "Qin";
datetime = time;
inp_Qin = [name, datetime, It];
TVGM_to_swmm(swmminp, inp_Qin);

%% Control model: simulation, decision-making, execution
% ------------------------------------------------------------------------------------------------------------
% Initial boundary
% ------------------------------------------------------------------------------------------------------------
td = 0;
h(1) = TankPara.Hmin;         
QIN(1) = 0;
QOUT(1) = 0;
Ot(1) = 0;
QOVER(1) = 0;
QUSE(1) = 0;
Hp(1) = 0;
orifice(1) = 0;
pump1(1) = 0;
pump2(1) = 0;
Aorifice = orifice(1);
Apump1 = pump1(1);
Apump2 = pump2(1);
Qmax = 0;
Qtarget = 0;
use = 0;
peakcal = 0;

% ------------------------------------------------------------------------------------------------------------
% Simulation start
% ------------------------------------------------------------------------------------------------------------
swmm = SWMM;
swmm.initialize(swmminp);
links = swmm.get_all(swmminp, swmm.LINK, swmm.NONE);
nodes = swmm.get_all(swmminp, swmm.NODE, swmm.NONE);
ts = 1; % number of simulation time steps
tc = 1; % number of control time steps
n = 1; % number of targete flow updating
k = 0; % SWMM iteration step
while (~swmm.is_over)
    k = k + 1;
    tnow(k,1) = swmm.run_step * 60;
    if fix(tnow(k,1)) / dt == ts 
        % read the current system states, which in practice are acquired through real-time monitoring of the system
        H1(ts) = swmm.get('Yb12', swmm.DEPTH, swmm.SI); 
        Qrunoff(ts) = swmm.get('Yb12', swmm.INFLOW, swmm.SI) * 1000;
        Qrunoff(ts) = min(Qrunoff(ts), It(ts+1));
        h(ts) = swmm.get('T1', swmm.DEPTH, swmm.SI);
        QIN(ts) = swmm.get('FM1', swmm.FLOW, swmm.SI) * 1000;
        QOVER(ts) = swmm.get('FM2', swmm.FLOW, swmm.SI) * 1000;
        QOUT(ts) = swmm.get('P1', swmm.FLOW, swmm.SI) * 1000;
        QUSE(ts) = swmm.get('P2', swmm.FLOW, swmm.SI) * 1000; 
        Ot(ts) = swmm.get('J23', swmm.FLOW, swmm.SI) * 1000;
        if fix(tnow(k,1)) / constep == tc 
        % ------------------------------------------------------------------------------------------------------------
        % A control time step
        % ------------------------------------------------------------------------------------------------------------
            timedate = strcat(string(time(ts,1)), " ", string(time(ts,2)));
            if contains(timedate, Tuse1) || contains(timedate, Tuse2) || (use > 0 && use < usetime/constep)
                ControlPara.Wuse = 1;
                use = use + 1;
            else
                ControlPara.Wuse = 0;
                use = 0;
            end
            % ------------------------------------------------------------------------------------------------------------
            % Save system state variables and control points
            % ------------------------------------------------------------------------------------------------------------
            StatePara = struct('It', It(ts+1), 'Qrunoff', Qrunoff(ts), 'H1', H1(ts), 'h', h(ts), 'Qin', QIN(ts), 'Qout', QOUT(ts), 'Qover', QOVER(ts), 'Quse', QUSE(ts), 'Q', Ot(ts));
            SetpointPara = struct('Aorifice', orifice(ts), 'Apump1', pump1(ts), 'Apump2', pump2(ts));
            
            fprintf("System states reading complete! Currently the %d control time step:", tc);
            % ------------------------------------------------------------------------------------------------------------
            % Forecast rainfall calculator
            % ------------------------------------------------------------------------------------------------------------
            PP = P(ts+1:min(ts+1+24*60/constep, N)); % rainfall process pediction for MPC
            QP = It(ts+1:min(ts+1+24*60/constep, N)); % inflow process pediction for MPC
            % realistic forecast
            % num = fix(ts*dt/60, 0) + 1;
            % PP6 = PPx(num, 3);
            % PP24 = PPx(num, 4);
            
            % idealized forecast
            PP24(ts) = Function_perfectPP(ControlPara.Tp, P, ts, dt, N);
            PP6(ts) = Function_perfectPP(ControlPara.Tc, P, ts, dt, N);
            
            PredictPara = struct('PP', PP, 'QP', QP, 'PP24', PP24(ts), 'PP6', PP6(ts));
            % ------------------------------------------------------------------------------------------------------------
            % Weather condition discrimination for triggering a control algorithm
            % ------------------------------------------------------------------------------------------------------------
            Hp(ts) = Function_runoffpredict(PP24(ts), area, TankPara.As, TankPara.Hmax, TankPara.Hmin, ControlPara.RA);
            if StatePara.It > Qmin
                % ------------------------------------------------------------------------------------------------------------
                % Wet period - predictive FLC
                % ------------------------------------------------------------------------------------------------------------
                td(ts) = 0;
                Qmax(ts) = Function_peakpredict_flc(PP24(ts), Vmax, TankPara.Hmin, Hp(ts), h(ts));
                if Qmax(ts) > Qtarget(n, 1)
                    Qtarget(n+1, 1) = Qmax(ts); 
                    Qtarget(n+1, 2) = ts;
                    ControlPara.Qmax = Qtarget(n+1, 1);
                    n = n + 1;
                    peakcal = 0;
                end
                [Aorifice, Apump1, Apump2] = control_model_wet(RTCmodel, StatePara, TankPara, ControlPara, PredictPara, SetpointPara, timedate);
                
            else
                % ------------------------------------------------------------------------------------------------------------
                % Dry period - predictive RBC
                % ------------------------------------------------------------------------------------------------------------
                if StatePara.h <= TankPara.Hmin + 0.0001
                    td(ts) = 0;
                elseif ts == 1
                    td(ts) = 1;
                else
                    td(ts) = td(ts-1) + 1;
                end
                
                if peakcal == 0
                    Qtarget(n+1, 1) = 0;
                    Qtarget(n+1, 2) = ts;
                    ControlPara.Qmax = Qtarget(n+1, 1);
                    n = n + 1;
                    peakcal = 1;
                end
                [Aorifice, Apump1, Apump2] = control_model_dry(RTCmodel, StatePara, TankPara, ControlPara, PredictPara, SetpointPara, td(ts), timedate);

            end
            QT(ts) = ControlPara.Qmax;
            fprintf("决策结果为:%f,%f,%f\n", Aorifice, Apump1, Apump2);
            % ------------------------------------------------------------------------------------------------------------
            % Adjust the setpoints in SWMM
            % ------------------------------------------------------------------------------------------------------------
            swmm.modify_setting('FM1', Aorifice);
            swmm.modify_setting('P1', Apump1);
            swmm.modify_setting('P2', Apump2);                  
            tc = tc + 1;
        end
        % ------------------------------------------------------------------------------------------------------------
        % Adjust the setpoints in SWMM
        % ------------------------------------------------------------------------------------------------------------
        % H(ts+1) = h;
        % QIN(ts+1) = Qin;
        % QOUT(ts+1) = Qout;
        % QOVER(ts+1) = Qover;
        % QUSE(ts+1) = Quse;
        % Ot(ts+1) = Q;
        swmm.modify_setting('FM1', Aorifice);
        swmm.modify_setting('P1', Apump1);
        swmm.modify_setting('P2', Apump2);
        orifice(ts+1) = Aorifice;
        pump1(ts+1) = Apump1;
        pump2(ts+1) = Apump2;
                      
        % ------------------------------------------------------------------------------------------------------------
        % State variables are corrected based on the water balance equations
        % ------------------------------------------------------------------------------------------------------------
        if orifice(ts+1)==1
            Qin_new(ts) = Qrunoff(ts);           
        else
            Qin_new(ts) = min(Function_valve_to_flow(H1(ts), TankPara.d, TankPara.Z0, orifice(ts+1)), Qrunoff(ts));
        end
        QOVER(ts) = min(QOVER(ts),Qin_new(ts));
        Ot(ts) = min(Qrunoff(ts)+QOUT(ts)+QOVER(ts), Ot(ts));      
        Ot_new(ts) = Qrunoff(ts)+QOUT(ts)+QOVER(ts)-Qin_new(ts);       
        ts = ts + 1;
    end
end
Ot_new(Ot_new < 0.1) = 0;
[errors,duration] = swmm.finish;

%% Performance evaluation
% ------------------------------------------------------------------------------------------------------------
% Rainfall event division
% ------------------------------------------------------------------------------------------------------------
% PTmin = 10, minimum cumulative rainfall, mm
% PDmin = 30, minimum rainfall duration, min
% PImin = 0.1, minimum rainfall intensity, mm/5min
% interval = 2, time interval, h
% deltaT = 5, time period length, 5min
% Qmin = 0.1, minimum flow threshold, L/s
PQdiviPara = struct('PTmin', 10, 'Pmaxmin', 1, 'PDmin', 30, 'PImin', 0.1, 'interval', 2, 'deltaT', 5, 'Qmin', 0.1);
[Pst, Pstnew, Ped, Pednew, Qst, Qstnew, Qed, Qednew, PQ_features] = Function_storm_division(P(2:end), It(2:end), PQdiviPara, time(2:end, :));
eventnew = length(Pstnew);
Ot2 = Ot_new - QOUT(1:N-2);
for e = 1 : eventnew
    a = Pstnew(e, 1); 
    b = Qednew(e, 1); 

    next = find(Pst == a)+1;
    if next >= size(Pst, 1)
        c = N-2;
    else
        c = Pst(next, 1)-1; 
    end
    
    if e < eventnew
        d = Pstnew(min(e+1,eventnew), 1); 
    else
        d = N-2;
    end     
    
    % ------------------------------------------------------------------------------------------------------------
    % Peak flow control of a single event
    % ------------------------------------------------------------------------------------------------------------
    Event_inpeak(e, 1) = max(It(a+1:b+1));
    % Event_outpeak(e, 1) = max(max(Ot(a:b)-QOUT(a:b)), 0);
    Event_outpeak(e, 1) = max(Ot2(a:b));
    Event_Rpeak(e, 1) = (Event_inpeak(e, 1) - Event_outpeak(e, 1)) / Event_inpeak(e, 1);
    
    % ------------------------------------------------------------------------------------------------------------
    % Runoff volume control of a single event
    % ------------------------------------------------------------------------------------------------------------
     Event_runoffin(e, 1) = (0.5*(It(a+1)+It(b+1))+sum(It(a+2:b)) * dt * 60 / 1000);
    Event_runoffout(e, 1) = (0.5*(Ot2(a)+Ot2(b))+sum(Ot2(a+1:b-1)) * dt * 60 / 1000);
    Event_Rrunoffin(e, 1) = 1 - Event_runoffin(e, 1) / (sum(P(a:b)) * area / 1000); %径流总量控制率-源头
    Event_Rrunoffout(e, 1) = 1 - Event_runoffout(e, 1) / (sum(P(a:b)) * area / 1000); %径流总量控制率-系统
    Event_Rrunoff(e, 1) = 1 -  Event_runoffout(e, 1) / Event_runoffin(e, 1);   
    % ------------------------------------------------------------------------------------------------------------
    % Detention control of a single event
    % ------------------------------------------------------------------------------------------------------------
    Event_Rtd(e, 1) = max(td(a:d));
    
    % ------------------------------------------------------------------------------------------------------------
    % Reuse control of a single event
    % ------------------------------------------------------------------------------------------------------------
    Event_reuse(e, 1) = (0.5*(QUSE(a)+QUSE(d))+sum(QUSE(a+1:d-1)) * dt * 60 / 1000);
          
end

% ------------------------------------------------------------------------------------------------------------
% Average peak flow reduction rate - Rpeak;
% Average runoff volume reduction rate - Rrunoff;
% Average detention time - Rtd;
% Average reuse ratio - Ruse
% ------------------------------------------------------------------------------------------------------------
Rpeak = sum(Event_Rpeak) / eventnew;
Rrunoffout = sum(Event_Rrunoff) / eventnew;
Rtd = sum(Event_Rtd) / eventnew;
RuseS =sum(QUSE * 300 / 1000) / sum(It * 300 / 1000);%系统雨水回用率，回用量/系统总入流


R_quantity = 0.5 * (Rpeak + Rrunoffout);
R_quanlity = Rtd / (60 / constep);
R_economy = RuseS;

result_event = [Event_Rpeak, Event_inpeak, Event_outpeak, Event_Rrunoff, Event_Rtd, Event_reuse];
result_total = [Rpeak, Rrunoffout,R_quantity, R_quanlity, R_economy, W_use];

