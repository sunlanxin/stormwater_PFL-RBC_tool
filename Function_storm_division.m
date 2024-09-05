function [Pst, Pstnew, Ped, Pednew, Qst, Qstnew, Qed, Qednew, PQ_features] = Function_storm_division(P, Q, PQdiviPara, time)
% ------------------------------------------------------------------------------------------------------------
% Rainfall classification criteria: 
% minimum rainfall intensity of 1.2 mm/h (0.1mm/5min), 
% minimum cumulative rainfall of 10mm, 
% time interval of 2h
% ------------------------------------------------------------------------------------------------------------
PTmin = PQdiviPara.PTmin; % minimum cumulative rainfall, mm
PDmin = PQdiviPara.PDmin; % minimum duration, 5min
PImin = PQdiviPara.PImin; % minimum rainfall intensity, mm/5min
interval = PQdiviPara.interval; % time interval, h
deltaT = PQdiviPara.deltaT;% time period length, 5min
Qmin = PQdiviPara.Qmin; % minimum flow threshold, L/s
Pmaxmin = PQdiviPara.Pmaxmin; % minimum rain peak threshold, L/s

%% Rainfall series division
N = size(P,1);
TdP = zeros(N,1);
event = 1;
Pst = zeros(N,1);
Ped = zeros(N,1);
for n = 1 : N
    if n == 1
        if P(n,1) >= PImin
            Pst(event, 1) = n;
        end
    else
        if P(n, 1) >= PImin && P(n-1, 1) < PImin && ~Pst(event, 1)
            Pst(event, 1) = n;
        end
    end
    
    if P(n,1) < PImin  && n > 1
        TdP(n, 1) = TdP(n-1, 1) + 1;
        if n < N
            if TdP(n, 1) >= interval * 60 / deltaT && P(n+1, 1) >= PImin && Pst(event, 1) 
                Ped(event, 1) = n - TdP(n, 1);
                event = event + 1;
            end
        end
    else
        TdP(n,1) = 0;
    end
end
if Pst(event, 1) && ~Ped(event, 1)
    Ped(event, 1) = n - TdP(n, 1);
end
Pst = Pst(1:event, 1);
Ped = Ped(1:event, 1);

%% Flow series division in reference to rainfall
Qst = Pst;
Qed = Ped;
for i = 1:event
    a = Ped(i, 1);
    if i < event
        b = Qst(i+1, 1) - 1;
    else
        b = N;
    end
    
    for j = b:-1:a
        if sum(Q(j:b)) <= Qmin
            Qed(i, 1) = j;
        end
    end
end

%% Screening of effective rainfall-flow events and corresponding
Pstnew = zeros(event, 1);
Pednew = zeros(event, 1);
Qstnew = zeros(event, 1);
Qednew = zeros(event, 1);
eventnew = 0;
for e = 1 : event
   
    duration = Ped(e, 1) - Pst(e, 1);
    Total = sum(P(Pst(e, 1):Ped(e, 1)));
    Pmax = max(P(Pst(e, 1):Ped(e, 1)));
    if duration > PDmin/deltaT && Total > PTmin  && Pmax >= Pmaxmin
        eventnew = eventnew + 1;
        Pstnew(eventnew, 1) =  Pst(e, 1);
        Pednew(eventnew, 1) =  Ped(e, 1);
        
        % Characteristic values
        PDuration(eventnew, 1) = duration;
        PTotal(eventnew, 1) = Total;
        Ppeak(eventnew, 1) = max(P(Pst(e, 1):Ped(e, 1)));
        id1 = Pst(e, 1)-1+find(P(Pst(e, 1):Ped(e, 1)) == Ppeak(eventnew, 1),1);
        Peaktime(eventnew, 1) = strcat(string(time(id1,1)), "	", string(time(id1,2)));
        
        QTotal(eventnew, 1) = sum(Q(Qst(e, 1):Qed(e, 1))) * deltaT * 60 / 1000;
        Qpeak(eventnew, 1) = max(Q(Qst(e, 1):Qed(e, 1)));
        id2 = Qst(e, 1)-1+find(Q(Qst(e, 1):Qed(e, 1)) == Qpeak(eventnew, 1),1);
        QPeaktime(eventnew, 1) = strcat(string(time(id2,1)),"	",string(time(id2,2)));
        
        Qstnew(eventnew, 1) =  Qst(e, 1);
        Qednew(eventnew, 1) =  Qed(e, 1);
    end
end

Pstnew = Pstnew(1:eventnew, 1);
Pednew = Pednew(1:eventnew, 1);
Qstnew = Qstnew(1:eventnew, 1);
Qednew = Qednew(1:eventnew, 1);
PQ_features = [PDuration, PTotal, Ppeak, Peaktime, QTotal, Qpeak, QPeaktime];