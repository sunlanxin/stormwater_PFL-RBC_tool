function PP = Function_perfectPP(hour, Pcs, tnow, deltat, te)
% ------------------------------------------------------------------------------------------------------------
% This function is used to calculate the idealized forecast rainfall within a specified time range
% tnow: start time number
% hour: time range, hours
% Pcs: monitored rainfall sequence
% deltat: time step of rainfall sequence, min
% te: length of rainfall sequence
% ------------------------------------------------------------------------------------------------------------
startTime =  tnow;
endTime = min(startTime + hour*60/deltat - 1, te);
PP = sum(Pcs(startTime:endTime));

end