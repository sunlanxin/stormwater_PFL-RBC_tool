function PP = Function_actualPP(num, Psdata, hour)
% ------------------------------------------------------------------------------------------------------------
% This function is used to get the realized forecast rainfall from the given Psdata
% ------------------------------------------------------------------------------------------------------------
if hour == 3
    PP = str2num(Psdata(num, 2)); 
elseif hour == 6
    PP = str2num(Psdata(num, 3));
elseif hour == 24
    PP = str2num(Psdata(num, 4));
else
   error('Error. \n The entered foresight period is outside the data range.')
end
    
end