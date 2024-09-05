function Qmax = Function_peakpredict_flc(PP24, Vmax, Hmin, Hp, h)
% ------------------------------------------------------------------------------------------------------------
% % This function is used to read the trained FLC controller and calculate the target flow rate
% input: 24-hour forecast rainfall depth, storage capacity
% output: FLC-calculated target flow
% ------------------------------------------------------------------------------------------------------------
h = max(h, Hmin);
if Hp >= h + 0.0001
    Qmax = 0;
else    
a2 = readfis('target flow');
Qmax = evalfis([PP24, Vmax],a2);   
end