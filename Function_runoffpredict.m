function Hp = Function_runoffpredict(PP24, A, As, Hmax, Hmin, RA)
% ------------------------------------------------------------------------------------------------------------
% This function is used to calculate the water level for pre-storm release, based on the prediction the total inflow volume in the future 24 hours.
% A multi-level rainfall-runoff coefficient method is employed to make the prediction.
% ------------------------------------------------------------------------------------------------------------
if PP24 > 100
    RA_new = RA;
elseif  PP24 > 50
    RA_new = 0.75 * RA;
elseif  PP24 > 10
    RA_new = 0.5 * RA;
elseif  PP24 > 3
    RA_new = 0.25 * RA;
else
    RA_new = 0;
end

Wnext = (PP24 / 1000) * RA_new * A;
Hp = max(Hmax - Wnext / As, Hmin);

end