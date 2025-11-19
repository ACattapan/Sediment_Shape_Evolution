function mu=mu_of_L_func(L,mu_max,Ka)
%--------------------------------------------------------------------------
% THIS CODE COMPUTES THE RELATIVE MASS LOSS AS A FUNCTION OF TRAVEL
% DISTANCE.
%--------------------------------------------------------------------------
arguments
    L double
    mu_max double
    Ka double
end
%--------------------------------------------------------------------------
mu=mu_max-exp(-Ka.*L);