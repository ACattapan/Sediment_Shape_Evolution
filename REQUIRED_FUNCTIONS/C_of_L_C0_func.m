function C=C_of_L_C0_func(L,mu_max,C0,C_max,a,Ka)
%--------------------------------------------------------------------------
% THIS CODE COMPUTES THE RELATIVE MASS LOSS AS A FUNCTION OF TRAVEL
% DISTANCE.
%--------------------------------------------------------------------------
arguments
    L double
    mu_max double
    C0 double
    C_max double
    a double
    Ka double
end
%--------------------------------------------------------------------------
mu=mu_of_L_func(L,mu_max,Ka);
C=C0+(C_max-C0).*(1-(1-mu./mu_max).^a).^(1/a);