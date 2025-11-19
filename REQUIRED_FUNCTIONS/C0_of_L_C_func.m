function C0=C0_of_L_C_func(L,mu_max,C,C_max,a,Ka)
%--------------------------------------------------------------------------
% THIS CODE COMPUTES THE RELATIVE MASS LOSS AS A FUNCTION OF TRAVEL
% DISTANCE.
%--------------------------------------------------------------------------
arguments
    L double
    mu_max double
    C double
    C_max double
    a double
    Ka double
end
%--------------------------------------------------------------------------
mu=mu_of_L_func(L,mu_max,Ka);
C0=(C-C_max.*(1-(1-mu./mu_max).^a).^(1/a))./(1-(1-(1-mu./mu_max).^a).^(1/a));