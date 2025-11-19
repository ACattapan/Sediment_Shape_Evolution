function C=C_of_mu_func(mu,mu_max,C0,C_max,a)
%--------------------------------------------------------------------------
% THIS CODE COMPUTES THE RELATIVE MASS LOSS AS A FUNCTION OF TRAVEL
% DISTANCE.
%--------------------------------------------------------------------------
arguments
    mu double
    mu_max double
    C0 double
    C_max double
    a double
end
%--------------------------------------------------------------------------
C=C0+(C_max-C0).*(1-(1-mu./mu_max).^a).^(1/a);