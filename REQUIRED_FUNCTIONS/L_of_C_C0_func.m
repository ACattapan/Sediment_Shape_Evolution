function L=L_of_C_C0_func(C,mu_max,C0,C_max,a,Ka)
%--------------------------------------------------------------------------
% THIS CODE IMPORTS .tif IMAGES OF TRAVEL DISTANCES FOR EACH SAMPLING LOCATION
% AND SAVES
%--------------------------------------------------------------------------
arguments
    C       double
    mu_max  double
    C0      double
    C_max   double
    a       double
    Ka      double
end
%--------------------------------------------------------------------------
if Ka>0.0
    % for i=1:length(mu)
    %     L(i,1)=-1/Ka*log(1-mu(i));
    % end
    mu=mu_func(C,mu_max,C0,C_max,a);
    L=-1/Ka.*log(1-mu);
end
% L=Length_func(mu_func(C,mu_max,C0,C_max,a),K);