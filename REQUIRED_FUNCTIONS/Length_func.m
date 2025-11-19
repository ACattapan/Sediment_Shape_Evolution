function L=Length_func(mu,K)
%--------------------------------------------------------------------------
% THIS CODE IMPORTS .tif IMAGES OF TRAVEL DISTANCES FOR EACH SAMPLING LOCATION
% AND SAVES
%--------------------------------------------------------------------------
arguments
    mu double
    K double
end
%--------------------------------------------------------------------------
if K>0.0
    % for i=1:length(mu)
    %     L(i,1)=-1/K*log(1-mu(i));
    % end
    L=-1/K.*log(1-mu);
end
% L=Length_func(mu_func(C,mu_max,C0,C_max,a),K);