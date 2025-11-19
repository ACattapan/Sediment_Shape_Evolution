function mu=mu_func(C,mu_max,C0,C_max,a)
%--------------------------------------------------------------------------
% THIS CODE IMPORTS .tif IMAGES OF TRAVEL DISTANCES FOR EACH SAMPLING LOCATION
% AND SAVES
%--------------------------------------------------------------------------
arguments
    C double
    mu_max double
    C0 double
    C_max double
    a double
end
%--------------------------------------------------------------------------
for i=1:length(C)
    if and(C(i)>=C0,C(i)<=C_max)
        mu(i,1)=mu_max*(1-(1-((C(i)-C0)/(C_max-C0))^a)^(1/a));
    elseif C(i)<C0
        mu(i,1)=0;
    elseif C(i)>C_max
        mu(i,1)=1;
    end
end