a=[1];ka=[1];L=[0:0.001:8];C0=[0.8];%, 0.85 , 3, 0.1
mu_max=1;
C_max=1;
figure;hold on;
ax_GS=gca;
% xscale(ax_GS,"log")
pause('on')
for k=1:size(C0,2)
    for i=1:size(a,2)
        % if i==1
        %     LineColor=''
        for j=1:size(ka,2)
            mu=mu_of_L_func(L,mu_max,ka(j));
            C{i,j,k}=C_of_mu_func(mu,mu_max,C0(k),C_max,a(i));
            plot(L,C{i,j,k});
            % pause(5)
        end
    end
end
a=[6];ka=[0.1];C0=[0.8];
mu=mu_of_L_func(L,mu_max,ka);
C_plot=C_of_mu_func(mu,mu_max,C0,C_max,a);
plot(L,C_plot);